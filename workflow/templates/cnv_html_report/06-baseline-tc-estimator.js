// Auto-suggest the baseline offset (the CN=2 anchor) and tumor content (TC),
// replacing the old PureCN-ploidy-driven "Adjust to ploidy" button.
//
// Heuristic (see docs/cnv_report.md and the design plan for the full
// derivation): find the lowest-log2, sufficiently large, BAF-balanced
// (non-LOH) segment genome-wide for the active caller and treat it as the
// CN=2 anchor - this sidesteps ploidy estimation entirely. Then, among
// segments that are clearly below that anchor and BAF-skewed (a real
// single-allele loss), solve for the TC that would resolve the *observed*
// (purity-diluted) BAF skew back to the theoretical pure-tumor extreme (0 or
// 1), reusing the same purity-dilution model as transformBAF/
// #toAbsoluteCopyNumber elsewhere in the report.

const MIN_SNPS_PER_SEGMENT = 10; // minimum BAF points to trust a segment's balance/skew call
const MIN_SNPS_PER_BIN = 3; // minimum BAF points for a position-bin to count
const N_POSITION_BINS = 10; // bins per segment for the position-weighted balance check
const EXTREME_BAF_MARGIN = 0.05; // drop raw baf < this or > 1-this before aggregating
const BALANCED_BAF_THRESHOLD = 0.1; // median-of-bin-medians |baf-0.5| below this => balanced
const MIN_SEGMENT_LENGTH = 5000000; // "relatively large" - fixed, configurable, default 5Mb
const CLUSTER_CAUTION_FRACTION = 0.5; // fraction of the TC-implied one-copy-state gap to use as
// clustering tolerance (see #findBaselineAnchor and findDeletionSegment) - configurable
const MIN_CLUSTER_TOTAL_LENGTH = 20000000; // min combined length to trust an anchor cluster

function median(values) {
  if (values.length === 0) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 !== 0
    ? sorted[mid]
    : (sorted[mid - 1] + sorted[mid]) / 2;
}

function bafPointPosition(point) {
  return point.pos !== undefined && point.pos !== null
    ? point.pos
    : Math.floor((point.start + point.end) / 2);
}

// Position-binned balance/skew statistic for a segment: bins matched BAF
// points by genomic position first, then takes the median of per-bin
// medians, so a small SNP-dense sub-region can't dominate a segment that's
// actually balanced (or skewed) across most of its genomic extent. Returns
// null when there isn't enough reliable data to trust a call either way.
function segmentBalance(segment, bafPointsForChrom) {
  const matched = bafPointsForChrom.filter((p) => {
    const pos = bafPointPosition(p);
    return pos >= segment.start && pos <= segment.end;
  });

  const filtered = matched.filter(
    (p) => p.baf >= EXTREME_BAF_MARGIN && p.baf <= 1 - EXTREME_BAF_MARGIN
  );
  if (filtered.length < MIN_SNPS_PER_SEGMENT) return null;

  const binWidth = (segment.end - segment.start) / N_POSITION_BINS;
  const bins = Array.from({ length: N_POSITION_BINS }, () => []);
  filtered.forEach((p) => {
    const pos = bafPointPosition(p);
    let idx = Math.floor((pos - segment.start) / binWidth);
    if (idx < 0) idx = 0;
    if (idx >= N_POSITION_BINS) idx = N_POSITION_BINS - 1;
    bins[idx].push(p);
  });

  const binMedians = bins
    .filter((bin) => bin.length >= MIN_SNPS_PER_BIN)
    .map((bin) => median(bin.map((p) => Math.abs(p.baf - 0.5))));
  if (binMedians.length === 0) return null;

  // balance = median-of-bin-medians of |baf-0.5|. Note this is NOT the same
  // as folding the pooled BAF values to their minor-allele side and taking
  // one plain median afterwards - without haplotype phasing, which allele
  // (ref/alt) is "lost" at a given SNP is arbitrary, so raw BAF values in a
  // real LOH/deletion segment split into two mirror-image clusters around
  // 0.5 (e.g. one near 0.3, one near 0.7); pooling them unfolded before
  // taking a median can land anywhere between the two modes and hide real
  // skew. Folding each point first (|baf-0.5|), as done here, avoids that.
  return { balance: median(binMedians) };
}

// Sex chromosomes are excluded from both the anchor and deletion-candidate
// pools: chrX's BAF pattern is confounded by X-inactivation (a normal,
// non-tumor-specific phenomenon that can look like allelic skew regardless
// of any real copy-number event), and chrX/Y copy number itself depends on
// the patient's sex rather than following the same autosomal 2-copies
// -unless-altered model everything else here assumes.
function isAutosome(chromosome) {
  return chromosome !== "chrX" && chromosome !== "chrY" && chromosome !== "X" && chromosome !== "Y";
}

// Every segment for the given caller, genome-wide (autosomes only - see
// isAutosome), with its length and BAF balance statistic attached. Segments
// where the balance call is unreliable (segmentBalance returned null) are
// dropped entirely.
//
// Known limitation: this reliability check runs per raw segment, BEFORE
// clusterByLog2Proximity groups same-event neighbors together. A caller
// whose segmentation happens to fragment part of a real event into pieces
// too small to individually clear MIN_SNPS_PER_SEGMENT will lose that part
// of the event entirely, even though the same underlying BAF data (shared
// across callers) might have been dense enough if pooled across the whole
// event first. This can make two callers' depth-based TC estimates for
// what's really the same deletion diverge, if one caller's segmentation is
// coarser (keeping the event whole) than another's (splitting it up).
// Moving this check to run on each cluster's pooled points, after
// clustering, would fix this, but is a real restructuring - not done here.
function collectCandidates(cnvData, callerIndex) {
  const candidates = [];
  cnvData.forEach((chromData) => {
    if (!isAutosome(chromData.chromosome)) return;
    const caller = chromData.callers[callerIndex];
    if (!caller || !caller.segments) return;
    caller.segments.forEach((segment) => {
      const stats = segmentBalance(segment, chromData.baf);
      if (stats === null) return;
      candidates.push({
        ...segment,
        length: segment.end - segment.start,
        balance: stats.balance,
      });
    });
  });
  return candidates;
}

function weightedMedianLog2(segments) {
  const sorted = [...segments].sort((a, b) => a.log2 - b.log2);
  const totalLength = sorted.reduce((sum, s) => sum + s.length, 0);
  let cumulative = 0;
  for (const s of sorted) {
    cumulative += s.length;
    if (cumulative >= totalLength / 2) return s.log2;
  }
  return sorted[sorted.length - 1].log2;
}

// The log2 gap to the next real copy-number state (CN=3) at the currently
// -assumed TC, scaled down by CLUSTER_CAUTION_FRACTION to stay well clear of
// that boundary - see findBaselineAnchor's comment for the full derivation.
// Shared between anchor clustering (how far a same-state neighbor can sit
// above the seed) and findDeletionSegment (how far a same-state,
// LOH-skewed neighbor can sit above the anchor).
function clusterTolerance(currentTc) {
  const gapToNextState = Math.log2(1 + currentTc / 2);
  return CLUSTER_CAUTION_FRACTION * gapToNextState;
}

// Find the CN=2 baseline anchor: the lowest-log2, balanced, sufficiently
// large segment, grown into a cluster of same-state neighbors using a
// TC-derived (not arbitrary) tolerance. Starting from the lowest untried
// balanced/large segment as seed, grows a cluster upward only (never
// downward - the seed is by construction the lowest untried candidate);
// every member must already be a balanced+large survivor, never merely
// close in log2. If a seed's cluster doesn't reach MIN_CLUSTER_TOTAL_LENGTH,
// that seed is discarded (this is what correctly skips an isolated low
// outlier rather than pulling the anchor down to it) and the next-lowest
// untried survivor is tried instead. Returns the anchor's log2 value, or
// null if no balanced/large segment exists at all.
function findBaselineAnchor(cnvData, callerIndex, currentTc) {
  const large = collectCandidates(cnvData, callerIndex).filter(
    (s) => s.length >= MIN_SEGMENT_LENGTH
  );
  const balanced = large.filter((s) => s.balance <= BALANCED_BAF_THRESHOLD);
  if (balanced.length === 0) return null;

  balanced.sort((a, b) => a.log2 - b.log2);
  const tolerance = clusterTolerance(currentTc);

  for (const seed of balanced) {
    const cluster = balanced.filter(
      (s) => s.log2 >= seed.log2 && s.log2 <= seed.log2 + tolerance
    );
    const totalLength = cluster.reduce((sum, s) => sum + s.length, 0);
    if (totalLength >= MIN_CLUSTER_TOTAL_LENGTH) {
      return weightedMedianLog2(cluster);
    }
  }

  // Last resort: no cluster anywhere reached the reliability bar - fall back
  // to the single largest balanced/large segment rather than fail outright.
  const largestSurvivor = balanced.reduce((best, s) =>
    s.length > best.length ? s : best
  );
  return largestSurvivor.log2;
}

// Greedily groups pre-sorted-by-log2 segments into clusters: a new segment
// joins the current cluster if it's within `tolerance` of the cluster's
// current LAST (i.e. highest-so-far) member, otherwise it starts a new
// cluster. This is what correctly merges same-event segments that
// segmentation split apart (e.g. two adjacent, near-identical-log2 segments
// separated only by a small intervening outlier segment) into one larger,
// more reliable combined candidate, rather than comparing each in isolation.
function clusterByLog2Proximity(segments, tolerance) {
  if (segments.length === 0) return [];
  const sorted = [...segments].sort((a, b) => a.log2 - b.log2);
  const clusters = [[sorted[0]]];
  for (let i = 1; i < sorted.length; i++) {
    const current = sorted[i];
    const currentCluster = clusters[clusters.length - 1];
    const previous = currentCluster[currentCluster.length - 1];
    if (current.log2 - previous.log2 <= tolerance) {
      currentCluster.push(current);
    } else {
      clusters.push([current]);
    }
  }
  return clusters;
}

// Among segments no higher than the anchor's own cluster range (raw log2
// space, scale-invariant - doesn't need TC to already be known beyond what
// clusterTolerance already uses) and BAF-skewed (unlike the balanced anchor
// candidates), group same-event neighbors together (see
// clusterByLog2Proximity) and pick the largest resulting cluster - preferring
// a genuine NET LOSS cluster (clearly below the anchor) over a SAME-LEVEL
// (copy-neutral LOH) cluster whenever both exist, regardless of size: a
// depth-confirmed net loss is a more reliable signal than BAF skew alone,
// which has no independent depth cross-check at the anchor's own level (see
// estimateTcFromNetLoss vs estimateTcFromDeletion).
//
// Deliberately has NO lower bound on the candidate pool itself: a segment
// doesn't need to sit below the anchor to be a valid single-allele-loss
// candidate. A segment at roughly the SAME log2 as the anchor but with
// skewed BAF is a copy-neutral LOH (e.g. 2+0 - one allele lost, the other
// duplicated to compensate, so total copy number matches the background) -
// just as valid a source for the BAF-skew TC formula as a genuine
// lower-log2 deletion (e.g. 1+0), since both are fully monoallelic in the
// pure tumor. Segments CLEARLY ABOVE the anchor are excluded because a real
// gain more plausibly retains both alleles in some ratio (e.g. 3+1) rather
// than being fully monoallelic, which the zero-BAF-in-pure-tumor assumption
// doesn't hold for. Any segment length is eligible here - these events can
// be focal/small, unlike the anchor which needs to be large.
//
// Returns { cluster: { log2, balance, length }, isNetLoss } | null, where
// log2/balance are length-weighted across the cluster's members.
function findDeletionSegment(cnvData, callerIndex, anchorLog2, currentTc) {
  const tolerance = clusterTolerance(currentTc);
  const upperBound = anchorLog2 + tolerance;
  const lowerBoundForSameLevel = anchorLog2 - tolerance;

  const all = collectCandidates(cnvData, callerIndex);
  const skewed = all.filter(
    (s) => s.log2 <= upperBound && s.balance > BALANCED_BAF_THRESHOLD
  );
  if (skewed.length === 0) return null;

  const clusters = clusterByLog2Proximity(skewed, tolerance).map((members) => {
    const length = members.reduce((sum, s) => sum + s.length, 0);
    return {
      members,
      length,
      log2: weightedMedianLog2(members),
      balance: members.reduce((sum, s) => sum + s.balance * s.length, 0) / length,
    };
  });

  const netLoss = clusters.filter((c) => c.log2 < lowerBoundForSameLevel);
  const sameLevel = clusters.filter((c) => c.log2 >= lowerBoundForSameLevel);
  const pool = netLoss.length > 0 ? netLoss : sameLevel;

  const chosen = pool.reduce((best, c) => (c.length > best.length ? c : best));
  return { cluster: chosen, isNetLoss: netLoss.length > 0 };
}

// Solve transformBAF's purity-dilution model (observed = tc*pureTumorBAF +
// (1-tc)*0.5) for TC, assuming the LOH is fully resolved in the pure tumor
// (pureTumorBAF = 0 or 1): tc = 1 - 2*minorAlleleFraction. Reuses
// segmentBalance's already-correctly-folded "balance" statistic rather than
// a separately pooled BAF median - folding each point to |baf-0.5| first
// and taking the median afterwards is order-preserving, so
// minorAlleleFraction = 0.5 - balance exactly; a plain median of pooled,
// unfolded BAF values would not be (see segmentBalance's comment). Used for
// SAME-LEVEL (copy-neutral LOH) clusters, where depth is uninformative by
// construction (see estimateTcFromNetLoss).
function estimateTcFromDeletion(cluster) {
  const minorAlleleFraction = 0.5 - cluster.balance;
  const tc = 1 - 2 * minorAlleleFraction;
  return Math.min(1, Math.max(0, tc));
}

// Solves #toAbsoluteCopyNumber's model directly for TC, assuming the
// cluster is a genuine NET LOSS (exactly 1 absolute copy) relative to the
// anchor: adjCopies(anchorLog2)=2 and adjCopies(cluster.log2)=1 hold
// simultaneously for the same TC. Since K = tc*psi + 2*(1-tc) = 2^(1-anchorLog2)
// is fixed by the anchor equation alone (tc-independent - see
// baselineOffsetFromAnchor/clusterTolerance's derivation), substituting into
// the second equation and solving directly for tc gives a closed form with
// no need to know psi/baselineOffset first:
//   tc = 2 - 2^(1 + cluster.log2 - anchorLog2)
// Only valid when the cluster sits clearly BELOW the anchor - AT the
// anchor's own level this is degenerate (cluster.log2≈anchorLog2 gives
// tc≈0, meaningless), which is exactly why same-level LOH clusters use
// estimateTcFromDeletion (BAF-based) instead. Depth is generally a more
// reliable signal than a single cluster's BAF here (better SNR), which is
// also why net-loss clusters are preferred over same-level ones in
// findDeletionSegment whenever both exist.
function estimateTcFromNetLoss(cluster, anchorLog2) {
  const tc = 2 - 2 ** (1 + cluster.log2 - anchorLog2);
  return Math.min(1, Math.max(0, tc));
}

// Converts the anchor's raw log2 value into the baseline-offset slider value
// (psi = 2 * 2^baselineOffset is the assumed neutral copy number in the pure
// -tumor frame - see #toAbsoluteCopyNumber). The anchor's raw log2 is the
// *observed*, TC-diluted value, so recovering the pure-tumor-frame psi that
// makes this segment resolve to exactly 2 absolute copies requires knowing
// TC: solving adjCopies(anchorLog2) = 2 in #toAbsoluteCopyNumber's model for
// psi (not assuming tc=1) gives psi = (K - 2*(1-tc)) / tc, where
// K = 2^(1-anchorLog2) is tc-independent (the same K used to derive the
// clustering tolerance above). Using tc=1 here (i.e. baselineOffset =
// -anchorLog2) would only be correct when "Simulate purity" stays off -
// since this estimator also sets/enables TC, the offset must be solved
// jointly with whichever TC value will actually end up active.
function baselineOffsetFromAnchor(anchorLog2, tc) {
  const K = 2 ** (1 - anchorLog2);
  const psi = (K - 2 * (1 - tc)) / tc;
  return Math.log2(psi / 2);
}

// Exact inverse of baselineOffsetFromAnchor: given a baseline offset that is
// (or is about to be) active at a given TC, what raw log2 value does it
// currently treat as exactly 2 absolute copies? Used to keep a
// baseline/anchor pair self-consistent whenever TC changes after the fact
// (either via the TC slider, or via toggling "Simulate purity", which
// switches the TC that's actually applied between 1 and the real value) -
// see 05-main.js's reapplyBaselineFromAnchor.
function anchorLog2FromBaselineOffset(baselineOffset, tc) {
  const psi = 2 * 2 ** baselineOffset;
  const K = tc * psi + 2 * (1 - tc);
  return 1 - Math.log2(K);
}

// Returns { baselineOffset, tc, anchorLog2 } | null. tc is itself null when
// a baseline anchor was found but no qualifying deletion segment was - the
// baseline can still be applied, but no TC suggestion is possible. Returns
// null overall only when no anchor could be found at all.
function estimateBaselineAndTc(cnvData, callerIndex, currentTc) {
  const anchorLog2 = findBaselineAnchor(cnvData, callerIndex, currentTc);
  if (anchorLog2 === null) return null;

  const found = findDeletionSegment(cnvData, callerIndex, anchorLog2, currentTc);
  const tc =
    found === null
      ? null
      : found.isNetLoss
        ? estimateTcFromNetLoss(found.cluster, anchorLog2)
        : estimateTcFromDeletion(found.cluster);

  // Solve the baseline offset against whichever TC will actually end up
  // active once this result is applied: the newly estimated TC if one was
  // found (since that's what gets set), otherwise the TC that's already
  // active (since TC is left unchanged in that case) - so the anchor
  // segment resolves to exactly 2 absolute copies in the final state.
  const baselineOffset = baselineOffsetFromAnchor(anchorLog2, tc !== null ? tc : currentTc);

  return { baselineOffset, tc, anchorLog2 };
}
