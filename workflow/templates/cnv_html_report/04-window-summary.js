const MAX_POINTS = 1000;
const MIN_COPY_NUMBER = 0.001;
const MIN_LOG2_RATIO = Math.log2(MIN_COPY_NUMBER / 2);
const CN_VIEW_Y_MAX = 5;
const MIN_Y_ZOOM_FACTOR = 0.1; // prevents a degenerate zero-height Y domain
// Separate max zoom-out factors per view mode - a single shared max would
// force an awkward compromise, since "reasonable" zoom-out differs hugely
// between the two: copy number needs to reach real amplifications (which can
// exceed 100 copies), while log2 ratio's meaningful range is only a handful
// of units wide, and the same large factor would make it mostly blank space.
// 40x the default [0, CN_VIEW_Y_MAX] range gives a max of 200 copies at full
// zoom-out (with the floor-clamp in zoomYDomain keeping the bottom at 0).
const MAX_Y_ZOOM_FACTOR_COPY_NUMBER = 40;
// 5x the default [-2, 2] log2 range gives a max of roughly +/-10 - still
// generous (2^10 = 1024-fold), but far short of copy number's 40x.
const MAX_Y_ZOOM_FACTOR_LOG2 = 5;

function maxYZoomFactorFor(isCopyNumberView) {
  return isCopyNumberView ? MAX_Y_ZOOM_FACTOR_COPY_NUMBER : MAX_Y_ZOOM_FACTOR_LOG2;
}

// Shared axis tick-label helpers: pick tick POSITIONS so their (baseline
// -relabeled) LABEL values are always nice integers, rather than fixing
// positions at nice round numbers and letting labels drift to arbitrary
// decimals once baselineOffset != 0. `lo`/`hi` are already in LABEL space
// (i.e. post-relabeling); callers invert each returned label back to a
// position themselves (label - baselineOffset for log2, label / 2^baselineOffset
// for copy number).
function integerTickStep(lo, hi, targetCount) {
  const rawStep = (hi - lo) / Math.max(targetCount, 1);
  return Math.max(1, Math.round(rawStep));
}

function integerLabelsInRange(lo, hi, targetCount) {
  const step = integerTickStep(lo, hi, targetCount);
  const start = Math.ceil(lo / step) * step;
  const labels = [];
  for (let l = start; l <= hi; l += step) labels.push(l);
  return labels;
}

// Applies a Y-zoom factor around [yMin, yMax]'s center. In copy-number view,
// zooming out symmetrically would push the bottom below 0 - a copy number
// can never actually be negative, so instead of moving the "0" row up into
// the plot (wasting half the zoom-out budget on meaningless space), the
// bottom is clamped to 0 and the would-be-negative amount is redirected
// upward instead, so the bottom of the visible range always stays at 0.
// Named zoomYDomain (not applyYZoom) to avoid colliding with 05-main.js's
// distinct, UI-level applyYZoom(factor) - all report JS files concatenate
// into one script with no module namespacing, so a same-named top-level
// function declared later silently overwrites an earlier one.
function zoomYDomain(yMin, yMax, yZoomFactor, isCopyNumberView) {
  const mid = (yMin + yMax) / 2;
  const halfSpan = ((yMax - yMin) / 2) * yZoomFactor;
  let newYMin = mid - halfSpan;
  let newYMax = mid + halfSpan;
  if (isCopyNumberView && newYMin < 0) {
    newYMax += -newYMin;
    newYMin = 0;
  }
  return [newYMin, newYMax];
}

function* generateWindowSlices(points, scale, posAttr, windowSize = 5) {
  let offset = scale.domain()[0];
  let currentWindow = [];
  for (let p of points) {
    // Compute position from start/end if posAttr is missing
    const pos = p[posAttr] !== undefined ? p[posAttr] : (p.start + p.end) / 2;
    if (pos < offset) continue;
    if (offset > scale.domain()[1]) break;

    if (pos < offset + windowSize) {
      currentWindow.push(p);
    } else {
      while (pos >= offset + windowSize) {
        yield { window: currentWindow, start: offset };
        offset += windowSize;
        currentWindow = [];
      }
      currentWindow = [p];
    }
  }
  if (currentWindow.length > 0) {
    yield { window: currentWindow, start: offset };
  }
}

function summariseWindow(points, windowStart, windowSize, valAttr, minValue = undefined, offset = 0) {
  if (!points || points.length === 0) return null;
  let hasOutliers = false;
  let filtered = points;
  if (minValue !== undefined) {
    // Compare against a threshold shifted once, rather than adding `offset`
    // back onto every point — the latter is a subtract-then-add round trip
    // (the caller already subtracted `offset` into x[valAttr]) that is not
    // guaranteed to be bit-exact for arbitrary offset values, which could
    // silently let floor-clamped points escape the outlier filter.
    const adjustedMinValue = minValue - offset;
    hasOutliers = points.some((x) => (x[valAttr] ?? x.mean ?? 0) <= adjustedMinValue);
    filtered = points.filter((x) => (x[valAttr] ?? x.mean ?? 0) > adjustedMinValue);
  }
  if (filtered.length === 0) return null;

  let sum = filtered.reduce((a, b) => a + (b[valAttr] ?? b.mean ?? 0), 0);
  let mean = sum / filtered.length;

  // Calculate pooled variance: Between-group variance + Within-group variance
  // Between-group: Variance of the means
  let betweenVar = filtered.reduce((a, b) => a + Math.pow((b[valAttr] ?? b.mean ?? 0) - mean, 2), 0) / filtered.length;
  
  // Within-group: Average of the variances (sd^2) of the input points (if available)
  // If input is raw, sd is undefined/0, so withinVar contributes 0.
  let withinVar = filtered.reduce((a, b) => a + Math.pow(b.sd ?? 0, 2), 0) / filtered.length;

  let sd = Math.sqrt(betweenVar + withinVar);
  
  return {
    start: windowStart,
    end: windowStart + windowSize,
    pos: windowStart + windowSize / 2,
    baf: mean,
    mean: mean,
    sd: sd,
    hasOutliers: hasOutliers
  };
}

function slidingPixelWindowBAF(points, scale, posAttr = "pos", pixelWindowSize = 5, force = false) {
  const [d0, d1] = scale.domain();
  // Compute position from start/end if posAttr is missing
  points = points.filter(p => {
    const pos = p[posAttr] !== undefined ? p[posAttr] : (p.start + p.end) / 2;
    const endPos = p.end !== undefined ? p.end : pos;
    return endPos >= d0 && (p.start !== undefined ? p.start : pos) <= d1;
  });
  if (points.length === 0) return [];

  // Sort by position to ensure correct binning
  points.sort((a, b) => {
    const posA = a[posAttr] !== undefined ? a[posAttr] : (a.start + a.end) / 2;
    const posB = b[posAttr] !== undefined ? b[posAttr] : (b.start + b.end) / 2;
    return posA - posB;
  });

  if (points[0]?.baf instanceof Array) {
    points = points.flatMap(p => p.baf.map(v => ({ ...p, baf: v })));
  }
  
  let windowSize = Math.max(1, Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]));
  if (!force && (windowSize < 4 || points.length <= MAX_POINTS)) return points;

  let reducedPoints = [];
  for (let slice of generateWindowSlices(points, scale, posAttr, windowSize)) {
    let { window, start } = slice;
    if (window.length === 0) continue;

    let hist = Array(5).fill(0);
    window.forEach(p => {
      let b = p.baf ?? p.mean;
      if (b !== undefined) hist[Math.min(4, Math.floor(b * 5))] += 1;
    });

    if (hist.indexOf(Math.max(...hist)) !== 2) {
      let hi = window.filter(p => (p.baf ?? p.mean) >= 0.5);
      let lo = window.filter(p => (p.baf ?? p.mean) < 0.5);
      let sHi = summariseWindow(hi, start, windowSize, "baf");
      let sLo = summariseWindow(lo, start, windowSize, "baf");
      if (sHi) reducedPoints.push(sHi);
      if (sLo) reducedPoints.push(sLo);
    } else {
      let s = summariseWindow(window, start, windowSize, "baf");
      if (s) reducedPoints.push(s);
    }
  }
  return reducedPoints;
}



function slidingPixelWindow(
  points,
  scale,
  posAttr,
  valAttr,
  offset = 0,
  pixelWindowSize = 5,
  force = false,
  minValue = MIN_LOG2_RATIO
) {
  const [d0, d1] = scale.domain();
  points = points.filter(p => (p.end ?? p[posAttr]) >= d0 && (p.start ?? p[posAttr]) < d1);
  if (points.length === 0) return [];

  let windowSize = Math.max(1, Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]));
  if (!force && (windowSize < 4 || points.length <= MAX_POINTS)) {
    return points;
  }

  let reducedPoints = [];

  for (let slice of generateWindowSlices(points, scale, posAttr, windowSize)) {
    let window = slice.window;
    let windowStart = slice.start;
    if (window.length === 0) {
      continue;
    }
    let s = summariseWindow(
      window,
      windowStart,
      windowSize,
      valAttr,
      minValue,
      offset
    );
    if (s) reducedPoints.push(s);
  }

  return reducedPoints;
}

