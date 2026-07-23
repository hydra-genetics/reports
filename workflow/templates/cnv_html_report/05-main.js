const getTextDimensions = function (text, fontSize) {
  const canvas =
    getTextDimensions.canvas ||
    (getTextDimensions.canvas = document.createElement("canvas"));
  const context = canvas.getContext("2d");
  context.font = `${fontSize} sans-serif`;
  const metrics = context.measureText(text);
  return [
    metrics.width,
    metrics.actualBoundingBoxAscent + metrics.actualBoundingBoxDescent,
  ];
};

// Modal handling start
function hideModalOnClick(evt) {
  const dimensions = this.getBoundingClientRect();
  if (
    evt.clientX < dimensions.left ||
    evt.clientX > dimensions.right ||
    evt.clientY < dimensions.top ||
    evt.clientY > dimensions.bottom
  ) {
    this.close();
  }
}

for (const dialog of document.querySelectorAll("dialog .close")) {
  dialog.addEventListener("click", (e) => {
    e.preventDefault();
    e.currentTarget.parentNode.close();
  });
}

const messageModal = document.getElementById("chromosome-view-dialog");
messageModal.addEventListener("click", hideModalOnClick);

function setModalMessage(msg, className) {
  let message = document.createElement("p");
  const messageText = document.createTextNode(msg);

  let icon = document.createElement("i");

  if (className === "error") {
    icon.className = "bi-exclamation-circle-fill";
  } else if (className === "warning") {
    icon.className = "bi-exclamation-circle-fill";
  } else if (className === "info") {
    icon.className = "bi-info-circle-fill";
  }

  message.appendChild(icon);
  message.appendChild(messageText);

  messageModal.className = className ? className : "";
  messageModal.firstChild?.remove();
  messageModal.prepend(message);
}

const helpModal = document.getElementById("help-modal");
helpModal.addEventListener("click", hideModalOnClick);

function showHelp() {
  helpModal.showModal();
}
// Modal handling end

const cnFromRatio = function (ratio, refPloidy) {
  if (refPloidy === undefined) {
    refPloidy = 2;
  }
  return refPloidy * 2 ** ratio;
};

// Sort callers to put Jumble first
const sortedCallers = [...cnvData[0].callers].sort((a, b) => {
  if (a.name.toLowerCase() === 'jumble') return -1;
  if (b.name.toLowerCase() === 'jumble') return 1;
  return a.label.localeCompare(b.label);
});

d3.select("#dataset-picker")
  .selectAll("div")
  .data(sortedCallers)
  .join("div")
  .call((e) => {
    e.append("input")
      .attr("type", "radio")
      .property("checked", (_, i) => i === 0)
      .attr("value", (d) => {
        // Find the original index in cnvData[0].callers
        return cnvData[0].callers.findIndex(c => c.name === d.name);
      })
      .attr("id", (d) => `dataset-${d.name}`)
      .attr("name", "dataset");
    return e;
  })
  .call((e) => {
    e.append("label")
      .attr("for", (d) => `dataset-${d.name}`)
      .text((d) => d.label);
    return e;
  });

const initialCallerIndex = cnvData[0].callers.findIndex(c => c.name.toLowerCase() === 'jumble');
const chromosomePlot = new ChromosomePlot({
  element: document.querySelector("#chromosome-view"),
  data: cnvData[0],
  tc: originalTc,
  caller: initialCallerIndex !== -1 ? initialCallerIndex : 0,
  widePlotWidth: widePlotWidth,
});

const genomePlot = new GenomePlot({
  element: document.querySelector("#genome-view"),
  data: cnvData,
  tc: originalTc,
  caller: initialCallerIndex !== -1 ? initialCallerIndex : 0,
  widePlotWidth: widePlotWidth,
});

const resultsTable = new ResultsTable(d3.select("#cnv-table"), {
  data: cnvData,
  tc: originalTc,
  caller: initialCallerIndex !== -1 ? initialCallerIndex : 0,
  filter: d3.select("#table-filter-toggle").node().checked,
});

// Show binned data warning if any part of the dataset was downsampled by the backend
if (cnvData[0].is_baf_binned || cnvData[0].is_log2_binned) {
  d3.select(".binned-data-warning").classed("hidden", false);
  
  // Update tooltip to be specific about what was binned
  let warningTitle = "Data has been downsampled for performance. 'Show all points' only shows the downsampled data.";
  if (cnvData[0].is_baf_binned && cnvData[0].is_log2_binned) {
    warningTitle = "Both BAF and Log2 data have been downsampled. 'Show all points' only shows the binned points.";
  } else if (cnvData[0].is_baf_binned) {
    warningTitle = "BAF data has been downsampled. 'Show all points' only shows the binned BAF points.";
  } else if (cnvData[0].is_log2_binned) {
    const callers = cnvData[0].binned_callers.join(", ");
    warningTitle = `Log2 data for ${callers} has been downsampled. 'Show all points' only shows the binned Log2 points.`;
  }
  d3.select(".binned-data-warning").attr("title", warningTitle);
}

chromosomePlot.addEventListener("zoom", (e) => {
  d3.selectAll(".data-range-warning").classed(
    "hidden",
    !e.detail.dataOutsideRange
  );
});

chromosomePlot.addEventListener("max-zoom-reached", () => {
  setModalMessage(
    "Trying to zoom in too far. " +
    `Current lower limit is ${chromosomePlot.minZoomRange} bp.`,
    "error"
  );
  messageModal.showModal();
});

genomePlot.addEventListener("chromosome-change", (e) => {
  chromosomePlot.setData(
    cnvData[e.detail.chromosome],
    e.detail.start,
    e.detail.end
  );
});

genomePlot.addEventListener("chromosome-zoom", (e) => {
  chromosomePlot.setData(null, e.detail.start, e.detail.end);
});

resultsTable.addEventListener("zoom-to-region", (e) => {
  genomePlot.selectChromosome(
    e.detail.chromosome,
    e.detail.start,
    e.detail.start + e.detail.length
  );
});

d3.select("#table-filter-toggle").on("change", (event) => {
  resultsTable.filter = event.target.checked;
});

d3.select("#chromosome-fit-to-data").on("change", (e) => {
  chromosomePlot.fitToData = e.target.checked;
});

d3.select("#chromosome-show-all-datapoints").on("change", (e) => {
  chromosomePlot.showAllData = e.target.checked;
  genomePlot.showAllData = e.target.checked;
});

const cancerGeneColoringToggle = d3.select("#chromosome-cancer-gene-coloring");

cancerGeneColoringToggle.on("change", (e) => {
  chromosomePlot.cancerGeneColoring = e.target.checked;
  const legend = d3.select("#cancer-gene-legend");
  legend.classed("hidden", !e.target.checked);
  if (e.target.checked) {
    updateCancerGeneLegend(cnvData);
  }
});

function updateCancerGeneLegend(data) {
  const legend = d3.select("#cancer-gene-legend");
  legend.html(""); // Clear existing

  const roles = new Map();

  // Collect unique roles and their colors from all chromosomes
  data.forEach((chrom) => {
    if (chrom.annotations) {
      chrom.annotations.forEach((ann) => {
        if (ann.role && ann.color && !roles.has(ann.role)) {
          roles.set(ann.role, ann.color);
        }
      });
    }
  });

  if (roles.size === 0) {
    legend.classed("hidden", true);
    return;
  }

  // Sort roles alphabetically
  const sortedRoles = Array.from(roles.keys()).sort();

  // Initialize active roles if empty
  if (chromosomePlot.activeCancerGeneRoles.size === 0) {
    chromosomePlot.activeCancerGeneRoles = new Set(sortedRoles);
  }

  sortedRoles.forEach((role) => {
    const color = roles.get(role);
    const item = legend.append("div").attr("class", "legend-item");
    const label = item.append("label").attr("class", "legend-item-label");

    label
      .append("input")
      .attr("type", "checkbox")
      .attr("class", "legend-checkbox")
      .property("checked", chromosomePlot.activeCancerGeneRoles.has(role))
      .on("change", (e) => {
        const activeRoles = new Set(chromosomePlot.activeCancerGeneRoles);
        if (e.target.checked) {
          activeRoles.add(role);
        } else {
          activeRoles.delete(role);
        }
        chromosomePlot.activeCancerGeneRoles = activeRoles;
      });

    label
      .append("div")
      .attr("class", "legend-color")
      .style("background-color", color);
    label.append("span").attr("class", "legend-label").text(role);
  });
}

const baselineOffsetSlider = d3.select("#chromosome-baseline-offset");
const currentBaselineOffset = d3.select("#current-baseline-offset");
const baselineOffsetReset = d3.select("#reset-baseline-offset");

// Converts the anchor's raw log2 value into the baseline-offset slider value
// (psi = 2 * 2^baselineOffset is the assumed neutral copy number in the pure
// -tumor frame - see #toAbsoluteCopyNumber). The anchor's raw log2 is the
// *observed*, TC-diluted value, so recovering the pure-tumor-frame psi that
// makes this segment resolve to exactly 2 absolute copies requires knowing
// TC: solving adjCopies(anchorLog2) = 2 in #toAbsoluteCopyNumber's model for
// psi (not assuming tc=1) gives psi = (K - 2*(1-tc)) / tc, where
// K = 2^(1-anchorLog2) is tc-independent. Using tc=1 here (i.e. baselineOffset
// = -anchorLog2) would only be correct when "Simulate purity" stays off -
// since TC can change independently after a baseline is set, the offset must
// be solved jointly with whichever TC value is currently effective.
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
// see reapplyBaselineFromAnchor below.
function anchorLog2FromBaselineOffset(baselineOffset, tc) {
  const psi = 2 * 2 ** baselineOffset;
  const K = tc * psi + 2 * (1 - tc);
  return 1 - Math.log2(K);
}

// The raw log2 value the currently active baseline offset treats as exactly
// 2 absolute copies (see baselineOffsetFromAnchor/anchorLog2FromBaselineOffset
// above). baselineOffset alone isn't TC-independent once "Simulate purity" is
// active - the same anchor point needs a different baselineOffset at a
// different TC to keep resolving to exactly 2 copies - so this is tracked
// separately and used to keep the two in sync whenever TC changes (or
// "Simulate purity" is toggled, which switches the TC that's actually applied
// between 1 and the real value). Starts at 0 to match the default
// baselineOffset=0, which is itself TC-independent (see derivation).
let lastAnchorLog2 = 0;

// simulatePurity/currentTc are declared further below, but this is only
// ever called from within event handlers, which all run after the entire
// script (including those declarations) has finished executing.
function effectiveTc() {
  return simulatePurity.node().checked ? parseFloat(currentTc.node().value) : 1;
}

// Applies a baseline offset value to the slider/display/views only - does
// not touch lastAnchorLog2. Callers decide separately whether this offset
// should update the tracked anchor (a manual/estimated change should) or not
// (reapplyBaselineFromAnchor recomputes an offset *from* the existing
// anchor, so it must leave it untouched).
function setBaselineOffsetUI(dy) {
  const strdy = dy.toLocaleString("en-US", { minimumFractionDigits: 2 });
  baselineOffsetSlider.node().value = dy;
  currentBaselineOffset.node().value = strdy;
  currentBaselineOffset.node().classList.remove("invalid");
  currentBaselineOffset.node().title = "";
  baselineOffsetReset.property("disabled", dy === 0);
  chromosomePlot.setBaselineOffset(dy);
  genomePlot.setBaselineOffset(dy);
  resultsTable.setBaselineOffset(dy);
}

// Recomputes the baseline offset from the tracked anchor at whatever TC is
// currently effective, keeping that anchor pinned at exactly 2 absolute
// copies as TC changes. No-op if no anchor has ever been established.
function reapplyBaselineFromAnchor() {
  if (lastAnchorLog2 === null) return;
  const minDy = parseFloat(baselineOffsetSlider.node().min);
  const maxDy = parseFloat(baselineOffsetSlider.node().max);
  const dy = Math.min(
    maxDy,
    Math.max(minDy, baselineOffsetFromAnchor(lastAnchorLog2, effectiveTc()))
  );
  setBaselineOffsetUI(dy);
}

baselineOffsetSlider.on("change", () => {
  currentBaselineOffset.node().dispatchEvent(new Event("change"));
});

baselineOffsetSlider.on("input", (e) => {
  const dy = parseFloat(e.target.value);
  const strdy = dy.toLocaleString("en-US", { minimumFractionDigits: 2 });
  currentBaselineOffset.node().value = strdy;
});

currentBaselineOffset.on("change", (e) => {
  const minDy = e.target.min ? parseFloat(e.target.min) : -2.0;
  const maxDy = e.target.max ? parseFloat(e.target.max) : 2.0;
  const dy = parseFloat(e.target.value);

  if (dy < minDy || dy > maxDy) {
    e.target.classList.add("invalid");
    e.target.title = `Value outside the valid range [${minDy}, ${maxDy}]`;
    console.error(
      `baseline offset outside the valid range [${minDy}, ${maxDy}]`
    );
    return;
  }

  setBaselineOffsetUI(dy);
  // A manual edit defines a fresh anchor at whichever raw log2 this offset
  // currently corresponds to, at the currently-effective TC.
  lastAnchorLog2 = anchorLog2FromBaselineOffset(dy, effectiveTc());
});

baselineOffsetReset.on("click", () => {
  baselineOffsetSlider.node().value = 0;
  baselineOffsetReset.property("disabled", true);
  currentBaselineOffset.node().value = "0.00";
  baselineOffsetSlider.node().dispatchEvent(new Event("change"));
});

const simulatePurity = d3.select("#simulate-purity");
const tcAdjustSlider = d3.select("#tc-adjuster");
const currentTc = d3.select("#current-tc");
const tcAdjustReset = d3.select("#reset-tc");

const roundSegmentsToInteger = d3.select("#round-segments-integer");
const viewModeInputs = d3.selectAll("input[name=view-mode]");

function updateRoundSegmentsEnablement() {
  const enable = simulatePurity.node().checked;
  roundSegmentsToInteger.property("disabled", !enable);
  if (!enable) {
    roundSegmentsToInteger.property("checked", false);
    roundSegmentsToInteger.node().dispatchEvent(new Event("change"));
  }
}

function applyViewMode(mode) {
  viewModeInputs.property("checked", function () {
    return this.value === mode;
  });
  chromosomePlot.setViewMode(mode);
  genomePlot.setViewMode(mode);
  updateRoundSegmentsEnablement();
}

simulatePurity.on("change", (e) => {
  const checked = e.target.checked;
  tcAdjustSlider.property("disabled", !checked);
  currentTc.property("disabled", !checked);
  if (!checked) {
    // resultsTable always applies TC (unlike the plots, which assume 100%
    // purity while "Simulate purity" is off) — reset the (now-disabled)
    // TC input back to the default so the table doesn't keep using a
    // stale, no-longer-adjustable TC the user tried earlier.
    tcAdjustSlider.node().value = originalTc;
    currentTc.node().value = originalTc;
    tcAdjustReset.property("disabled", true);
  }
  currentTc.node().dispatchEvent(new Event("change"));
  chromosomePlot.setSimulatePurity(checked);
  genomePlot.setSimulatePurity(checked);
  if (checked) {
    applyViewMode("copyNumber");
  } else {
    applyViewMode("log2");
  }
});

roundSegmentsToInteger.on("change", (e) => {
  const checked = e.target.checked;
  chromosomePlot.setRoundSegmentsToInteger(checked);
  genomePlot.setRoundSegmentsToInteger(checked);
  resultsTable.setRoundSegmentsToInteger(checked);
});

viewModeInputs.on("change", (e) => {
  applyViewMode(e.target.value);
});

tcAdjustSlider.on("change", () => {
  currentTc.node().dispatchEvent(new Event("change"));
});

tcAdjustSlider.on("input", (e) => {
  const dv = parseFloat(e.target.value);
  const strdv = dv.toLocaleString("en-US", { minimumFractionDigits: 2 });
  currentTc.node().value = strdv;
});

currentTc.on("change", (e) => {
  const minTc = e.target.min ? parseFloat(e.target.min) : 0;
  const maxTc = e.target.max ? parseFloat(e.target.max) : 1;
  const tc = parseFloat(e.target.value);

  if (tc < minTc || tc > maxTc) {
    e.target.classList.add("invalid");
    e.target.title = `Value outside the valid range [${minTc}, ${maxTc}]`;
    console.error(
      `tumor cell content outside the valid range [${minTc}, ${maxTc}]`
    );
    return;
  }

  e.target.classList.remove("invalid");
  e.target.title = "";

  tcAdjustReset.property("disabled", true);
  if (tc != originalTc) {
    tcAdjustReset.property("disabled", false);
  }

  const strtc = tc.toLocaleString("en-US", { minimumFractionDigits: 2 });
  tcAdjustSlider.node().value = tc;
  currentTc.node().value = strtc;
  chromosomePlot.setTc(tc);
  genomePlot.setTc(tc);
  resultsTable.setTc(tc);

  // Keep whichever anchor is currently tracked pinned at exactly 2 copies
  // under the new TC. This also fires whenever "Simulate purity" is
  // toggled, since that handler unconditionally dispatches "change" on
  // currentTc to apply the effective-TC switch between 1 and the real value.
  reapplyBaselineFromAnchor();
});

tcAdjustReset.on("click", () => {
  tcAdjustSlider.node().value = originalTc;
  tcAdjustReset.property("disabled", true);
  currentTc.node().value = originalTc;
  tcAdjustSlider.node().dispatchEvent(new Event("change"));
});

const Y_ZOOM_STEP = 0.8; // zoom-in shrink factor; zoom-out is 1/Y_ZOOM_STEP

const yZoomIn = d3.select("#y-zoom-in");
const yZoomOut = d3.select("#y-zoom-out");
const yZoomReset = d3.select("#reset-y-zoom");
let currentYZoomFactor = 1;

function applyYZoom(factor) {
  // Clamp here too (setYZoom clamps internally) so the tracked factor never
  // drifts past what's actually rendered - otherwise repeated zoom-in clicks
  // past the floor would require extra zoom-out clicks before anything
  // visibly changes again. The max depends on the active view mode (copy
  // number vs log2 have very different useful zoom-out ranges - see
  // maxYZoomFactorFor); both plots share the same view mode.
  currentYZoomFactor = Math.min(
    maxYZoomFactorFor(chromosomePlot.viewMode === "copyNumber"),
    Math.max(MIN_Y_ZOOM_FACTOR, factor)
  );
  yZoomReset.property("disabled", currentYZoomFactor === 1);
  chromosomePlot.setYZoom(currentYZoomFactor);
  genomePlot.setYZoom(currentYZoomFactor);
}

yZoomIn.on("click", () => applyYZoom(currentYZoomFactor * Y_ZOOM_STEP));
yZoomOut.on("click", () => applyYZoom(currentYZoomFactor / Y_ZOOM_STEP));
yZoomReset.on("click", () => applyYZoom(1));

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomePlot.activeCaller = parseInt(e.target.value);
  genomePlot.activeCaller = parseInt(e.target.value);
  resultsTable.activeCaller = parseInt(e.target.value);
});

d3.select("#chromosome-equal-distance").on("change", (e) => {
  chromosomePlot.equalDistance = e.target.checked;
});

const genes = new Map();

function addGeneLocation(gene, chromosome, start, end) {
  if (!genes.has(gene)) {
    genes.set(gene, []);
  }
  const locations = genes.get(gene);
  // Avoid exact duplicates
  const exists = locations.some(
    (l) => l.chromosome === chromosome && l.start === start && l.end === end
  );
  if (!exists) {
    locations.push({ chromosome, start, end });
  }
}

// Add comprehensive gene index first (from RefSeq)
if (cnvData[0].gene_search_index) {
  Object.entries(cnvData[0].gene_search_index).forEach(([gene, info]) => {
    addGeneLocation(gene, info.chrom, info.start, info.end);
  });
}

// Add annotation genes (overwrites if duplicates, but usually same or better resolution on annotations)
cnvData.forEach((chromData) => {
  chromData.annotations.forEach((anno) => {
    addGeneLocation(anno.name, chromData.chromosome, anno.start, anno.end);
  });
});

const geneList = d3.select("#gene-list");
Array.from(genes.keys())
  .sort()
  .forEach((gene) => {
    geneList.append("option").attr("value", gene);
  });

let lastSearchTerm = "";
let searchIndex = 0;

d3.select("#gene-search").on("change", (e) => {
  const geneName = e.target.value;
  const geneLocations = genes.get(geneName);
  const errorIcon = d3.select("#gene-search-error");

  if (geneLocations && geneLocations.length > 0) {
    errorIcon.style("display", "none");

    // Group locations by chromosome
    const locsByChrom = new Map();
    geneLocations.forEach(loc => {
        if (!locsByChrom.has(loc.chromosome)) {
            locsByChrom.set(loc.chromosome, []);
        }
        locsByChrom.get(loc.chromosome).push(loc);
    });
    const chromosomes = Array.from(locsByChrom.keys()); // e.g. ["1", "17"]

    // Cycle through chromosomes if same search
    if (geneName === lastSearchTerm) {
      searchIndex = (searchIndex + 1) % chromosomes.length;
    } else {
      searchIndex = 0;
      lastSearchTerm = geneName;
    }

    const targetChrom = chromosomes[searchIndex];
    const targetRegions = locsByChrom.get(targetChrom);
    const firstRegion = targetRegions[0]; // For determining start/end zoom if desired? 
    // Usually we just switch to chromosome. The plot handles zoom? 
    // No, standard `selectChromosome` resets zoom to full chromosome unless specified.
    // The previous code didn't zoom to gene, it just switched chromosome and highlighted.
    
    // We need to add "name" to each region for the plotter to use it as tooltip/label
    // Pass ALL locations so highlights persist if user switches chromosome manually
    const regionsWithName = geneLocations.map(r => ({...r, name: geneName}));

    const chromIndex = cnvData.findIndex(
      (d) => d.chromosome === targetChrom
    );
    if (chromIndex !== -1) {
      // Switch chromosome if needed
      if (genomePlot.selectedChromosomeIndex !== chromIndex) {
        genomePlot.selectChromosome(cnvData[chromIndex].chromosome);
      }

      // Wait for chromosome switch to complete, then highlight ALL regions
      setTimeout(() => {
        chromosomePlot.highlightRegions(regionsWithName);
      }, 50);

      // Keep input populated but blur to allow "Enter" to trigger change again if desired?
      // "change" only triggers on blur or enter committed change. 
      // To allow re-triggering with Enter, we might need to keep focus or handle "keyup".
      // But standard implementation: 
      // If user types "TP53" [Enter], highlights 1. 
      // If user hits [Enter] again, does "change" fire? No, only if value changed.
      // So we clear value? If we clear value, next search is "fresh".
      // User wants to see multiple hits.
      
      // Adaptation: Don't clear value. Let user press Enter again? 
      // D3 "change" event only fires if value changed and lost focus or enter.
      // If value is same, "change" might not fire on second Enter.
      // We might need to listen to 'keydown' for Enter.
      
      // For now, let's keep the value. 
      // e.target.value = ""; // Don't clear
       e.target.blur(); // Blur to show we are done. 
       // If they click and hit enter again, it might not trigger change if value is same.
       // Let's add specific key listener for cycling.
    }
  } else {
    if (geneName !== "") {
      errorIcon.style("display", "inline-block");
    } else {
      errorIcon.style("display", "none");
    }
  }
});

// Add keydown listener to force cycle on Enter even if value hasn't changed
d3.select("#gene-search").on("keydown", (e) => {
    if (e.key === "Enter") {
        e.preventDefault();
        // Manually trigger the change logic if value is same as last search
        if (e.target.value === lastSearchTerm) {
             d3.select("#gene-search").dispatch("change");
        } else {
            // standard behavior will trigger change on Enter anyway
             e.target.blur(); // Trigger change
        }
    }
});

// Reconcile "Absolute copy number"'s enabled/checked state with whatever
// "Simulate purity" was rendered as, in case a site configures
// default_absolute_copy_number=true without default_simulate_purity=true
// (that combination is invalid — absolute copy number requires purity
// simulation — and this forces the checkbox back to unchecked+disabled
// before anything below can act on its template-rendered "checked" state).
updateRoundSegmentsEnablement();

// Apply config-driven default-checked state for controls whose "checked"
// attribute alone doesn't apply its side effects (enabling other controls,
// informing the plot objects) — dispatch "change" as if the user had just
// clicked it, going through the exact same handlers/cascade above.
// Order matters: simulate-purity first, since round-segments-integer's own
// enablement (and default-checked state, if misconfigured without purity
// simulation) depends on that cascade having already run.
[
  "#simulate-purity",
  "#round-segments-integer",
  "#chromosome-equal-distance",
  "#chromosome-show-all-datapoints",
].forEach((selector) => {
  const node = d3.select(selector).node();
  if (node.checked) {
    node.dispatchEvent(new Event("change"));
  }
});
