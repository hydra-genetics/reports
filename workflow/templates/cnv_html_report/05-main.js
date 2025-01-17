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

const cnFromRatio = function (ratio, refPloidy) {
  if (refPloidy === undefined) {
    refPloidy = 2;
  }
  return refPloidy * 2 ** ratio;
};

d3.select("#dataset-picker")
  .selectAll("div")
  .data(cnvData[0].callers)
  .join("div")
  .call((e) => {
    e.append("input")
      .attr("type", "radio")
      .property("checked", (_, i) => i === 0)
      .attr("value", (_, i) => i)
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

const chromosomePlot = new ChromosomePlot({
  element: document.querySelector("#chromosome-view"),
  data: cnvData[0],
});

const genomePlot = new GenomePlot({
  element: document.querySelector("#genome-view"),
  data: cnvData,
});

const resultsTable = new ResultsTable(d3.select("#cnv-table"), {
  data: cnvData,
  filter: d3.select("#table-filter-toggle").node().checked,
});

const messageModal = document.querySelector("dialog");
messageModal.addEventListener("click", (e) => {
  const modalDimensions = messageModal.getBoundingClientRect();
  if (
    e.clientX < modalDimensions.left ||
    e.clientX > modalDimensions.right ||
    e.clientY < modalDimensions.top ||
    e.clientY > modalDimensions.bottom
  ) {
    messageModal.close();
  }
});

document.querySelector("dialog button.close").addEventListener("click", (e) => {
  e.currentTarget.parentNode.close();
});

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
});

const baselineOffsetSlider = d3.select("#chromosome-baseline-offset");
const currentBaselineOffset = d3.select("#current-baseline-offset");
const baselineOffsetReset = d3.select("#reset-baseline-offset");

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

  e.target.classList.remove("invalid");
  e.target.title = "";

  baselineOffsetReset.property("disabled", true);
  if (dy != 0) {
    baselineOffsetReset.property("disabled", false);
  }

  const strdy = dy.toLocaleString("en-US", { minimumFractionDigits: 2 });
  baselineOffsetSlider.node().value = dy;
  currentBaselineOffset.node().value = strdy;
  chromosomePlot.setBaselineOffset(dy);
  genomePlot.setBaselineOffset(dy);
});

baselineOffsetReset.on("click", () => {
  baselineOffsetSlider.node().value = 0;
  baselineOffsetReset.property("disabled", true);
  currentBaselineOffset.node().value = "0.00";
  baselineOffsetSlider.node().dispatchEvent(new Event("change"));
});

const tcAdjustSlider = d3.select("#tc-adjuster");
const currentTc = d3.select("#current-tc");
const tcAdjustReset = d3.select("#reset-tc");

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
  // chromosomePlot.setTc(tc);
  // genomePlot.setTc(tc);
});

tcAdjustReset.on("click", () => {
  tcAdjustSlider.node().value = originalTc;
  tcAdjustReset.property("disabled", true);
  currentTc.node().value = originalTc;
  tcAdjustSlider.node().dispatchEvent(new Event("change"));
});

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomePlot.activeCaller = parseInt(e.target.value);
  genomePlot.activeCaller = parseInt(e.target.value);
  resultsTable.activeCaller = parseInt(e.target.value);
});
