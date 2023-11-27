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

d3.select("#reports-version")
  .call((p) =>
    p
      .select("a")
      .attr(
        "href",
        `https://github.com/hydra-genetics/reports/tree/v${VERSION}`
      )
  )
  .select("span")
  .text(VERSION);

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
    icon.className = "fa-solid fa-circle-exclamation";
  } else if (className === "warning") {
    icon.className = "fa-solid fa-triangle-exclamation";
  } else if (className === "info") {
    icon.className = "fa-solid fa-circle-info";
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

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomePlot.activeCaller = parseInt(e.target.value);
  genomePlot.activeCaller = parseInt(e.target.value);
  resultsTable.activeCaller = parseInt(e.target.value);
});
