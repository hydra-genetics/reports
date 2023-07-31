const getTextDimensions = function (text, fontSize) {
  let div = document.createElement("div");

  div.innerText = text;
  div.style.position = "absolute";
  div.style.float = "left";
  div.style.fontSize = fontSize;
  div.style.whiteSpace = "nowrap";
  div.style.visibility = "hidden";

  document.body.append(div);
  let width = div.clientWidth;
  let height = div.clientHeight;
  div.remove();

  return [width, height];
};

// Populate dataset picker
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

// Chromosome plot
const chromosomeView = new ChromosomePlot({
  element: document.querySelector("#chromosome-view"),
  data: cnvData[0],
});

// Genome plot
const genomeView = new GenomePlot({
  element: document.querySelector("#genome-view"),
  data: cnvData,
});

// CNV table
const resultsTable = new ResultsTable(d3.select("#cnv-table"), {
  data: cnvData,
  filter: d3.select("#table-filter-toggle").node().checked,
});

// Event listeners
chromosomeView.addEventListener("zoom", (e) => {
  d3.selectAll(".data-range-warning").classed(
    "hidden",
    !e.detail.dataOutsideRange
  );
});

genomeView.addEventListener("chromosome-change", (e) => {
  chromosomeView.data = cnvData[e.detail.chromosome];
});

resultsTable.addEventListener("zoom-to-region", (e) => {
  genomeView.selectChromosome(e.detail.chromosome);
  chromosomeView.zoomTo(e.detail.start, e.detail.start + e.detail.length);
});

d3.select("#table-filter-toggle").on("change", (event) => {
  resultsTable.filter = event.target.checked;
});

d3.select("#chromosome-fit-to-data").on("change", (e) => {
  chromosomeView.fitToData = e.target.checked;
});

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomeView.activeCaller = parseInt(e.target.value);
  genomeView.activeCaller = parseInt(e.target.value);
  resultsTable.activeCaller = parseInt(e.target.value);
});
