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

const getActiveCaller = () => {
  return +d3.select("input[name=dataset]:checked").node().value;
};

const plotGenomeView = () => {
  const totalLength = d3.sum(cnvData.map((d) => d.length));
  const height = 400;
  const width = 800;
  const margin = {
    top: 10,
    right: 30,
    bottom: 60,
    left: 60,
    between: 20,
  };
  const panelWidths = cnvData.map(
    (d) => ((width - margin.left - margin.right) * d.length) / totalLength
  );
  const panelHeight =
    (height - margin.top - margin.bottom - margin.between) / 2;

  const xScales = cnvData.map((d, i) =>
    d3.scaleLinear().domain([0, d.length]).range([0, panelWidths[i]])
  );

  const yScaleRange = 2;
  const yScale = d3
    .scaleLinear()
    .domain([-yScaleRange, yScaleRange])
    .range([panelHeight, 0]);

  const yScaleVAF = d3.scaleLinear().domain([0, 1]).range([panelHeight, 0]);

  const yAxis = (g) => g.call(d3.axisLeft(yScale).ticks(5));

  const yAxisVAF = (g) => g.call(d3.axisLeft(yScaleVAF).ticks(5));

  const svg = d3
    .select("#genome-view")
    .attr("preserveAspectRatio", "xMinYMin meet")
    .attr("viewBox", [0, 0, width, height])
    .attr("style", "max-width: 100%; max-height: 500px; height: auto;");

  const plotArea = svg
    .append("g")
    .attr("transform", `translate(${margin.left}, ${margin.top})`);

  const lrArea = plotArea.append("g").attr("class", "genome-view-area");
  const vafArea = plotArea
    .append("g")
    .attr("class", "genome-view-area")
    .attr("transform", `translate(0,${panelHeight + margin.between})`);

  const addPanels = function (g) {
    const panels = g
      .selectAll(".chromosome-panel")
      .data(cnvData)
      .join("g")
      .attr("data-index", (d, i) => i)
      .attr("class", "chromosome-panel")
      .attr(
        "transform",
        (d, i) =>
          `translate(${i === 0 ? 0 : d3.sum(panelWidths.slice(0, i))}, 0)`
      );

    // Panel backgrounds
    panels
      .append("rect")
      .attr("class", "bg-rect")
      .attr("width", (d, i) => panelWidths[i])
      .attr("height", panelHeight)
      .attr("fill", "#FFF")
      .attr("stroke", "#333");

    return panels;
  };

  const lrPanels = addPanels(lrArea);
  const vafPanels = addPanels(vafArea);

  const lrGrid = lrPanels
    .append("g")
    .attr("class", "grid")
    .attr("data-index", (d, i) => i);
  const vafGrid = vafPanels
    .append("g")
    .attr("class", "grid")
    .attr("data-index", (d, i) => i);

  const ratioPanels = lrPanels
    .append("g")
    .attr("class", "regions")
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("data-index", (d, i) => i);

  const segmentPanels = lrPanels
    .append("g")
    .attr("class", "segments")
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("data-index", (d, i) => i);

  const plotRegions = () => {
    ratioPanels
      .selectAll(".point")
      // Only plot every fifth point for performance
      .data(
        (d) =>
          d.callers[getActiveCaller()].ratios.filter((r, i) => i % 5 === 0),
        (d) => d.start
      )
      .join(
        (enter) =>
          enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d, i, g) =>
              xScales[g[i].parentNode.dataset.index](d.start)
            )
            .attr("cy", yScale(-yScaleRange - 0.2))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0.3)
            .call((enter) =>
              enter.transition().attr("cy", (d) => yScale(d.log2))
            ),
        (update) =>
          update.call((update) =>
            update
              .transition()
              .attr("cx", (d, i, g) =>
                xScales[g[i].parentNode.dataset.index](d.start)
              )
              .attr("cy", (d) => yScale(d.log2))
          ),
        (exit) =>
          exit
            .transition()
            .attr("cy", yScale(yScaleRange + 0.2))
            .remove()
      );
  };

  const plotSegments = () => {
    segmentPanels
      .selectAll(".segment")
      // Only draw segments that will actually be visible
      .data(
        (d) =>
          d.callers[getActiveCaller()].segments.filter(
            (s) => s.end - s.start > totalLength / width
          ),
        (d) => [d.start, d.end, d.log2]
      )
      .join(
        (enter) =>
          enter
            .append("path")
            .attr("class", "segment")
            .attr("d", (d, i, g) => {
              let j = g[i].parentNode.dataset.index;
              let xScale = xScales[j];
              return `M${xScale(d.start)} ${yScale(
                -yScaleRange - 0.2
              )} L ${xScale(d.end)} ${yScale(-yScaleRange - 0.2)}`;
            })
            .attr("stroke", "orange")
            .attr("stroke-width", 2)
            .call((enter) =>
              enter.transition().attr("d", (d, i, g) => {
                let j = g[i].parentNode.dataset.index;
                let xScale = xScales[j];
                return `M${xScale(d.start)} ${yScale(d.log2)} L ${xScale(
                  d.end
                )} ${yScale(d.log2)}`;
              })
            ),
        (update) =>
          update.attr("d", (d, i, g) => {
            let j = g[i].parentNode.dataset.index;
            let xScale = xScales[j];
            return `M${xScale(d.start)} ${yScale(d.log2)} L ${xScale(
              d.end
            )} ${yScale(d.log2)}`;
          }),
        (exit) =>
          exit
            .transition()
            .attr("d", (d, i, g) => {
              let j = g[i].parentNode.dataset.index;
              let xScale = xScales[j];
              return `M${xScale(d.start)} ${yScale(
                yScaleRange + 0.2
              )} L ${xScale(d.end)} ${yScale(yScaleRange + 0.2)}`;
            })
            .remove()
      );
  };

  // Log ratio grid lines
  lrGrid
    .selectAll(".gridline")
    .data(yScale.ticks())
    .join("line")
    .attr("x1", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[0])
    .attr("x2", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[1])
    .attr("y1", (d) => yScale(d))
    .attr("y2", (d) => yScale(d))
    .attr("class", "gridline");

  // VAF grid lines
  vafGrid
    .selectAll(".gridline")
    .data(yScaleVAF.ticks())
    .join("line")
    .attr("x1", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[0])
    .attr("x2", (d, i, g) => xScales[g[i].parentNode.dataset.index].range()[1])
    .attr("y1", (d) => yScaleVAF(d))
    .attr("y2", (d) => yScaleVAF(d))
    .attr("class", "gridline");

  // VAF
  vafPanels
    .append("g")
    .attr("class", "vaf")
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("data-index", (d, i) => i)
    .selectAll(".point")
    .data((d) => d.vaf)
    .join("circle")
    .attr("cx", (d, i, g) => xScales[g[i].parentNode.dataset.index](d.pos))
    .attr("cy", (d) => yScaleVAF(d.vaf))
    .attr("r", 2)
    .attr("fill", "#333")
    .attr("fill-opacity", 0.3);

  // Clip path to create inner stroke
  const overlayClip = d3.selectAll(".genome-view-area").append("g");
  overlayClip
    .selectAll(".panel-overlay-clip")
    .data(panelWidths)
    .join("clipPath")
    .attr("class", "panel-overlay-clip")
    .attr("id", (d, i) => `panel-${i}-overlay-clip`)
    .append("rect")
    .attr("width", (d) => d)
    .attr("height", panelHeight);

  let selectedChromosomeIndex = 0;
  const selectChromosome = (chromosome) => {
    previousChromosomeIndex = selectedChromosomeIndex;
    selectedChromosomeIndex = cnvData.findIndex(
      (d) => d.chromosome === chromosome
    );
    if (previousChromosomeIndex === selectedChromosomeIndex) {
      return;
    }
    plotArea.selectAll(".panel-overlay").classed("selected", false);
    plotArea
      .selectAll(`.panel-${selectedChromosomeIndex}-overlay`)
      .classed("selected", true);
    chromosomeView.update(cnvData[selectedChromosomeIndex]);
  };

  const overlays = d3.selectAll(".genome-view-area").append("g");
  overlays
    .selectAll(".panel-overlay")
    .data(panelWidths)
    .join("rect")
    .attr(
      "class",
      (d, i) => `panel-overlay panel-${i}-overlay${i === 0 ? " selected" : ""}`
    )
    .attr(
      "transform",
      (d, i) => `translate(${i === 0 ? 0 : d3.sum(panelWidths.slice(0, i))}, 0)`
    )
    .attr("data-index", (d, i) => i)
    .attr("width", (d) => d)
    .attr("height", panelHeight)
    .attr("clip-path", (d, i) => `url(#panel-${i}-overlay-clip)`)
    .attr("fill", "#000")
    .attr("fill-opacity", 0)
    .attr("stroke", "forestgreen")
    .on("mouseenter", (e) => {
      plotArea.selectAll(".panel-overlay").attr("fill-opacity", 0);
      d3.selectAll(`.panel-${e.target.dataset.index}-overlay`).attr(
        "fill-opacity",
        0.2
      );
    })
    .on("mouseout", (e) => {
      d3.selectAll(`.panel-${e.target.dataset.index}-overlay`).attr(
        "fill-opacity",
        0
      );
    })
    .on("click", (e, d, i) =>
      selectChromosome(cnvData[e.target.dataset.index].chromosome)
    );

  // Y axes
  svg
    .append("g")
    .attr("transform", `translate(${margin.left}, ${margin.top})`)
    .attr("class", "y-axis")
    .call(yAxis);

  svg
    .append("g")
    .attr(
      "transform",
      `translate(${margin.left}, ${margin.top + panelHeight + margin.between})`
    )
    .attr("class", "y-axis")
    .call(yAxisVAF);

  // Labels
  vafPanels
    .append("text")
    .attr(
      "transform",
      (d, i) =>
        `translate(${panelWidths[i] / 2},${panelHeight + 10}) rotate(-90)`
    )
    .attr("class", "x-label")
    .text((d, i) => d.label)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "central");

  svg
    .append("text")
    .attr(
      "transform",
      `translate(0,${margin.top + panelHeight / 2}) rotate(-90)`
    )
    .attr("class", "y-label")
    .text("log2 ratio")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  svg
    .append("text")
    .attr(
      "transform",
      `translate(0,${
        margin.top + margin.between + (3 * panelHeight) / 2
      }) rotate(-90)`
    )
    .attr("class", "y-label")
    .text("VAF")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  const update = () => {
    plotRegions();
    plotSegments();
  };

  const getSelectedChromosome = () => selectedChromosomeIndex;

  update();

  return {
    update: update,
    getSelectedChromosome: getSelectedChromosome,
    selectChromosome: selectChromosome,
  };
};

const populateTable = () => {
  const table = d3.select("#cnv-table");
  if (table.empty()) {
    return null;
  }

  const tableHeader = table.append("thead").append("tr");
  const tableBody = table.append("tbody");
  let applyFilter = d3.select("#table-filter-toggle").node().checked;

  let tableData = null;

  const getColumnDef = (col) => {
    switch (col) {
      // Integers
      case "start":
      case "length":
      case "end":
        return {
          class: "right",
          format: (x) => x.toLocaleString(undefined, {}),
        };

      // Floating point, allow for missing numbers
      case "copyNumber":
      case "corrCopyNumber":
      case "baf":
      case "cn":
        return {
          class: "right tooltip-trigger",
          format: (x) => {
            if (x !== null && !isNaN(x)) {
              return x.toLocaleString(undefined, { minimumFractionDigits: 2 });
            }
            return "NA";
          },
        };

      case "view":
        return {
          class: "view-region-link",
          format: (x) => x,
        };

      case "type":
        return {
          class: "left",
          format: (x) => x,
        };

      // Strings
      case "caller":
      case "chromosome":
      case "gene":
      default:
        return {
          class: "left",
          format: (x) => x,
        };
    }
  };

  const getColumnLabel = (col) => {
    const columns = {
      view: "View",
      caller: "Caller",
      chromosome: "Chr",
      genes: "Genes",
      start: "Start",
      end: "End",
      length: "Length",
      type: "Type",
      cn: "CN",
      baf: "BAF",
    };

    if (columns[col]) {
      return columns[col];
    }

    return col;
  };

  const tooltip = d3
    .select(".container")
    .append("div")
    .attr("class", "copy-number-tooltip hidden")
    .call((d) =>
      d
        .append("table")
        .call((t) =>
          t
            .append("thead")
            .append("tr")
            .selectAll("th")
            .data(["caller", "type", "cn"])
            .join("th")
            .text(getColumnLabel)
        )
        .call((t) => t.append("tbody"))
    );

  const copyNumberTooltip = (data, x, y) => {
    tooltip.classed("hidden", false);
    tooltip
      .select("tbody")
      .selectAll("tr")
      .data(data)
      .join("tr")
      .selectAll("td")
      .data((d) => Object.entries(d))
      .join("td")
      .text((d) => getColumnDef(d[0]).format(d[1]));
  };

  const showCopyNumberTooltip = (e) => {
    const tableRow = e.target.parentNode.dataset.index;
    const others = tableData[tableRow].others;

    copyNumberTooltip(others);
  };

  const hideCopyNumberTooltip = (e) => {
    tooltip.classed("hidden", true);
  };

  const update = () => {
    const callerIndex = getActiveCaller();
    tableData = cnvData
      .map((d) =>
        d.callers[callerIndex].cnvs
          .filter((di) => !applyFilter || (applyFilter && di.passed_filter))
          .map((di) => {
            const allCols = { view: "🔍", chromosome: d.chromosome, ...di };
            // Don't display caller and filter status in table
            const { caller, passed_filter, ...cols } = allCols;
            return cols;
          })
      )
      .flat();

    if (tableData.length === 0) {
      tableData = [{ "No data to display": [] }];
    } else {
      // Find the corresponding copy numbers from the other caller(s)
      for (cnv of tableData) {
        // Same chromosome
        let chromData = cnvData.filter((d) => d.chromosome === cnv.chromosome);
        // ... different caller
        let callerData = chromData[0].callers.filter(
          (_d, i) => i !== callerIndex
        );
        // ... same gene
        let otherCnvs = callerData
          .map((d) =>
            d.cnvs
              .filter(
                (c) => cnv.genes.filter((g) => c.genes.includes(g)).length > 0
              )
              .map((c) => {
                return { caller: d.name, type: c.type, cn: c.cn };
              })
          )
          .flat();
        cnv.others = otherCnvs;
      }
    }

    tableHeader
      .selectAll("th")
      .data(
        Object.keys(tableData[0]).filter((k) => k !== "others"),
        (d) => d
      )
      .join("th")
      .text(getColumnLabel)
      .attr("class", (d) => getColumnDef(d).class);

    tableBody
      .selectAll("tr")
      .data(tableData)
      .join("tr")
      .attr("data-chromosome", (d) => d.chromosome)
      .attr("data-start", (d) => d.start)
      .attr("data-length", (d) => d.length)
      .attr("data-index", (_, i) => i)
      .selectAll("td")
      .data((d) => Object.entries(d).filter(([k, _]) => k !== "others"))
      .join("td")
      .text(([key, value]) => getColumnDef(key).format(value))
      .attr("class", ([key, _]) => getColumnDef(key).class);

    tableBody.selectAll(".view-region-link").on("click", (e) => {
      const dataset = e.target.parentElement.dataset;
      zoomToRegion(
        dataset.chromosome,
        Number(dataset.start),
        Number(dataset.start) + Number(dataset.length)
      );
    });

    tableBody
      .selectAll(".tooltip-trigger")
      .on("mouseenter", showCopyNumberTooltip)
      .on("mouseout", hideCopyNumberTooltip)
      .on("mousemove", (e) => {
        const dims = tooltip.node().getBoundingClientRect();
        const offset = d3.select(".container").node().getBoundingClientRect().y;
        const maxHeight =
          document.documentElement.clientHeight - offset - dims.height;
        tooltip
          .style("top", `${Math.min(e.layerY, maxHeight)}px`)
          .style("left", `${e.layerX - dims.width - 20}px`);
      });
  };

  d3.select("#table-filter-toggle").on("change", (event) => {
    applyFilter = event.target.checked;
    update();
  });

  const getData = () => tableData;

  update();

  return {
    getData: getData,
    update: update,
  };
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

const zoomToRegion = (chromosome, start, end, padding = 0.05) => {
  genomeView.selectChromosome(chromosome);
  let bpPadding = (end - start) * padding;
  chromosomeView.zoomTo(start - bpPadding, end + bpPadding);
};

const chromosomeView = new ChromosomePlot({
  element: document.querySelector("#chromosome-view"),
  data: cnvData[0],
});
chromosomeView.addEventListener("zoom", (e) => {
  d3.selectAll(".data-range-warning").classed(
    "hidden",
    !e.detail.dataOutsideRange
  );
});

const genomeView = plotGenomeView();
const resultsTable = populateTable();

d3.select("#chromosome-fit-to-data").on("change", (e) => {
  chromosomeView.fitToData = e.target.checked;
});

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomeView.activeCaller = e.target.value;

  genomeView.update();
  if (resultsTable) {
    resultsTable.update();
  }
});
