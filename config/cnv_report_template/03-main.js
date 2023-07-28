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

const genomeView = new GenomePlot({
  element: document.querySelector("#genome-view"),
  data: cnvData,
});

genomeView.addEventListener("chromosome-change", (e) => {
  chromosomeView.data = cnvData[e.detail.chromosome];
});

const resultsTable = populateTable();

d3.select("#chromosome-fit-to-data").on("change", (e) => {
  chromosomeView.fitToData = e.target.checked;
});

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomeView.activeCaller = e.target.value;
  genomeView.activeCaller = e.target.value;

  if (resultsTable) {
    resultsTable.update();
  }
});
