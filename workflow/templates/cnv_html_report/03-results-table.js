class ResultsTable extends EventTarget {
  #table;
  #isFiltered;
  #header;
  #body;
  #data;
  #activeCaller;
  #tooltip;

  constructor(element, config) {
    super();

    if (!element) {
      throw Error("no d3 selection supplied");
    }

    this.#table = element;
    this.#isFiltered = config?.filter;

    this.#header = this.#table.append("thead").append("tr");
    this.#body = this.#table.append("tbody");

    if (!config?.data) {
      throw Error("no data supplied");
    }

    // This array controls what columns are displayed in the table,
    // and in what order.
    this.visibleColumns = [
      "view",
      "position",
      "length",
      "genes",
      "type",
      "cn",
      "baf",
    ];

    this.#data = config?.data;
    this.#activeCaller = config?.caller ? config.caller : 0;

    this.#tooltip = this.initTooltip();

    this.update();
  }

  set activeCaller(index) {
    this.#activeCaller = index;
    this.update();
  }

  set filter(isFiltered) {
    this.#isFiltered = isFiltered;
    this.update();
  }

  columnDef(col) {
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
          class: "",
          format: (x) => `<i class="view-region-link bi bi-search">${x}</i>`,
        };

      case "type":
        return {
          class: "left",
          format: (x) => x,
        };

      case "genes":
        return {
          class: "left",
          format: (x) => x.join(", "),
        };

      case "position":
        return {
          class: "left",
          format: (x) =>
            `<span class="clipboard-copy" title="Copy to clipboard">${x}</span>`,
        };

      // Strings
      case "caller":
      case "chromosome":
      default:
        return {
          class: "left",
          format: (x) => x,
        };
    }
  }

  columnLabel(col) {
    const columns = {
      view: "View",
      caller: "Caller",
      chromosome: "Chr",
      genes: "Genes",
      start: "Start",
      end: "End",
      length: "Length",
      position: "Position",
      type: "Type",
      cn: "CN",
      baf: "BAF",
    };

    if (columns[col]) {
      return columns[col];
    }

    return col;
  }

  initTooltip() {
    return d3
      .select("body")
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
              .text(this.columnLabel)
          )
          .call((t) => t.append("tbody"))
      );
  }

  copyNumberTooltip(data) {
    this.#tooltip.classed("hidden", false);
    this.#tooltip
      .select("tbody")
      .selectAll("tr")
      .data(data)
      .join("tr")
      .selectAll("td")
      .data((d) => Object.entries(d))
      .join("td")
      .text((d) => this.columnDef(d[0]).format(d[1]));
  }

  showCopyNumberTooltip(rowIndex, data) {
    this.copyNumberTooltip(data[rowIndex].others);
  }

  hideCopyNumberTooltip() {
    this.#tooltip.classed("hidden", true);
  }

  update() {
    let tableData = this.#data
      .map((d) =>
        d.callers[this.#activeCaller].cnvs
          .filter(
            (di) => !this.#isFiltered || (this.#isFiltered && di.passed_filter)
          )
          .map((di) => {
            const cols = { view: "", chromosome: d.chromosome, ...di };
            const end = (di.start + di.length - 1).toLocaleString();
            cols.position = `${
              d.chromosome
            }:${di.start.toLocaleString()}-${end}`;
            return cols;
          })
      )
      .flat();

    if (tableData.length === 0) {
      tableData = [{ "No data to display": [] }];
    } else {
      // Find the corresponding copy numbers from the other caller(s)
      for (let cnv of tableData) {
        // Same chromosome
        let chromData = this.#data.filter(
          (d) => d.chromosome === cnv.chromosome
        );
        // ... different caller
        let callerData = chromData[0].callers.filter(
          (_, i) => i !== this.#activeCaller
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

    this.#header
      .selectAll("th")
      .data(this.visibleColumns)
      .join("th")
      .text(this.columnLabel)
      .attr("class", (d) => this.columnDef(d).class);

    this.#body
      .selectAll("tr")
      .data(tableData)
      .join("tr")
      .attr("data-chromosome", (d) => d.chromosome)
      .attr("data-start", (d) => d.start)
      .attr("data-length", (d) => d.length)
      .attr("data-index", (_, i) => i)
      .selectAll("td")
      .data((d) =>
        this.visibleColumns.map((c) => {
          return { column: c, value: d[c] };
        })
      )
      .join("td")
      .html((d) => this.columnDef(d.column).format(d.value))
      .attr("class", (d) => this.columnDef(d.column).class);

    this.#body.selectAll(".clipboard-copy").on("click", (e) => {
      const text = e.target.innerHTML;
      navigator.permissions
        .query({ name: "clipboard-write" })
        .then((result) => {
          if (result.state == "granted" || result.state == "prompt") {
            navigator.clipboard
              .writeText(text)
              .catch((res) =>
                console.error("failed to write to clipboard: ", res)
              );
          } else {
            console.warn("permission denied: writing to clipboard");
          }
        });
    });

    this.#body.selectAll(".view-region-link").on("click", (e) => {
      const rowData = e.target.parentNode.parentElement.dataset;
      this.dispatchEvent(
        new CustomEvent("zoom-to-region", {
          detail: {
            chromosome: rowData.chromosome,
            start: Number(rowData.start),
            length: Number(rowData.length),
          },
        })
      );
    });

    this.#body
      .selectAll(".tooltip-trigger")
      .on("mouseenter", (e) =>
        this.showCopyNumberTooltip(e.target.parentNode.dataset.index, tableData)
      )
      .on("mouseout", () => this.hideCopyNumberTooltip())
      .on("mousemove", (e) => {
        const dims = this.#tooltip.node().getBoundingClientRect();
        const offset = d3.select(".container").node().getBoundingClientRect().y;
        const maxHeight =
          document.documentElement.clientHeight - offset - dims.height;
        this.#tooltip
          .style("top", `${Math.min(e.layerY, maxHeight)}px`)
          .style("left", `${e.layerX - dims.width - 20}px`);
      });
  }
}
