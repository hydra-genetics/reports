class ResultsTable extends EventTarget {
  #table;
  #isFiltered;
  #header;
  #body;
  #data;
  #activeCaller;
  #tooltip;
  #tc;
  #baselineOffset;
  #roundSegmentsToInteger;
  #filterConfig;

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
      "adjustedCn",
      "baf",
    ];

    this.#data = config?.data;
    this.#activeCaller = config?.caller ?? 0;
    this.#tc = config?.tc ? config.tc : 1;
    this.#baselineOffset = config?.baselineOffset ? config.baselineOffset : 0;
    this.#roundSegmentsToInteger = config?.roundSegmentsToInteger ? config.roundSegmentsToInteger : false;
    this.#filterConfig = this.#data[0]?.table_filter_config;

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

  setTc(tc) {
    this.#tc = tc;
    this.update();
  }

  setBaselineOffset(dy) {
    this.#baselineOffset = dy;
    this.update();
  }

  setRoundSegmentsToInteger(checked) {
    this.#roundSegmentsToInteger = checked;
    this.update();
  }

  // Same psi-scaling formula as ChromosomePlot/GenomePlot's
  // #toAbsoluteCopyNumber (see 01-chromosome-plot.js for the derivation),
  // except purity is always applied here — unlike the plots' Log2 view,
  // which defaults to assuming 100% purity until "Simulate purity" is
  // checked, the table's Adjusted CN should reflect the sample's actual
  // estimated TC by default (config.tc, i.e. originalTc), not an
  // artificial 100%-purity assumption.
  #toAbsoluteCopyNumber(x) {
    const tc = this.#tc;
    const psi = 2 * 2 ** this.#baselineOffset;
    let adjCopies = (2 ** x * (tc * psi + 2 * (1 - tc)) - 2 * (1 - tc)) / tc;
    if (this.#roundSegmentsToInteger) {
      adjCopies = Math.round(adjCopies);
    }
    return Math.max(adjCopies, MIN_COPY_NUMBER);
  }

  // Segments are a contiguous, non-overlapping partition per chromosome, so
  // usually exactly one matches; the max-overlap tiebreak only matters for a
  // call straddling a segment boundary.
  #findMatchingSegment(segments, start, end) {
    let best = null;
    let bestOverlap = -Infinity;
    for (const seg of segments) {
      if (start <= seg.end && end >= seg.start) {
        const overlap = Math.min(end, seg.end) - Math.max(start, seg.start);
        if (overlap > bestOverlap) {
          bestOverlap = overlap;
          best = seg;
        }
      }
    }
    return best;
  }

  // Recompute a CNV's live absolute copy number from its own caller's raw
  // segments, matching #toAbsoluteCopyNumber's psi-scaling. Falls back to
  // the call's static (server-computed) cn when no segment matches.
  #liveCn(cnv, segments) {
    const matchedSegment = this.#findMatchingSegment(segments, cnv.start, cnv.start + cnv.length - 1);
    return matchedSegment ? this.#toAbsoluteCopyNumber(matchedSegment.log2) : cnv.cn;
  }

  // Recursively evaluate a structured any_of/all_of/not/leaf filter condition
  // (same shape and same three-valued semantics as merge_cnv_json.py's
  // evaluate_filter_condition) against a plain {adjustedCn, baf, length,
  // caller, extra} fields object. A leaf on a field with no data returns
  // null ("unknown") rather than false, so that not(unknown) stays unknown
  // instead of spuriously becoming true - e.g. not(BAF-in-neutral-range)
  // must not pass just because BAF wasn't measured for this call.
  // #firstMatchingGroup only treats an exact `true` result as a match, so
  // null (like false) never matches - no special-casing needed there. A leaf
  // may set `default` to a value used in place of missing data instead of
  // returning null - e.g. a caller that only reports a co-located variant's
  // AF when one exists can set `default: 0` on that leaf, so "no matching
  // variant" is treated as "no evidence of a problem" rather than blocking
  // classification entirely (see config_table_filter.yaml's Twist_AF leaves).
  #evaluateCondition(fields, condition) {
    if ("any_of" in condition) {
      const results = condition.any_of.map((c) => this.#evaluateCondition(fields, c));
      if (results.some((r) => r === true)) return true;
      if (results.some((r) => r === null)) return null;
      return false;
    }
    if ("all_of" in condition) {
      const results = condition.all_of.map((c) => this.#evaluateCondition(fields, c));
      if (results.some((r) => r === false)) return false;
      if (results.some((r) => r === null)) return null;
      return true;
    }
    if ("not" in condition) {
      const result = this.#evaluateCondition(fields, condition.not);
      return result === null ? null : !result;
    }

    const { field, operator, value } = condition;
    const fixed = {
      adjusted_cn: fields.adjustedCn,
      baf: fields.baf,
      length: fields.length,
      caller: fields.caller,
    };
    let cnvValue = field in fixed ? fixed[field] : fields.extra?.[field];
    if (cnvValue === null || cnvValue === undefined) {
      if (!("default" in condition)) return null;
      cnvValue = condition.default;
    }

    switch (operator) {
      case "=": return cnvValue === value;
      case "!=": return cnvValue !== value;
      case ">": return cnvValue > value;
      case "<": return cnvValue < value;
      case ">=": return cnvValue >= value;
      case "<=": return cnvValue <= value;
      default: throw new Error(`Unknown operator: ${operator}`);
    }
  }

  // Returns the name of the first group in table_filter_config (in
  // definition order) whose condition matches these fields, or null if none
  // match. A group may set `filter: false` to be used for Type classification
  // only, excluded when filterOnly is true (the filter toggle's own check).
  #firstMatchingGroup(fields, filterOnly) {
    if (!this.#filterConfig) return null;
    if (fields.adjustedCn === null || fields.adjustedCn === undefined) return null;
    for (const [name, condition] of Object.entries(this.#filterConfig)) {
      if (filterOnly && condition.filter === false) continue;
      if (this.#evaluateCondition(fields, condition)) return name;
    }
    return null;
  }

  #passesTableFilter(fields) {
    return this.#firstMatchingGroup(fields, true) !== null;
  }

  // Live "Type" classification, reusing the same named groups as the filter
  // toggle (including filter:false groups, which only affect Type) so the two
  // stay consistent for whichever groups are shared. Falls back to the call's
  // static (VCF SVTYPE) classification when no table_filter_config is
  // configured, or when there's no usable live CN to classify — necessary
  // since caller SVTYPE vocabularies aren't standardized (cnvkit: DUP/DEL/
  // COPY_NORMAL, GATK: <COPY_GAIN>/<COPY_LOSS>/<COPY_NORMAL>). Once a live CN
  // is available, "no group matched" is trusted as genuine "Copy neutral" -
  // any leaf that would otherwise block a match on missing data (e.g. a
  // caller reporting no BAF) is expected to declare a `default` if it
  // shouldn't be a blocker (see #evaluateCondition).
  #liveType(fields, staticType) {
    if (!this.#filterConfig) return staticType;
    if (fields.adjustedCn === null || fields.adjustedCn === undefined) return staticType;
    const group = this.#firstMatchingGroup(fields, false);
    if (!group) return "Copy neutral";
    const label = this.#filterConfig[group].label;
    if (label) return label;
    return group.charAt(0).toUpperCase() + group.slice(1).replace(/_/g, " ");
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
      case "adjustedCn":
        return {
          class: "right tooltip-trigger",
          format: (x) => {
            if (x !== null && !isNaN(x)) {
              return x.toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2 });
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
          format: (x) => {
            if (!x) return "";
            return x.join(", ");
          },
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
      adjustedCn: "Adjusted CN",
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

    // Find the corresponding copy numbers from the other caller(s), recompute
    // this row's live Adjusted CN, and decide whether it (or an overlapping
    // call from another caller) currently passes the live table filter.
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
      let overlappingOtherCnvs = callerData
        .map((d) =>
          d.cnvs
            .filter(
              (c) => cnv.genes.filter((g) => c.genes.includes(g)).length > 0
            )
            .map((c) => ({ caller: d.name, type: c.type, cn: c.cn, di: c, segments: d.segments }))
        )
        .flat();
      cnv.others = overlappingOtherCnvs.map(({ caller, type, cn }) => ({ caller, type, cn }));

      // Live-recomputed copy number, matching the plots: find this call's
      // own raw segment (same caller, overlapping position) and reapply
      // the current purity/ploidy/rounding settings to its raw log2.
      const segments = chromData[0].callers[this.#activeCaller].segments;
      const matchedSegment = this.#findMatchingSegment(
        segments,
        cnv.start,
        cnv.start + cnv.length - 1
      );
      cnv.adjustedCn = matchedSegment
        ? this.#toAbsoluteCopyNumber(matchedSegment.log2)
        : null;

      const ownFields = {
        adjustedCn: cnv.adjustedCn ?? cnv.cn,
        baf: cnv.baf,
        length: cnv.length,
        caller: cnv.caller,
        extra: cnv.extra,
      };
      cnv.type = this.#liveType(ownFields, cnv.type);

      const ownPasses = this.#passesTableFilter(ownFields);
      // Cross-caller rescue: if an overlapping call from another caller
      // currently passes under its own live-recomputed CN, this row should
      // show as passing too (mirrors merge_cnv_json.py's filter_chr_cnvs).
      const rescuedByOther = overlappingOtherCnvs.some((other) => {
        const otherFields = {
          adjustedCn: this.#liveCn(other.di, other.segments),
          baf: other.di.baf,
          length: other.di.length,
          caller: other.di.caller,
          extra: other.di.extra,
        };
        return this.#passesTableFilter(otherFields);
      });
      cnv.passesLiveFilter = ownPasses || rescuedByOther;
    }

    if (this.#isFiltered && this.#filterConfig) {
      tableData = tableData.filter((cnv) => cnv.passesLiveFilter);
    }

    const hasData = tableData.length > 0;

    this.#header
      .selectAll("th")
      .data(this.visibleColumns)
      .join("th")
      .text(this.columnLabel)
      .attr("class", (d) => this.columnDef(d).class);

    this.#body.selectAll("tr").remove();

    if (!hasData) {
      this.#body
        .append("tr")
        .append("td")
        .text("No data to display")
        .attr("colspan", this.visibleColumns.length);
      return;
    }

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
