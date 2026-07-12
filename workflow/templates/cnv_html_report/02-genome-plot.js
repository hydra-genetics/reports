const MAX_POINTS_GENOME = 400;

class GenomePlot extends EventTarget {
  #data;
  #activeCaller;
  #fitToData;
  #plotArea;
  #lrArea;
  #bafArea;
  #ratios;
  #segments;
  #selectedChromosome;
  #canvas;
  #ctx;
  #showAllData;
  #oorPanels;

  constructor(config) {
    super();

    this.element = config?.element ? config.element : document.body;
    this.height = config?.height ? config.height : 400;
    this.widePlotWidth = config?.widePlotWidth ? config.widePlotWidth : false;
    this.width = this.widePlotWidth ? this.widePlotWidth : (config?.width ? config.width : 800);
    this.#data = config?.data;
    this.#showAllData = config?.showAllData ? config.showAllData : false;
    this.baselineOffset = config?.baselineOffset ? config.baselineOffset : 0;
    this.#activeCaller = config?.caller ? config.caller : 0;
    this.#selectedChromosome = config?.selectedChromosome
      ? config.selectedChromosome
      : 0;
    this.animationDuration = config?.animationDuration
      ? config.animationDuration
      : 500;
    this.margin = config?.margin
      ? config.margin
      : {
        top: 10,
        right: 60,
        bottom: 60,
        left: 60,
        between: 20,
      };

    this.simulatePurity = config?.simulatePurity
      ? config.simulatePurity
      : false;
    this.roundSegmentsToInteger = config?.roundSegmentsToInteger
      ? config.roundSegmentsToInteger
      : false;
    this.viewMode = config?.viewMode ? config.viewMode : "log2";
    this.tc = config?.tc ? config.tc : 1;

    this.totalLength = d3.sum(this.#data.map((d) => d.length));

    this.panelWidths = this.#data.map(
      (d) =>
        ((this.width - this.margin.left - this.margin.right) * d.length) /
        this.totalLength
    );
    this.panelHeight =
      (this.height -
        this.margin.top -
        this.margin.bottom -
        this.margin.between) /
      2;

    this.xScales = this.#data.map((d, i) =>
      d3.scaleLinear().domain([0, d.length]).range([0, this.panelWidths[i]])
    );

    this.ratioYScaleRange = 2;
    this.ratioYScale = d3.scaleLinear().range([this.panelHeight, 0]);
    this.#updateRatioRange();

    this.bafYScale = d3
      .scaleLinear()
      .domain([0, 1])
      .range([this.panelHeight, 0]);
    this.ratioYAxis = (g) => g.call(d3.axisLeft(this.ratioYScale).ticks(5));
    this.bafYAxis = (g) => g.call(d3.axisLeft(this.bafYScale).ticks(5));
    this.bafYAxisRight = (g) => g.call(d3.axisRight(this.bafYScale).ticks(5));

    this.cnYScale = d3.scaleLog().base(2).range([this.panelHeight, 0]);
    this.cnYAxis = (g) =>
      g.call(
        d3
          .axisRight(this.viewMode === "copyNumber" ? this.ratioYScale : this.cnYScale)
          .ticks(5)
      );

    this.svg = d3.select("#genome-view");
    if (this.widePlotWidth) {
      this.svg
        .attr("width", this.width)
        .attr("height", this.height)
        .style("max-width", "none")
        .style("height", "auto");
    } else {
      this.svg
        .attr("preserveAspectRatio", "xMinYMin meet")
        .attr("viewBox", [0, 0, this.width, this.height])
        .attr("style", "height: auto;");
    }

    this.#plotArea = this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);

    const lrArea = this.#plotArea.append("g").attr("class", "genome-view-area");
    const bafArea = this.#plotArea
      .append("g")
      .attr("class", "genome-view-area")
      .attr(
        "transform",
        `translate(0,${this.panelHeight + this.margin.between})`
      );

    this.lrPanels = this.addPanels(lrArea);
    this.bafPanels = this.addPanels(bafArea)
      .append("g")
      .attr("clip-path", (_, i) => `url(#panel-${i}-overlay-clip)`);

    this.ratioPanels = this.lrPanels
      .append("g")
      .attr("class", "regions")
      .attr("clip-path", (_, i) => `url(#panel-${i}-overlay-clip)`)
      .attr("data-index", (_, i) => i)
      .attr("data-caller", this.#activeCaller);

    this.segmentPanels = this.lrPanels
      .append("g")
      .attr("class", "segments")
      .attr("clip-path", (_, i) => `url(#panel-${i}-overlay-clip)`)
      .attr("data-index", (_, i) => i)
      .attr("data-caller", this.#activeCaller);

    // Unclipped group for out-of-range indicators (so arrows can poke beyond panel edge)
    this.#oorPanels = this.lrPanels
      .append("g")
      .attr("class", "oor-panels")
      .attr("data-index", (_, i) => i);

    const overlayClip = d3.selectAll(".genome-view-area").append("g");
    overlayClip
      .selectAll(".panel-overlay-clip")
      .data(this.panelWidths)
      .join("clipPath")
      .attr("class", "panel-overlay-clip")
      .attr("id", (_, i) => `panel-${i}-overlay-clip`)
      .append("rect")
      .attr("width", (d) => d)
      .attr("height", this.panelHeight);

    const overlays = d3.selectAll(".genome-view-area").append("g");
    overlays
      .selectAll(".panel-overlay")
      .data(this.panelWidths)
      .join("rect")
      .attr(
        "class",
        (_, i) =>
          `panel-overlay panel-${i}-overlay${i === 0 ? " selected" : ""}`
      )
      .attr(
        "transform",
        (_, i) =>
          `translate(${i === 0 ? 0 : d3.sum(this.panelWidths.slice(0, i))}, 0)`
      )
      .attr("data-index", (_, i) => i)
      .attr("width", (d) => d)
      .attr("height", this.panelHeight)
      .attr("clip-path", (_, i) => `url(#panel-${i}-overlay-clip)`)
      .attr("fill", "#000")
      .attr("fill-opacity", 0)
      .attr("stroke", "forestgreen")
      .on("mouseenter", (e) => {
        this.#plotArea.selectAll(".panel-overlay").attr("fill-opacity", 0);
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
      .on("click", (e) =>
        this.selectChromosome(this.#data[e.target.dataset.index].chromosome)
      );

    this.drawAxes();
    this.#updateCnAxis();
    this.drawGridLines();
    this.setLabels();
    this.#setupCanvas();
    this.plotRatios();
    this.plotSegments();
    this.plotOutOfRangeIndicators();
    this.plotBAF();
  }

  set activeCaller(caller) {
    this.#activeCaller = caller;
    this.ratioPanels.attr("data-caller", caller);
    this.segmentPanels.attr("data-caller", caller);
    this.update();
  }

  get activeCaller() {
    return this.#activeCaller;
  }

  setSimulatePurity(active) {
    this.simulatePurity = active;
    this.#updateRatioRange();
    this.#updateCnAxis();
    this.update();
  }

  /**
   * When simulating purity, widen the ratio row's visible range down to the
   * copy-number floor (instead of the fixed -2..2 log2-ratio window) so
   * near-zero-copy segments/points are always drawn in place rather than
   * relying on the out-of-range indicator. In Copy-number view, use a plain
   * linear [0, CN_VIEW_Y_MAX] range instead.
   */
  #updateRatioRange() {
    if (this.viewMode === "copyNumber") {
      this.ratioYMin = 0;
      this.ratioYMax = CN_VIEW_Y_MAX;
    } else {
      this.ratioYMin = this.simulatePurity ? MIN_LOG2_RATIO : -this.ratioYScaleRange;
      this.ratioYMax = this.ratioYScaleRange;
    }
    this.ratioYScale.domain([this.ratioYMin, this.ratioYMax]);
    // Called once before `drawAxes()` sets up the DOM (during construction,
    // where `drawAxes()`'s own initial render picks up this domain) and
    // again later from setSimulatePurity/setViewMode, where the axis DOM
    // already exists and needs to be refreshed explicitly. Not transitioned:
    // setSimulatePurity(true) triggers this and setViewMode("copyNumber")
    // back-to-back, and a second transition scheduled before the first has
    // started can interrupt it before its tick join ever applies, leaving
    // stale ticks on screen.
    if (this.svg) {
      this.svg.select(".ratio-y-axis").call(this.ratioYAxis);
    }
  }

  setRoundSegmentsToInteger(active) {
    this.roundSegmentsToInteger = active;
    this.update();
  }

  setViewMode(mode) {
    this.viewMode = mode;
    this.#updateRatioRange();
    this.#updateCnAxis();
    this.svg
      .select("#primary-y-label")
      .text(mode === "copyNumber" ? "Copy number" : "log2 ratio");
    this.update();
  }

  setTc(tc) {
    if (tc != this.tc) {
      this.tc = tc;
      this.update();
    }
  }

  #toAbsoluteCopyNumber(x, isSegment) {
    const tc = this.simulatePurity ? this.tc : 1;
    // Ploidy (psi) scales the sample's own average copy content, since coverage
    // ratios are normalised relative to the sample's own median, not a fixed
    // diploid reference. baselineOffset already encodes the assumed ploidy.
    const psi = 2 * 2 ** this.baselineOffset;
    let adjCopies = (2 ** x * (tc * psi + 2 * (1 - tc)) - 2 * (1 - tc)) / tc;
    if (isSegment && this.roundSegmentsToInteger) {
      adjCopies = Math.round(adjCopies);
    }
    return Math.max(adjCopies, MIN_COPY_NUMBER);
  }

  get #slidingWindowOffset() {
    return this.viewMode === "copyNumber" ? 0 : this.baselineOffset;
  }

  get #slidingWindowMinValue() {
    return this.viewMode === "copyNumber" ? MIN_COPY_NUMBER : MIN_LOG2_RATIO;
  }

  transformLog2Ratio(x, isSegment = false) {
    if (x === undefined || x === null || isNaN(x)) return 0;
    const tx = Math.log2(this.#toAbsoluteCopyNumber(x, isSegment) / 2);
    const res = tx - this.baselineOffset;
    return isFinite(res) ? res : (tx < 0 ? -10 : 10);
  }

  transformCopyNumber(x, isSegment = false) {
    if (x === undefined || x === null || isNaN(x)) return 2;
    return this.#toAbsoluteCopyNumber(x, isSegment);
  }

  transformValue(x, isSegment = false) {
    return this.viewMode === "copyNumber"
      ? this.transformCopyNumber(x, isSegment)
      : this.transformLog2Ratio(x, isSegment);
  }

  transformBAF(x) {
    if (x === undefined || x === null || isNaN(x)) return 0.5;
    let tx = x;
    if (this.simulatePurity) {
      tx = (tx - 0.5 * (1 - this.tc)) / this.tc;
      if (tx < 0) {
        tx = 0;
      } else if (tx > 1) {
        tx = 1;
      }
    }
    return isFinite(tx) ? tx : 0.5;
  }
  addPanels(g) {
    const panels = g
      .selectAll(".chromosome-panel")
      .data(this.#data)
      .join("g")
      .attr("data-index", (_, i) => i)
      .attr("class", "chromosome-panel")
      .attr(
        "transform",
        (_, i) =>
          `translate(${i === 0 ? 0 : d3.sum(this.panelWidths.slice(0, i))}, 0)`
      );

    // Panel backgrounds
    panels
      .append("rect")
      .attr("class", "bg-rect")
      .attr("width", (_, i) => this.panelWidths[i])
      .attr("height", this.panelHeight)
      .attr("fill", "#FFF")
      .attr("stroke", "#444");

    return panels;
  }

  drawAxes() {
    this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`)
      .attr("class", "y-axis ratio-y-axis")
      .call(this.ratioYAxis);

    this.svg
      .append("g")
      .attr(
        "transform",
        `translate(${this.width - this.margin.right}, ${this.margin.top})`
      )
      .attr("class", "y-axis cn-y-axis")
      .call(this.cnYAxis);

    this.svg
      .append("g")
      .attr(
        "transform",
        `translate(${this.margin.left}, ${this.margin.top + this.panelHeight + this.margin.between
        })`
      )
      .attr("class", "y-axis")
      .call(this.bafYAxis);

    this.svg
      .append("g")
      .attr(
        "transform",
        `translate(${this.width - this.margin.right}, ${this.margin.top + this.panelHeight + this.margin.between
        })`
      )
      .attr("class", "y-axis")
      .call(this.bafYAxisRight);
  }

  #updateCnAxis() {
    this.cnYScale.domain(
      this.viewMode === "copyNumber"
        ? [this.ratioYMin, this.ratioYMax]
        : [
            cnFromRatio(this.ratioYMin + this.baselineOffset),
            cnFromRatio(this.ratioYMax + this.baselineOffset),
          ]
    );
    // Not transitioned — see #updateRatioRange for why (this is invoked
    // back-to-back with it from setSimulatePurity/setViewMode).
    this.svg.select(".cn-y-axis").call(this.cnYAxis);
  }

  drawGridLines() {
    // Inserted right before .regions (rather than .append()'d, which would
    // render on top / z-order-last) so gridlines sit behind segments/points
    // but still above the panel's bg-rect background — otherwise the bolder
    // integer-CN gridlines visually obscure segments that land exactly on
    // an integer copy number.
    this.lrPanels
      .insert("g", ".regions")
      .attr("class", "grid")
      .attr("data-index", (d, i) => i);

    const bafGrid = this.bafPanels
      .append("g")
      .attr("class", "grid")
      .attr("data-index", (d, i) => i)
      .lower();

    bafGrid
      .selectAll(".gridline")
      .data(this.bafYScale.ticks())
      .join("line")
      .attr(
        "x1",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[0]
      )
      .attr(
        "x2",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[1]
      )
      .attr("y1", (d) => this.bafYScale(d))
      .attr("y2", (d) => this.bafYScale(d))
      .attr("class", "gridline");

    this.#updateLrGridLines();
  }

  /**
   * Redraw the ratio-row gridlines. Normally these mark log2-ratio ticks;
   * in "round segments to integer" mode they instead mark whole-number
   * copy-number positions (mirroring ChromosomePlot's integer gridlines),
   * since log2-spaced lines are not meaningful once segments are snapped
   * to integer CN.
   */
  #updateLrGridLines() {
    const integerMode = this.viewMode === "log2" && this.simulatePurity && this.roundSegmentsToInteger;
    let values, yScale;

    if (integerMode) {
      const [cnMin, cnMax] = this.cnYScale.domain();
      const startN = Math.max(0, Math.ceil(cnMin));
      const endN = Math.floor(cnMax);
      values = [];
      for (let n = startN; n <= endN; n++) {
        values.push(n);
      }
      yScale = this.cnYScale;
    } else if (this.viewMode === "copyNumber") {
      // Whole-number gridlines only — ratioYScale.ticks() with no explicit
      // count can pick a fractional step (e.g. 0.5) that doesn't match the
      // axis's own whole-number ticks.
      const [cnMin, cnMax] = this.ratioYScale.domain();
      const startN = Math.max(0, Math.ceil(cnMin));
      const endN = Math.floor(cnMax);
      values = [];
      for (let n = startN; n <= endN; n++) {
        values.push(n);
      }
      yScale = this.ratioYScale;
    } else {
      values = this.ratioYScale.ticks();
      yScale = this.ratioYScale;
    }

    this.lrPanels
      .select(".grid")
      .selectAll(".gridline")
      .data(values, (d) => d)
      .join("line")
      .attr(
        "x1",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[0]
      )
      .attr(
        "x2",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[1]
      )
      .attr("y1", (d) => yScale(d))
      .attr("y2", (d) => yScale(d))
      .attr("class", (d) => {
        if (integerMode) return "integer-cn-gridline";
        return this.viewMode === "log2" && d === 0 ? "gridline baseline" : "gridline";
      });

    // Copy-number view: "Adjust to ploidy" / baseline-offset marks a
    // reference line at the ploidy-equivalent CN value instead of shifting
    // data (unlike Log2 view's baseline-shift behavior).
    const refCn = this.viewMode === "copyNumber" ? 2 * 2 ** this.baselineOffset : null;
    const [dMin, dMax] = this.ratioYScale.domain();
    const showRef = refCn !== null && refCn >= dMin && refCn <= dMax;

    this.lrPanels
      .select(".grid")
      .selectAll(".baseline-ref-line")
      .data(showRef ? [refCn] : [])
      .join("line")
      .attr("class", "gridline baseline baseline-ref-line")
      .attr(
        "x1",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[0]
      )
      .attr(
        "x2",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[1]
      )
      .attr("y1", (d) => this.ratioYScale(d))
      .attr("y2", (d) => this.ratioYScale(d));
  }

  setLabels() {
    // Labels
    d3.select(this.bafPanels.node().parentNode.parentNode)
      .selectAll(".chromosome-panel")
      .append("text")
      .attr(
        "transform",
        (_, i) =>
          `translate(${this.panelWidths[i] / 2},${this.panelHeight + 10
          }) rotate(-90)`
      )
      .attr("class", "x-label")
      .text((d) => d.label)
      .attr("text-anchor", "end")
      .attr("dominant-baseline", "central");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(0,${this.margin.top + this.panelHeight / 2}) rotate(-90)`
      )
      .attr("id", "primary-y-label")
      .attr("class", "y-label")
      .text(this.viewMode === "copyNumber" ? "Copy number" : "log2 ratio")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(${this.width},${this.margin.top + this.panelHeight / 2}) rotate(90)`
      )
      .attr("class", "y-label cn-y-label")
      .text("Copy number")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(0,${this.margin.top + this.margin.between + (3 * this.panelHeight) / 2
        }) rotate(-90)`
      )
      .attr("class", "y-label")
      .text("BAF")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(${this.width},${this.margin.top + this.margin.between + (3 * this.panelHeight) / 2
        }) rotate(90)`
      )
      .attr("class", "y-label")
      .text("BAF")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");
  }

  plotRatios() {
    const self = this;

    const nRatioPoints = self.#data.map(
      (d) => d.callers[self.#activeCaller].ratios.length
    );
    const ratioPointsPerChromosome =
      nRatioPoints.reduce((a, b) => a + b, 0) / self.#data.length;

    // Draw ratios on Canvas (scatter only)
    this.#ctx.save();
    this.#ctx.translate(this.margin.left, this.margin.top);
    this.#ctx.fillStyle = "#888";
    this.#ctx.globalAlpha = 1.0;

    this.#data.forEach((chromData, i) => {
      const xOffset = i === 0 ? 0 : d3.sum(this.panelWidths.slice(0, i));
      const xScale = this.xScales[i];
      let panelRatios = chromData.callers[this.#activeCaller].ratios.map((d) => {
        let td = { ...d };
        td.log2 = self.transformValue(td.log2);
        return td;
      });

      if (ratioPointsPerChromosome > MAX_POINTS_GENOME) {
        panelRatios = slidingPixelWindow(
          panelRatios,
          xScale,
          "start",
          "log2",
          this.#slidingWindowOffset,
          3,
          true,
          this.#slidingWindowMinValue
        );
      }

      panelRatios.forEach((d) => {
        if (d.mean === undefined) {
          const domain = this.ratioYScale.domain();
          if (d.log2 < domain[0] || d.log2 > domain[1]) return;

          const rawX = xOffset + xScale(d.pos !== undefined ? d.pos : (d.start + d.end) / 2);
          const rawY = this.ratioYScale(d.log2);
          
          const x = Math.round(rawX);
          const y = Math.round(rawY);
          
          if (x >= xOffset && x <= xOffset + this.panelWidths[i]) {
            this.#ctx.beginPath();
            this.#ctx.arc(x, y, 1.0, 0, 2 * Math.PI);
            this.#ctx.fill();
          }
        }
      });
    });
    this.#ctx.restore();

    this.ratioPanels.each(function (panelData, i) {
      let panelRatios = panelData.callers[self.#activeCaller].ratios.filter(
        (p) => {
          const [x0, x1] = self.xScales[i].domain();
          return (p.end ?? p.start) >= x0 && p.start <= x1;
        }
      ).map((d) => {
        let td = { ...d };
        td.log2 = self.transformValue(td.log2);
        return td;
      });

      if (ratioPointsPerChromosome > MAX_POINTS_GENOME) {
        panelRatios = slidingPixelWindow(
          panelRatios,
          self.xScales[i],
          "start",
          "log2",
          self.#slidingWindowOffset,
          3,
          true,
          self.#slidingWindowMinValue
        );
      }
      const svgData = panelRatios.filter(
        (p) => p.mean !== undefined
      );

      d3.select(this)
        .selectAll(".data-point")
        .data(svgData, (d) => `${i}-${d.start}-${d.end}`)
        .join(
          (enter) => {
            // opacity 1 directly, not faded in — see the note in
            // plotSegments() for why an enter-transition on fresh
            // construction can get stuck invisible.
            let g = enter.append("g").attr("class", "data-point").attr("opacity", 1);

            g.append("rect")
              .attr("class", "variance-rect")
              .attr("x", (d) => self.xScales[i](d.start))
              .attr("y", (d) => self.ratioYScale(d.mean + d.sd))
              .attr(
                "width",
                (d) => self.xScales[i](d.end) - self.xScales[i](d.start)
              )
              .attr("height", (d) =>
                self.ratioYScale(self.ratioYScale.domain()[1] - 2 * d.sd)
              )
              .attr("fill", "#888")
              .attr("opacity", 0.3);

            g.append("line")
              .attr("class", "mean-line")
              .attr("x1", (d) => self.xScales[i](d.start))
              .attr("x2", (d) => self.xScales[i](d.end))
              .attr("y1", (d) => self.ratioYScale(d.mean))
              .attr("y2", (d) => self.ratioYScale(d.mean))
              .attr("stroke", "#444")
              .attr("opacity", 0.5);

            g.append("polygon")
              .attr("class", "outlier")
              .attr("points", (d) => {
                const start = self.xScales[i](d.start);
                const x0 = start + (self.xScales[i](d.end) - start) / 2;
                const x1 = x0 - 2;
                const x2 = x0 + 2;
                const y0 = self.ratioYScale.range()[0] - 3;
                const y1 = self.ratioYScale.range()[0] - 6;
                return `${x0},${y0},${x1},${y1},${x2},${y1}`;
              })
              .attr("fill", "red")
              .attr("opacity", (d) => (d.hasOutliers ? 1 : 0));

            return g;
          },
          (update) => {
            update
              .selectAll(".variance-rect")
              .data((d) => [d])
              .transition()
              .duration(self.animationDuration)
              .attr("y", (d) =>
                isNaN(d.mean)
                  ? self.ratioYScale.range()[0]
                  : self.ratioYScale(d.mean + d.sd)
              )
              .attr("height", (d) =>
                isNaN(d.sd)
                  ? 0
                  : self.ratioYScale(self.ratioYScale.domain()[1] - 2 * d.sd)
              );

            update
              .selectAll(".mean-line")
              .data((d) => [d])
              .transition()
              .duration(self.animationDuration)
              .attr("y1", (d) =>
                isNaN(d.mean) ? self.ratioYScale.range()[0] : self.ratioYScale(d.mean)
              )
              .attr("y2", (d) =>
                isNaN(d.mean) ? self.ratioYScale.range()[0] : self.ratioYScale(d.mean)
              );

            update
              .selectAll(".outlier")
              .data((d) => [d])
              .transition()
              .duration(self.animationDuration)
              .attr("opacity", (d) => (d.hasOutliers ? 1 : 0));

            return update;
          },
          (exit) => {
            exit.transition().duration(self.animationDuration).attr("opacity", 0).remove();
          }
        );
    });
  }

  /**
   * Draw red edge-line + arrow for every segment whose log2 ratio falls
   * outside the current static y-axis range ([ratioYMin, ratioYMax], widened
   * to the copy-number floor when simulating purity) in the genome overview plot.
   */
  plotOutOfRangeIndicators() {
    const self = this;
    const staticYMin = this.ratioYMin;
    const staticYMax = this.ratioYMax;
    const arrowSize = 4;
    const lineThickness = 2;

    const arrowPoints = (d, i) => {
      const cx = (self.xScales[i](d.start) + self.xScales[i](d.end)) / 2;
      const isAbove = d.log2 > staticYMax;
      const edgeY = isAbove ? 0 : self.panelHeight;
      const tipY  = isAbove ? -arrowSize * 2.5 : self.panelHeight + arrowSize * 2.5;
      return [
        [cx,             tipY ],
        [cx - arrowSize, edgeY],
        [cx + arrowSize, edgeY],
      ].map(p => p.join(",")).join(" ");
    };

    this.#oorPanels.each(function(panelData, i) {
      const xScale = self.xScales[i];

      const oorSegments = panelData.callers[self.activeCaller].segments
        .filter((d) => {
          const ratios = panelData.callers[self.activeCaller].ratios;
          let count = 0;
          for (let j = 0; j < ratios.length; j++) {
            const r = ratios[j];
            if (r.start >= d.start && r.end <= d.end) {
              count++;
              if (count >= 3) break;
            }
            if (r.start > d.end) break;
          }
          return count >= 3;
        })
        .map((d) => {
          let td = { ...d };
          td.log2 = self.transformValue(td.log2, true);
          return td;
        })
        .filter((d) => d.log2 < staticYMin || d.log2 > staticYMax);

      d3.select(this)
        .selectAll(".oor-indicator")
        .data(oorSegments, (d) => `${i}-${d.start}-${d.end}`)
        .join(
          (enter) => {
            const g = enter.append("g").attr("class", "oor-indicator");

            // Red line at the plot edge
            g.append("line")
              .attr("class", "oor-line")
              .attr("x1", (d) => {
                const x1 = xScale(d.start);
                const x2 = xScale(d.end);
                return x2 - x1 < 2 ? (x1 + x2) / 2 - 1 : x1;
              })
              .attr("x2", (d) => {
                const x1 = xScale(d.start);
                const x2 = xScale(d.end);
                return x2 - x1 < 2 ? (x1 + x2) / 2 + 1 : x2;
              })
              .attr("y1", (d) => d.log2 > staticYMax ? 0 : self.panelHeight)
              .attr("y2", (d) => d.log2 > staticYMax ? 0 : self.panelHeight)
              .attr("stroke", "red")
              .attr("stroke-width", lineThickness)
              .attr("stroke-linecap", "round");

            // Arrow pointing out of bounds
            g.append("polygon")
              .attr("class", "oor-arrow")
              .attr("points", (d) => arrowPoints(d, i))
              .attr("fill", "red");

            return g;
          },
          (update) => {
            update.select(".oor-line")
              .transition().duration(self.animationDuration)
              .attr("x1", (d) => {
                const x1 = xScale(d.start);
                const x2 = xScale(d.end);
                return x2 - x1 < 2 ? (x1 + x2) / 2 - 1 : x1;
              })
              .attr("x2", (d) => {
                const x1 = xScale(d.start);
                const x2 = xScale(d.end);
                return x2 - x1 < 2 ? (x1 + x2) / 2 + 1 : x2;
              })
              .attr("y1", (d) => d.log2 > staticYMax ? 0 : self.panelHeight)
              .attr("y2", (d) => d.log2 > staticYMax ? 0 : self.panelHeight);

            update.select(".oor-arrow")
              .transition().duration(self.animationDuration)
              .attr("points", (d) => arrowPoints(d, i));

            return update;
          },
          (exit) => exit.remove()
        );
    });
  }

  plotSegments() {
    const self = this;

    this.segmentPanels.each(function (panelData, i) {
      let panelSegments = panelData.callers[self.activeCaller].segments
        .filter((d) => d.end - d.start > self.totalLength / self.width)
        .map((d) => {
          let td = { ...d };
          td.log2 = self.transformValue(td.log2, true);
          td.caller = panelData.callers[self.activeCaller].name;
          return td;
        });

      d3.select(this)
        .selectAll(".segment")
        .data(panelSegments, (d) => `${d.caller}-${i}-${d.start}-${d.end}`)
        .join(
          (enter) => {
            return enter
              .append("line")
              .attr("class", "segment")
              .attr("x1", (d) => self.xScales[i](d.start))
              .attr("x2", (d) => self.xScales[i](d.end))
              .attr("y1", (d) => self.ratioYScale(d.log2))
              .attr("y2", (d) => self.ratioYScale(d.log2))
              .attr("stroke-width", 2)
              // No fade-in transition here (unlike the update/exit cases
              // below): a transition scheduled during a fresh page's initial
              // synchronous construction can get stuck barely past opacity 0
              // indefinitely (confirmed via d3.active() still reporting the
              // transition as running long after its duration should have
              // elapsed) — segments would stay invisible until some later,
              // unrelated update() call re-triggers plotSegments(). Setting
              // opacity directly avoids depending on that first frame at all.
              .attr("opacity", 1);
          },
          (update) => {
            return update
              .transition()
              .duration(self.animationDuration)
              .attr("y1", (d) => self.ratioYScale(d.log2))
              .attr("y2", (d) => self.ratioYScale(d.log2));
          },
          (exit) => {
            exit
              .transition()
              .duration(self.animationDuration)
              .attr("opacity", 0)
              .remove();
          }
        );
    });
  }

  plotBAF() {
    const self = this;

    const nBafPoints = self.#data.map((d) => d.baf.length);
    const bafPointsPerChromosome =
      nBafPoints.reduce((a, b) => a + b, 0) / self.#data.length;

    // Draw BAF on Canvas (scatter only)
    this.#ctx.save();
    this.#ctx.translate(
      this.margin.left,
      this.margin.top + this.panelHeight + this.margin.between
    );
    this.#ctx.fillStyle = "#888";
    this.#ctx.globalAlpha = 1.0;

    this.#data.forEach((chromData, i) => {
      const xOffset = i === 0 ? 0 : d3.sum(this.panelWidths.slice(0, i));
      const xScale = this.xScales[i];
      let bafData = chromData.baf.map((d) => {
        let td = { ...d };
        if (td.baf !== undefined) td.baf = self.transformBAF(td.baf);
        if (td.baf_min !== undefined) td.baf_min = self.transformBAF(td.baf_min);
        if (td.baf_max !== undefined) td.baf_max = self.transformBAF(td.baf_max);
        return td;
      });

      if (bafPointsPerChromosome > MAX_POINTS_GENOME) {
        bafData = slidingPixelWindowBAF(bafData, xScale, "pos", 3, true);
      }

      bafData.forEach((d) => {
        if (d.mean === undefined) {
          // Check if this is binned data (has baf_min/baf_max) or unbinned
          if (d.baf_min !== undefined && d.baf_max !== undefined) {
            // Binned data: draw as TWO rectangles (mirrored around 0.5)
            const xStart = Math.round(xScale(d.start !== undefined ? d.start : d.pos));
            const xEnd = Math.round(xScale(d.end !== undefined ? d.end : d.pos));
            const x = xOffset + xStart;

            const rawWidth = xEnd - xStart;
            const width = Math.max(2, rawWidth);
            const xAdjusted = rawWidth < 2 ? x - Math.floor((width - rawWidth) / 2) : x;

            // Draw upper rectangle (above 0.5)
            const bafMaxClamped = Math.min(1.0, d.baf_max);
            const bafMinClamped = Math.max(0.0, d.baf_min);
            const yMinUpper = Math.round(this.bafYScale(bafMaxClamped));
            const yMaxUpper = Math.round(this.bafYScale(bafMinClamped));
            const heightUpper = Math.max(2, yMaxUpper - yMinUpper);

            if (xAdjusted >= xOffset && xAdjusted <= xOffset + this.panelWidths[i]) {
              this.#ctx.fillRect(xAdjusted, yMinUpper, width, Math.max(2, heightUpper));
            }

            // Draw lower rectangle (mirrored below 0.5)
            const baf_min_mirrored = Math.max(0.0, 1 - d.baf_max);
            const baf_max_mirrored = Math.min(1.0, 1 - d.baf_min);
            const yMinLower = Math.round(this.bafYScale(baf_max_mirrored));
            const yMaxLower = Math.round(this.bafYScale(baf_min_mirrored));
            const heightLower = Math.max(1, yMaxLower - yMinLower);

            if (xAdjusted >= xOffset && xAdjusted <= xOffset + this.panelWidths[i]) {
              this.#ctx.fillRect(xAdjusted, yMinLower, width, Math.max(2, heightLower));
            }
          } else {
            // Unbinned data: draw as square
            const rawX = xOffset + xScale(d.pos !== undefined ? d.pos : (d.start + d.end) / 2);
            const rawY = this.bafYScale(Math.max(0.0, Math.min(1.0, d.baf)));
            
            const x = Math.round(rawX);
            const y = Math.round(rawY);
            
            if (x >= xOffset && x <= xOffset + this.panelWidths[i]) {
              this.#ctx.beginPath();
              this.#ctx.arc(x, y, 1.0, 0, 2 * Math.PI);
              this.#ctx.fill();
            }
          }
        }
      });
    });
    this.#ctx.restore();

    this.bafPanels.each(function (panelData, i) {
      let panelBaf = panelData.baf.map((d) => {
        let td = { ...d };
        if (td.baf !== undefined) td.baf = self.transformBAF(td.baf);
        if (td.baf_min !== undefined) td.baf_min = self.transformBAF(td.baf_min);
        if (td.baf_max !== undefined) td.baf_max = self.transformBAF(td.baf_max);
        return td;
      });
      if (bafPointsPerChromosome > MAX_POINTS_GENOME) {
        panelBaf = slidingPixelWindowBAF(panelBaf, self.xScales[i], "pos", 3, true);
      }

      const svgData = panelBaf.filter((p) => p.mean !== undefined);

      d3.select(this)
        .selectAll(".data-point")
        .data(svgData, (d) => `${i}-${d.start}-${d.end}:${d.mean < 0.5 ? "-" : "+"}`)
        .join(
          (enter) => {
            // opacity 1 directly, not faded in — see the note in
            // plotSegments() for why an enter-transition on fresh
            // construction can get stuck invisible.
            let g = enter.append("g").attr("class", "data-point").attr("opacity", 1);

            g.append("rect")
              .attr("class", "variance-rect")
              .attr("x", (d) => self.xScales[i](d.start))
              .attr("y", (d) => self.bafYScale(Math.min(1, d.mean + d.sd)))
              .attr(
                "width",
                (d) => self.xScales[i](d.end) - self.xScales[i](d.start)
              )
              .attr("height", (d) => {
                const yTop = self.bafYScale(Math.min(1, d.mean + d.sd));
                const yBottom = self.bafYScale(Math.max(0, d.mean - d.sd));
                return Math.max(2, yBottom - yTop);
              })
              .attr("fill", "#888")
              .attr("opacity", 0.3);

            g.append("line")
              .attr("class", "mean-line")
              .attr("x1", (d) => self.xScales[i](d.start))
              .attr("x2", (d) => self.xScales[i](d.end))
              .attr("y1", (d) => self.bafYScale(d.mean))
              .attr("y2", (d) => self.bafYScale(d.mean))
              .attr("stroke", "#444")
              .attr("stroke-width", 2)
              .attr("opacity", 0.8);

            return g;
          },
          (update) => {
            update
              .selectAll(".variance-rect")
              .data((d) => [d])
              .transition()
              .duration(self.animationDuration)
              .attr("y", (d) => self.bafYScale(Math.min(1, d.mean + d.sd)))
              .attr("height", (d) => {
                const yTop = self.bafYScale(Math.min(1, d.mean + d.sd));
                const yBottom = self.bafYScale(Math.max(0, d.mean - d.sd));
                return Math.max(2, yBottom - yTop);
              });

            update
              .selectAll(".mean-line")
              .data((d) => [d])
              .transition()
              .duration(self.animationDuration)
              .attr("y1", (d) => self.bafYScale(d.mean))
              .attr("y2", (d) => self.bafYScale(d.mean));

            return update;
          },
          (exit) => exit.transition().duration(self.animationDuration).attr("opacity", 0).remove()
        );
    });
  }

  selectChromosome(chromosome, start, end) {
    const previousChromosomeIndex = this.#selectedChromosome;
    const selectedChromosomeIndex = this.#data.findIndex(
      (d) => d.chromosome === chromosome
    );
    if (previousChromosomeIndex === selectedChromosomeIndex) {
      this.dispatchEvent(
        new CustomEvent("chromosome-zoom", {
          detail: {
            chromosome: this.#selectedChromosome,
            start: start,
            end: end,
          },
        })
      );
      return;
    }
    this.#selectedChromosome = selectedChromosomeIndex;
    this.#plotArea.selectAll(".panel-overlay").classed("selected", false);
    this.#plotArea
      .selectAll(`.panel-${selectedChromosomeIndex}-overlay`)
      .classed("selected", true);
    this.dispatchEvent(
      new CustomEvent("chromosome-change", {
        detail: {
          chromosome: this.#selectedChromosome,
          start: start,
          end: end,
        },
      })
    );
  }

  set showAllData(value) {
    this.#showAllData = value;
    this.update();
  }

  setBaselineOffset(dy) {
    this.baselineOffset = dy;
    this.#updateCnAxis();
    this.update();
  }

  #setupCanvas() {
    const container = document.querySelector("#genome-view-container");
    this.#canvas = document.createElement("canvas");
    this.#canvas.className = "plot-canvas print-visible";
    container.appendChild(this.#canvas);
    this.#ctx = this.#canvas.getContext("2d");

    // Handle high-DPI displays
    const dpr = window.devicePixelRatio || 1;
    this.#canvas.width = this.width * dpr;
    this.#canvas.height = this.height * dpr;
    if (this.widePlotWidth) {
      this.#canvas.style.width = this.width + "px";
    } else {
      this.#canvas.style.width = "100%";
    }
    this.#canvas.style.height = "auto";
    this.#ctx.scale(dpr, dpr);
  }

  #clearCanvas() {
    this.#ctx.clearRect(0, 0, this.width, this.height);
  }

  update() {
    this.#clearCanvas();
    this.#updateLrGridLines();
    this.plotRatios();
    this.plotSegments();
    this.plotOutOfRangeIndicators();
    this.plotBAF();
  }
}
