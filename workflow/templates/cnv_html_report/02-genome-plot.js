class GenomePlot extends EventTarget {
  #data;
  #activeCaller;
  #fitToData;
  #plotArea;
  #lrArea;
  #vafArea;
  #ratios;
  #segments;
  #selectedChromosome;

  constructor(config) {
    super();

    this.element = config?.element ? config.element : document.body;
    this.height = config?.height ? config.height : 400;
    this.width = config?.width ? config.width : 800;
    this.#data = config?.data;
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
          right: 30,
          bottom: 60,
          left: 60,
          between: 20,
        };

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
    this.ratioYScale = d3
      .scaleLinear()
      .domain([-this.ratioYScaleRange, this.ratioYScaleRange])
      .range([this.panelHeight, 0]);

    this.vafYScale = d3
      .scaleLinear()
      .domain([0, 1])
      .range([this.panelHeight, 0]);
    this.ratioYAxis = (g) => g.call(d3.axisLeft(this.ratioYScale).ticks(5));
    this.vafYAxis = (g) => g.call(d3.axisLeft(this.vafYScale).ticks(5));

    this.svg = d3
      .select("#genome-view")
      .attr("preserveAspectRatio", "xMinYMin meet")
      .attr("viewBox", [0, 0, this.width, this.height])
      .attr("style", "height: auto;");

    this.#plotArea = this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);

    const lrArea = this.#plotArea.append("g").attr("class", "genome-view-area");
    const vafArea = this.#plotArea
      .append("g")
      .attr("class", "genome-view-area")
      .attr(
        "transform",
        `translate(0,${this.panelHeight + this.margin.between})`
      );

    this.lrPanels = this.addPanels(lrArea);
    this.vafPanels = this.addPanels(vafArea);

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
    this.drawGridLines();
    this.setLabels();
    this.plotRatios();
    this.plotSegments();
    this.plotVAF();
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
    this.update();
  }

  setTc(tc) {
    if (tc != this.tc) {
      this.tc = tc;
      this.update();
    }
  }

  transformLog2Ratio(x) {
    let tx = x;
    if (this.simulatePurity) {
      const minCopyNumber = 1e-3;
      const adjCopies = (2 * 2 ** x - 2 * (1 - this.tc)) / this.tc;
      tx = Math.log2(Math.max(adjCopies, minCopyNumber) / 2);
    }
    return tx - this.baselineOffset;
  }

  transformVAF(x) {
    let tx = x;
    if (this.simulatePurity) {
      tx = (tx - 0.5 * (1 - this.tc)) / this.tc;
      if (tx < 0) {
        tx = 0;
      } else if (tx > 1) {
        tx = 1;
      }
    }
    return tx;
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
      .attr("stroke", "#333");

    return panels;
  }

  drawAxes() {
    this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`)
      .attr("class", "y-axis")
      .call(this.ratioYAxis);

    this.svg
      .append("g")
      .attr(
        "transform",
        `translate(${this.margin.left}, ${
          this.margin.top + this.panelHeight + this.margin.between
        })`
      )
      .attr("class", "y-axis")
      .call(this.vafYAxis);
  }

  drawGridLines() {
    const lrGrid = this.lrPanels
      .append("g")
      .attr("class", "grid")
      .attr("data-index", (d, i) => i);

    const vafGrid = this.vafPanels
      .append("g")
      .attr("class", "grid")
      .attr("data-index", (d, i) => i);

    lrGrid
      .selectAll(".gridline")
      .data(this.ratioYScale.ticks())
      .join("line")
      .attr(
        "x1",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[0]
      )
      .attr(
        "x2",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[1]
      )
      .attr("y1", (d) => this.ratioYScale(d))
      .attr("y2", (d) => this.ratioYScale(d))
      .attr("class", (d) => {
        return d === 0 ? "gridline baseline" : "gridline";
      });

    vafGrid
      .selectAll(".gridline")
      .data(this.vafYScale.ticks())
      .join("line")
      .attr(
        "x1",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[0]
      )
      .attr(
        "x2",
        (_, i, g) => this.xScales[g[i].parentNode.dataset.index].range()[1]
      )
      .attr("y1", (d) => this.vafYScale(d))
      .attr("y2", (d) => this.vafYScale(d))
      .attr("class", "gridline");
  }

  setLabels() {
    // Labels
    this.vafPanels
      .append("text")
      .attr(
        "transform",
        (_, i) =>
          `translate(${this.panelWidths[i] / 2},${
            this.panelHeight + 10
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
      .attr("class", "y-label")
      .text("log2 ratio")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(0,${
          this.margin.top + this.margin.between + (3 * this.panelHeight) / 2
        }) rotate(-90)`
      )
      .attr("class", "y-label")
      .text("VAF")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");
  }

  plotRatios() {
    const self = this;

    const ratioPointsPerChromosome =
      self.#data
        .map((d) => d.callers[self.activeCaller].ratios.length)
        .reduce((a, b) => a + b, 0) / self.#data.length;

    this.ratioPanels.each(function (panelData, i) {
      let panelRatios = panelData.callers[self.activeCaller].ratios.map((d) => {
        let td = { ...d };
        td.log2 = self.transformLog2Ratio(td.log2);
        return td;
      });

      if (ratioPointsPerChromosome > MAX_POINTS) {
        panelRatios = slidingPixelWindow(
          panelRatios,
          self.xScales[i],
          "start",
          "log2",
          3,
          true
        );
      }

      panelRatios = panelRatios.map((d) => {
        let td = { ...d };
        td.caller = panelData.callers[self.activeCaller].name;
        return td;
      });

      d3.select(this)
        .selectAll(".data-point")
        .data(panelRatios, (d) => {
          const suffix = d.mean ? "summary" : "point";
          return `${d.caller}-${i}-${d.start}-${d.end}-${suffix}`;
        })
        .join(
          (enter) => {
            if (enter.data()[0]?.mean) {
              let g = enter.append("g").attr("class", "data-point");

              g.append("rect")
                .attr("class", "variance-rect")
                .attr("x", (d) => self.xScales[i](d.start))
                .attr("y", (d) => self.ratioYScale(d.mean + d.sd))
                .attr("width", (d) => self.xScales[i](d.end - d.start))
                .attr("height", (d) =>
                  self.ratioYScale(self.ratioYScale.domain()[1] - 2 * d.sd)
                )
                .attr("fill", "#333")
                .attr("fill-opacity", 0.3);

              g.append("line")
                .attr("class", "mean-line")
                .attr("x1", (d) => self.xScales[i](d.start))
                .attr("x2", (d) => self.xScales[i](d.end))
                .attr("y1", (d) => self.ratioYScale(d.mean))
                .attr("y2", (d) => self.ratioYScale(d.mean))
                .attr("stroke", "#333")
                .attr("stroke-width", 1)
                .attr("opacity", 0.3);
              return g;
            } else {
              return enter
                .append("circle")
                .attr("class", "data-point")
                .attr("cx", (d) => self.xScales[i](d.start))
                .attr("cy", (d) => self.ratioYScale(d.log2))
                .attr("r", 2)
                .attr("fill", "#333")
                .attr("fill-opacity", 0.3);
            }
          },
          (update) => {
            if (update.data()[0]?.mean) {
              update
                .selectAll(".variance-rect")
                .data((d) => [d])
                .attr("y", (d) => self.ratioYScale(d.mean + d.sd))
                .attr("height", (d) =>
                  self.ratioYScale(self.ratioYScale.domain()[1] - 2 * d.sd)
                );

              update
                .selectAll(".mean-line")
                .data((d) => [d])
                .attr("y1", (d) => self.ratioYScale(d.mean))
                .attr("y2", (d) => self.ratioYScale(d.mean));
            } else {
              return update.attr("cy", (d) => self.ratioYScale(d.log2));
            }
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
          td.log2 = self.transformLog2Ratio(td.log2);
          return td;
        });

      d3.select(this)
        .selectAll(".segment")
        .data(panelSegments)
        .join(
          (enter) => {
            return enter
              .append("line")
              .attr("class", "segment")
              .attr("x1", (d) => self.xScales[i](d.start))
              .attr("x2", (d) => self.xScales[i](d.end))
              .attr("y1", (d) => self.ratioYScale(d.log2))
              .attr("y2", (d) => self.ratioYScale(d.log2))
              .attr("stroke-width", 2);
          },
          (update) => {
            return update
              .attr("y1", (d) => self.ratioYScale(d.log2))
              .attr("y2", (d) => self.ratioYScale(d.log2));
          },
          (exit) => {
            exit.remove();
          }
        );
    });
  }

  plotVAF() {
    const self = this;

    const nVafPoints = self.#data.map((d) => d.vaf.length);
    const vafPointsPerChromosome =
      nVafPoints.reduce((a, b) => a + b, 0) / self.#data.length;

    this.vafPanels.each(function (panelData, i) {
      let panelVaf = panelData.vaf.map((d) => {
        let td = { ...d };
        td.vaf = self.transformVAF(td.vaf);
        return td;
      });
      let summarisedData = false;
      if (vafPointsPerChromosome > MAX_POINTS) {
        summarisedData = true;
        panelVaf = slidingPixelWindowVAF(panelVaf, self.xScales[i], 3, true);
      }

      d3.select(this)
        .selectAll(".data-point")
        .data(panelVaf, (d) => {
          if (summarisedData) {
            return `${i}-${d.start}-${d.end}`;
          }
          return `${i}-${d.pos}`;
        })
        .join(
          (enter) => {
            if (summarisedData) {
              let g = enter.append("g").attr("class", "data-point");

              g.append("rect")
                .attr("class", "variance-rect")
                .attr("x", (d) => self.xScales[i](d.start))
                .attr("y", (d) => self.vafYScale(d.mean + d.sd))
                .attr("width", (d) => self.xScales[i](d.end - d.start))
                .attr("height", (d) =>
                  self.vafYScale(self.vafYScale.domain()[1] - 2 * d.sd)
                )
                .attr("fill", "#333")
                .attr("fill-opacity", 0.3);

              g.append("line")
                .attr("class", "mean-line")
                .attr("x1", (d) => self.xScales[i](d.start))
                .attr("x2", (d) => self.xScales[i](d.end))
                .attr("y1", (d) => self.vafYScale(d.mean))
                .attr("y2", (d) => self.vafYScale(d.mean))
                .attr("stroke", "#333")
                .attr("opacity", 0.3);

              return g;
            } else {
              return enter
                .append("circle")
                .attr("class", "data-point")
                .attr("cx", (d) => self.xScales[i](d.pos))
                .attr("cy", (d) => self.vafYScale(d.vaf))
                .attr("r", 2)
                .attr("fill", "#333")
                .attr("fill-opacity", 0.3);
            }
          },
          (update) => {
            if (summarisedData) {
              update
                .selectAll(".variance-rect")
                .data((d) => [d])
                .attr("y", (d) => self.vafYScale(d.mean + d.sd))
                .attr("height", (d) =>
                  self.vafYScale(self.vafYScale.domain()[1] - 2 * d.sd)
                );

              update
                .selectAll(".mean-line")
                .data((d) => [d])
                .attr("y1", (d) => self.vafYScale(d.mean))
                .attr("y2", (d) => self.vafYScale(d.mean));

              return update;
            } else {
              return update.attr("cy", (d) => self.vafYScale(d.vaf));
            }
          },
          (exit) => exit.remove()
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

  setBaselineOffset(dy) {
    this.baselineOffset = dy;
    this.update();
  }

  update() {
    this.plotRatios();
    this.plotSegments();
    this.plotVAF();
  }
}
