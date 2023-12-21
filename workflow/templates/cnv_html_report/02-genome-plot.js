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

    this.panelWidths = cnvData.map(
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
      .attr("style", "max-width: 100%; max-height: 500px; height: auto;");

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
      .attr("data-index", (_, i) => i);

    this.segmentPanels = this.lrPanels
      .append("g")
      .attr("class", "segments")
      .attr("clip-path", (_, i) => `url(#panel-${i}-overlay-clip)`)
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
    this.drawGridLines();
    this.setLabels();
    this.plotRatios();
    this.plotVAF();
    this.plotSegments();
  }

  set activeCaller(caller) {
    this.#activeCaller = caller;
    this.update();
  }

  get activeCaller() {
    return this.#activeCaller;
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

  set activeCaller(caller) {
    this.#activeCaller = caller;
    this.update();
  }

  get activeCaller() {
    return this.#activeCaller;
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
      .attr("class", "gridline");

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
    const nRatioPoints = cnvData.map(
      (d) => d.callers[this.#activeCaller].ratios.length
    );
    const ratioPointsPerChromosome =
      nRatioPoints.reduce((a, b) => a + b, 0) / cnvData.length;

    this.ratioPanels
      .selectAll(".data-point")
      .data(
        (d, i, g) => {
          if (ratioPointsPerChromosome < MAX_POINTS) {
            return d.callers[this.#activeCaller].ratios;
          } else {
            return slidingPixelWindow(
              d.callers[this.#activeCaller].ratios,
              this.xScales[g[i].parentNode.dataset.index],
              "start",
              "log2",
              3,
              true
            );
          }
        },
        (d) => d.start
      )
      .join(
        (enter) => {
          if (enter.data()[0]?.mean) {
            return enter
              .append("g")
              .attr("class", "data-point")
              .call((g) =>
                g
                  .append("line")
                  .attr("class", "mean")
                  .attr("x1", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.start
                    )
                  )
                  .attr("x2", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.end
                    )
                  )
                  .attr("y1", this.ratioYScale.range()[0])
                  .attr("y2", this.ratioYScale.range()[0])
                  .attr("stroke", "#333")
                  .attr("opacity", 0.5)
                  .transition()
                  .duration(this.animationDuration)
                  .attr("y1", (d) => this.ratioYScale(d.mean))
                  .attr("y2", (d) => this.ratioYScale(d.mean))
              )
              .call((g) =>
                g
                  .append("rect")
                  .attr("class", "variance-rect")
                  .attr("x", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.start
                    )
                  )
                  .attr("y", this.ratioYScale.range()[0])
                  .attr("width", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.end - d.start
                    )
                  )
                  .attr("height", 0)
                  .attr("fill", "#333")
                  .attr("fill-opacity", 0.3)
                  .transition()
                  .duration(this.animationDuration)
                  .attr("y", (d) => this.ratioYScale(d.mean + d.sd))
                  .attr("height", (d) =>
                    this.ratioYScale(this.ratioYScale.domain()[1] - 2 * d.sd)
                  )
              )
              .call((g) => g.transition().duration(this.animationDuration));
          }

          return enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d, i, g) =>
              this.xScales[g[i].parentNode.dataset.index](d.start)
            )
            .attr("cy", this.ratioYScale.range()[0])
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0.3)
            .call((enter) =>
              enter
                .transition()
                .duration(this.animationDuration)
                .attr("cy", (d) => this.ratioYScale(d.log2))
            );
        },
        (update) => {
          if (update.data()[0]?.mean) {
            return update
              .call((update) =>
                update
                  .selectAll(".mean")
                  .transition()
                  .duration(this.animationDuration)
                  .attr("x1", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.start
                    )
                  )
                  .attr("x2", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.end
                    )
                  )
                  .attr("y1", (d) => this.ratioYScale(d.mean))
                  .attr("y2", (d) => this.ratioYScale(d.mean))
              )
              .call((update) =>
                update
                  .selectAll(".variance-rect")
                  .transition()
                  .duration(this.animationDuration)
                  .attr("x", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.start
                    )
                  )
                  .attr("y", (d) => this.ratioYScale(d.mean + d.sd))
                  .attr("height", (d) =>
                    this.ratioYScale(this.ratioYScale.domain()[1] - 2 * d.sd)
                  )
              );
          }

          return update.call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr("cx", (d, i, g) =>
                this.xScales[g[i].parentNode.dataset.index](d.start)
              )
              .attr("cy", (d) => this.ratioYScale(d.log2))
          );
        },
        (exit) => {
          exit
            .call((exit) =>
              exit
                .selectAll(".variance-rect")
                .transition()
                .duration(this.animationDuration)
                .attr("y", 0)
                .attr("height", 0)
            )
            .call((exit) =>
              exit
                .selectAll(".mean")
                .transition()
                .duration(this.animationDuration)
                .attr("y1", 0)
                .attr("y2", 0)
            )
            .transition()
            .delay(this.animationDuration)
            .remove();
        }
      );
  }

  plotSegments() {
    this.segmentPanels
      .selectAll(".segment")
      // Only draw segments that will actually be visible
      .data(
        (d) =>
          d.callers[this.#activeCaller].segments.filter(
            (s) => s.end - s.start > this.totalLength / this.width
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
              let xScale = this.xScales[j];
              return `M${xScale(d.start)} ${
                this.ratioYScale.range()[0]
              } L ${xScale(d.end)} ${this.ratioYScale.range()[0]}`;
            })
            .attr("stroke-width", 2)
            .call((enter) =>
              enter
                .transition()
                .duration(this.animationDuration)
                .attr("d", (d, i, g) => {
                  let j = g[i].parentNode.dataset.index;
                  let xScale = this.xScales[j];
                  return `M${xScale(d.start)} ${this.ratioYScale(
                    d.log2
                  )} L ${xScale(d.end)} ${this.ratioYScale(d.log2)}`;
                })
            ),
        (update) =>
          update.attr("d", (d, i, g) => {
            let j = g[i].parentNode.dataset.index;
            let xScale = this.xScales[j];
            return `M${xScale(d.start)} ${this.ratioYScale(d.log2)} L ${xScale(
              d.end
            )} ${this.ratioYScale(d.log2)}`;
          }),
        (exit) =>
          exit
            .transition()
            .duration(this.animationDuration)
            .attr("d", (d, i, g) => {
              let j = g[i].parentNode.dataset.index;
              let xScale = this.xScales[j];
              return `M${xScale(d.start)} ${this.ratioYScale(
                this.ratioYScaleRange + 0.2
              )} L ${xScale(d.end)} ${this.ratioYScale(
                this.ratioYScaleRange + 0.2
              )}`;
            })
            .remove()
      );
  }

  plotVAF() {
    const nVafPoints = this.#data.map((d) => d.vaf.length);
    const vafPointsPerChromosome =
      nVafPoints.reduce((a, b) => a + b, 0) / this.#data.length;

    this.vafPanels
      .append("g")
      .attr("class", "vaf")
      .attr("clip-path", (_, i) => `url(#panel-${i}-overlay-clip)`)
      .attr("data-index", (_, i) => i)
      .selectAll(".data-point")
      .data((d, i, g) => {
        if (vafPointsPerChromosome < MAX_POINTS) {
          return d.vaf;
        } else {
          return slidingPixelWindowVAF(
            d.vaf,
            this.xScales[g[i].parentNode.dataset.index],
            3,
            true
          );
        }
      })
      .join(
        (enter) => {
          if (enter.data()[0]?.mean) {
            return enter
              .append("g")
              .attr("class", "data-point")
              .call((g) =>
                g
                  .append("rect")
                  .attr("class", "variance-rect")
                  .attr("x", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.start
                    )
                  )
                  .attr("y", this.vafYScale.range()[0])
                  .attr("height", 0)
                  .attr("width", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.end - d.start
                    )
                  )
                  .attr("fill", "#333")
                  .attr("fill-opacity", 0.3)
                  .transition()
                  .duration(this.animationDuration)
                  .attr("y", (d) => this.vafYScale(d.mean + d.sd))
                  .attr("height", (d) =>
                    this.vafYScale(this.vafYScale.domain()[1] - 2 * d.sd)
                  )
              )
              .call((g) =>
                g
                  .append("line")
                  .attr("class", "mean-line")
                  .attr("x1", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.start
                    )
                  )
                  .attr("x2", (d, i, g) =>
                    this.xScales[g[i].parentNode.parentNode.dataset.index](
                      d.end
                    )
                  )
                  .attr("y1", this.vafYScale.range()[0])
                  .attr("y2", this.vafYScale.range()[0])
                  .attr("stroke", "#333")
                  .attr("opacity", 0.5)
                  .transition()
                  .duration(this.animationDuration)
                  .attr("y1", (d) => this.vafYScale(d.mean))
                  .attr("y2", (d) => this.vafYScale(d.mean))
              );
          }

          return enter
            .append("circle")
            .attr("class", "data-point")
            .attr("cx", (d, i, g) =>
              this.xScales[g[i].parentNode.dataset.index](d.pos)
            )
            .attr("cy", (d) => this.vafYScale(d.vaf))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0.3);
        },
        (update) => update,
        (exit) => exit.remove()
      );
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

  update() {
    this.plotRatios();
    this.plotSegments();
  }
}
