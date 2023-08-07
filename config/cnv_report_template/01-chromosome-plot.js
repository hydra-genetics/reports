class ChromosomePlot extends EventTarget {
  #data;
  #activeCaller;
  #fitToData;
  #plotArea;
  #lrArea;
  #lrGrid;
  #vafArea;
  #vafGrid;
  #ratios;
  #segments;

  constructor(config) {
    super();

    this.element = config?.element ? config.element : document.body;
    this.name = config?.name ? config.name : "";
    this.#data = config?.data;
    this.#activeCaller = config?.caller ? config.caller : 0;
    this.length = this.#data?.length ? this.#data.length : 0;
    this.zoomRange = [0, this.length];
    this.#fitToData = config?.fitToData ? config.fitToData : false;
    this.animationDuration = config?.animationDuration
      ? config.animationDuration
      : 500;
    this.height = config?.height ? config.height : 400;
    this.width = config?.width ? config.width : 800;
    this.margin = config?.margin
      ? config.margin
      : {
          top: 10,
          right: 30,
          bottom: 40,
          left: 60,
          between: 40,
        };

    this.plotHeight =
      (this.height -
        this.margin.top -
        this.margin.between -
        this.margin.bottom) /
      2;

    this.xScale = d3
      .scaleLinear()
      .domain([0, this.length])
      .range([0, this.width - this.margin.left - this.margin.right]);
    this.ratioYScale = d3.scaleLinear().range([this.plotHeight, 0]);
    this.vafYScale = d3
      .scaleLinear()
      .domain([0, 1])
      .range([this.plotHeight, 0]);

    this.xAxis = (g) => g.call(d3.axisBottom(this.xScale).ticks(5));
    this.ratioYAxis = (g) =>
      g.call(
        d3
          .axisLeft(this.ratioYScale)
          .ticks(8)
          .tickFormat((y, i) => (i % 2 == 0 ? y : ""))
      );
    this.vafYAxis = (g) => g.call(d3.axisLeft(this.vafYScale).ticks(5));

    this.svg = d3
      .select(this.element)
      .attr("preserveAspectRatio", "xMinYMin meet")
      .attr("viewBox", [0, 0, this.width, this.height])
      .attr(
        "style",
        "max-width: 100%; height: auto; max-height: 500px; height: intrinsic;"
      );

    this.svg
      .append("clipPath")
      .attr("id", "lr-area-clip")
      .append("rect")
      .attr("width", this.width - this.margin.left - this.margin.right)
      .attr("height", this.plotHeight);

    this.svg
      .append("clipPath")
      .attr("id", "annotation-clip")
      .append("rect")
      .attr("width", this.width - this.margin.left - this.margin.right)
      .attr("height", 2 * this.plotHeight + this.margin.between);

    this.#plotArea = this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left},${this.margin.top})`)
      .attr("class", "plot-area");

    this.#lrArea = this.#plotArea
      .append("g")
      .attr("id", "lr-plot")
      .attr("clip-path", "url(#lr-area-clip)");
    this.#vafArea = this.#plotArea
      .append("g")
      .attr("id", "vaf-plot")
      .attr("clip-path", "url(#lr-area-clip)")
      .attr(
        "transform",
        `translate(0, ${this.plotHeight + this.margin.between})`
      );

    this.#lrGrid = this.#lrArea.append("g").attr("class", "grid");
    this.#vafGrid = this.#vafArea.append("g").attr("class", "grid");

    this.#ratios = this.#lrArea.append("g").attr("class", "ratios");
    this.#segments = this.#lrArea.append("g").attr("class", "segments");

    this.#drawAxes();
    this.#initializeZoom();
    this.#setLabels();
    this.#drawGridLines();
    this.update();
  }

  set activeCaller(caller) {
    this.#activeCaller = caller;
    this.update();
  }

  get activeCaller() {
    return this.#activeCaller;
  }

  set fitToData(fitToData) {
    this.#fitToData = fitToData;
    this.update();
  }

  get fitToData() {
    return this.#fitToData;
  }

  get data() {
    return this.#data;
  }

  set data(data) {
    this.#data = data;
    return this.update();
  }

  #plotRatios() {
    this.#ratios
      .selectAll(".point")
      .data(this.#data.callers[this.#activeCaller].ratios, (d) => [
        d.start,
        d.end,
        d.log2,
      ])
      .join(
        (enter) =>
          enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d) => this.xScale(d.start))
            .attr("cy", (d) => this.ratioYScale(d.log2))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0)
            .call((enter) =>
              enter
                .transition()
                .duration(this.animationDuration)
                .attr("fill-opacity", 0.3)
            ),
        (update) =>
          update.attr("fill-opacity", 0.3).call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr("cx", (d) => this.xScale(d.start))
              .attr("cy", (d) => this.ratioYScale(d.log2))
          ),
        (exit) => exit.transition().attr("fill-opacity", 0).remove()
      );
  }

  #plotSegments() {
    this.#segments
      .selectAll(".segment")
      .data(this.#data.callers[this.#activeCaller].segments, (d) => [
        d.start,
        d.end,
      ])
      .join(
        (enter) =>
          enter
            .append("path")
            .attr("class", "segment")
            .attr(
              "d",
              (d) =>
                `M${this.xScale(d.start)} ${this.ratioYScale(
                  d.log2
                )} L ${this.xScale(d.end)} ${this.ratioYScale(d.log2)}`
            )
            .attr("stroke", "orange")
            .attr("stroke-width", 2)
            .attr("stroke-opacity", 0)
            .call((enter) => enter.transition().attr("stroke-opacity", 1)),
        (update) =>
          update.attr("stroke-opacity", 1).call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr(
                "d",
                (d) =>
                  `M${this.xScale(d.start)} ${this.ratioYScale(
                    d.log2
                  )} L ${this.xScale(d.end)} ${this.ratioYScale(d.log2)}`
              )
          ),
        (exit) => exit.transition().attr("stroke-opacity", 0).remove()
      );
  }

  #plotVAF() {
    this.#vafArea
      .selectAll(".point")
      .data(this.#data.vaf, (d) => [this.#data.chromosome, d.pos, d.vaf])
      .join(
        (enter) =>
          enter
            .append("circle")
            .attr("class", "point")
            .attr("cx", (d) => this.xScale(d.pos))
            .attr("cy", (d) => this.vafYScale(d.vaf))
            .attr("r", 2)
            .attr("fill", "#333")
            .attr("fill-opacity", 0)
            .call((enter) =>
              enter
                .transition()
                .duration(this.animationDuration)
                .attr("fill-opacity", 0.3)
            ),
        (update) =>
          update.call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr("cx", (d) => this.xScale(d.pos))
          ),
        (exit) =>
          exit
            .transition()
            .duration(this.animationDuration)
            .attr("fill-opacity", 0)
            .remove()
      );
  }

  #plotAnnotations() {
    this.#plotArea
      .selectAll(".annotation")
      .attr("clip-path", "url(#annotation-clip)")
      .data(this.#data.annotations, (d) => [
        this.#data.chromosome,
        d.name,
        d.start,
        d.end,
      ])
      .join(
        (enter) => {
          let annotation_group = enter.append("g").attr("class", "annotation");
          annotation_group
            .append("rect")
            .attr("class", "annotation-marker")
            .attr("x", (d) => this.xScale(d.start))
            .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            .attr("height", this.height - this.margin.top - this.margin.bottom)
            .attr("stroke", "#000")
            .attr("stroke-width", 0.5)
            .attr("fill", "#333")
            .attr("fill-opacity", 0)
            .attr("pointer-events", "none")
            .call((enter) => enter.transition().attr("fill-opacity", 0.05));
          annotation_group
            .append("rect")
            .attr("class", "annotation-label-background")
            .attr("x", (d) => {
              let [labelWidth, _] = getTextDimensions(d.name, "0.8rem");
              return (
                this.xScale(d.start + (d.end - d.start) / 2) -
                labelWidth / 2 -
                5
              );
            })
            .attr("y", (d) => {
              let [_, labelHeight] = getTextDimensions(d.name, "0.8rem");
              return (
                this.plotHeight + this.margin.between / 2 - labelHeight / 2 - 2
              );
            })
            .attr("width", (d) => getTextDimensions(d.name, "0.8rem")[0] + 10)
            .attr("height", (d) => getTextDimensions(d.name, "0.8rem")[1] + 4)
            .attr("fill", "#EEE")
            .attr("rx", 4);
          annotation_group
            .append("text")
            .attr("class", "annotation-label")
            .text((d) => d.name)
            .attr("x", (d) => this.xScale(d.start + (d.end - d.start) / 2))
            .attr("y", this.plotHeight + this.margin.between / 2)
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "central");
          return annotation_group;
        },
        (update) => {
          update.selectAll(".annotation-label").call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr("x", (d) => this.xScale(d.start + (d.end - d.start) / 2))
          );
          update.selectAll(".annotation-label-background").call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr("x", (d) => {
                let [labelWidth, _] = getTextDimensions(d.name, "0.8rem");
                return (
                  this.xScale(d.start + (d.end - d.start) / 2) -
                  labelWidth / 2 -
                  5
                );
              })
          );
          update
            .selectAll(".annotation-marker")
            .attr("fill-opacity", 0.05)
            .call((update) =>
              update
                .transition()
                .duration(this.animationDuration)
                .attr("x", (d) => this.xScale(d.start))
                .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            );
          return update;
        },
        (exit) =>
          exit
            .transition()
            .duration(this.animationDuration)
            .attr("fill-opacity", 0)
            .attr("stroke-opacity", 0)
            .remove()
      );
  }

  #drawAxes() {
    this.svg
      .append("g")
      .attr(
        "transform",
        `translate(${this.margin.left},${this.height - this.margin.bottom})`
      )
      .attr("class", "x-axis")
      .call(this.xAxis);
    this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left},${this.margin.top})`)
      .attr("class", "y-axis")
      .call(this.ratioYAxis);
    this.svg
      .append("g")
      .attr(
        "transform",
        `translate(${this.margin.left},${
          this.margin.top + this.plotHeight + this.margin.between
        })`
      )
      .attr("class", "y-axis")
      .call(this.vafYAxis);
  }

  #initializeZoom() {
    this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`)
      .attr("class", "zoom-layer")
      .append("rect")
      .attr("width", this.width - this.margin.left - this.margin.right)
      .attr("height", this.height - this.margin.top - this.margin.bottom)
      .attr("fill", "transparent")
      .call(
        d3
          .drag()
          .on("start", (e) => {
            d3.select(".zoom-layer")
              .append("rect")
              .attr("class", "zoom-region")
              .attr("pointer-events", "none")
              .attr("x", e.x)
              .attr("width", 0)
              .attr(
                "height",
                this.height - this.margin.bottom - this.margin.top
              )
              .attr("stroke-width", 0)
              .attr("fill", "#333")
              .attr("fill-opacity", 0.1);
          })
          .on("drag", (e) => {
            let leftBound = Math.min(e.x, e.subject.x);
            let width = Math.abs(Math.max(0, e.x) - e.subject.x);
            if (leftBound + width > this.xScale.range()[1]) {
              width = this.xScale.range()[1] - leftBound;
            }
            d3.select(".zoom-region")
              .attr("x", Math.max(0, Math.min(e.x, e.subject.x)))
              .attr("width", width);
          })
          .on("end", (e) => {
            d3.select(".zoom-region").remove();
            let xMin = Math.max(0, Math.min(e.x, e.subject.x));
            let xMax = Math.min(
              this.xScale.range()[1],
              Math.max(e.x, e.subject.x)
            );
            if (xMax - xMin < 3) {
              return;
            }
            this.zoomTo(this.xScale.invert(xMin), this.xScale.invert(xMax));
          })
      )
      .on("click", () => {
        let [xMin, xMax] = this.xScale.domain();
        if (xMax - xMin !== this.#data.length) {
          // Only reset if something actually changed
          this.resetZoom();
        }
      });
  }

  #drawGridLines() {
    this.#lrGrid
      .selectAll(".gridline")
      .data(this.ratioYScale.ticks(), (d) => d)
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("class", "gridline")
            .attr("x1", 0)
            .attr("x2", this.xScale(this.#data.length))
            .attr("y1", (d) => this.ratioYScale(d))
            .attr("y2", (d) => this.ratioYScale(d)),
        (update) =>
          update
            .transition()
            .duration(this.animationDuration)
            .attr("y1", (d) => this.ratioYScale(d))
            .attr("y2", (d) => this.ratioYScale(d)),
        (exit) => exit.remove()
      );

    this.#vafGrid
      .selectAll(".gridline")
      .data(this.vafYScale.ticks(), (d) => d)
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("class", "gridline")
            .attr("x1", 0)
            .attr("x2", this.xScale(this.#data.length))
            .attr("y1", (d) => this.vafYScale(d))
            .attr("y2", (d) => this.vafYScale(d)),
        (update) =>
          update
            .transition()
            .duration(this.animationDuration)
            .attr("y1", (d) => this.vafYScale(d))
            .attr("y2", (d) => this.vafYScale(d)),
        (exit) => exit.remove()
      );
  }

  #setLabels() {
    this.svg
      .append("text")
      .attr("transform", `translate(${this.width / 2},${this.height})`)
      .attr("class", "x-label")
      .text(this.data.label)
      .attr("text-anchor", "middle");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(0,${this.margin.top + this.plotHeight / 2}) rotate(-90)`
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
          this.height - this.margin.bottom - this.plotHeight / 2
        }) rotate(-90)`
      )
      .attr("class", "y-label")
      .text("VAF")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");
  }

  #updateAxes() {
    const [staticYMin, staticYMax] = [-2, 2];

    const [xMin, xMax] = this.zoomRange;
    let yMin, yMax;

    this.xScale.domain([0, this.length]);

    if (xMin > 0 || xMax < this.length) {
      this.xScale.domain([Math.max(xMin, 0), Math.min(xMax, this.length)]);

      const yValues = this.#data.callers[this.#activeCaller].ratios.filter(
        (d) => d.start > xMin && d.start < xMax
      );

      if (yValues.length === 0) {
        yMin = staticYMin;
        yMax = staticYMax;
      } else {
        yMin = yValues.map((d) => d.log2).reduce((a, d) => (d < a ? d : a));
        yMax = yValues.map((d) => d.log2).reduce((a, d) => (d > a ? d : a));
      }
    } else {
      [yMin, yMax] = d3.extent(
        this.#data.callers[this.#activeCaller].ratios,
        (d) => d.log2
      );
      if (!yMin && !yMax) {
        yMin = staticYMin;
        yMax = staticYMax;
      }
    }

    const padding = (yMax - yMin) * 0.05;

    if (this.fitToData) {
      this.dispatchEvent(
        new CustomEvent("zoom", {
          detail: { dataOutsideRange: false },
        })
      );
      this.ratioYScale.domain([yMin - padding, yMax + padding]);
    } else {
      this.dispatchEvent(
        new CustomEvent("zoom", {
          detail: { dataOutsideRange: yMin < staticYMin || yMax > staticYMax },
        })
      );
      this.ratioYScale.domain([staticYMin - padding, staticYMax + padding]);
    }

    this.svg
      .transition()
      .select(".x-axis")
      .duration(this.animationDuration)
      .call(this.xAxis);
    this.svg
      .transition()
      .select(".y-axis")
      .duration(this.animationDuration)
      .call(this.ratioYAxis);
    this.svg.select(".x-label").text(this.#data.label);
  }

  getZoomRange() {
    return this.zoomRange;
  }

  zoomTo(start, end) {
    this.zoomRange = [start, end];
    return this.update();
  }

  resetZoom() {
    this.zoomRange = [0, this.length];
    return this.update();
  }

  update() {
    this.#updateAxes();
    this.#drawGridLines();
    this.#plotRatios();
    this.#plotSegments();
    this.#plotVAF();
    this.#plotAnnotations();
    return this;
  }
}
