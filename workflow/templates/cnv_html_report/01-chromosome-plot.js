class YCursor {
  constructor(config) {
    this.parent = config?.element ? config.element : document.body;
    this.fontSize = config?.fontSize ? config.fontSize : "0.8rem";
    this.width = config?.width ? config.width : 100;
    this.labelMargin = config?.labelMargin
      ? config.labelMargin
      : { top: 2, right: 5, bottom: 5, left: 5 };
    this.textHeight = getTextDimensions("0,", this.fontSize)[1];
    this.labelHeight =
      this.textHeight + this.labelMargin.top + this.labelMargin.bottom;
    this.yScale = config?.yScale ? config.yScale : null;
    this.secondaryYScale = config?.secondaryYScale
      ? config.secondaryYScale
      : null;
    this.hidden = config?.hidden ? config.hidden : true;

    if (this.yScale === null) {
      throw Error("scale cannot be null");
    }

    this.cursor = this.parent
      .append("g")
      .lower()
      .attr("class", "cursor")
      .attr("pointer-events", "none")
      .attr("opacity", this.hidden ? 0 : 1);

    this.cursor
      .append("line")
      .attr("class", "cursor-line")
      .attr("x1", 0)
      .attr("x2", this.width)
      .attr("stroke", "black");

    this.cursor
      .append("rect")
      .attr("class", "cursor-label primary")
      .attr("y", -this.labelHeight)
      .attr("height", this.labelHeight)
      .attr("rx", 3)
      .attr("stroke", "black")
      .attr("fill", "white");

    this.cursor
      .append("text")
      .attr("class", "primary")
      .attr("x", 5)
      .attr("fill", "black")
      .attr("font-size", this.fontSize)
      .attr("alignment-baseline", "baseline");

    if (this.secondaryYScale !== null) {
      this.cursor
        .append("rect")
        .attr("class", "cursor-label secondary")
        .attr("y", -this.labelHeight)
        .attr("height", this.labelHeight)
        .attr("rx", 3)
        .attr("stroke", "black")
        .attr("fill", "white");

      this.cursor
        .append("text")
        .attr("class", "secondary")
        .attr("fill", "black")
        .attr("font-size", this.fontSize)
        .attr("alignment-baseline", "baseline");
    }
  }

  show() {
    this.hidden = false;
    this.cursor.attr("opacity", 1);
  }

  hide() {
    this.hidden = true;
    this.cursor.attr("opacity", 0);
  }

  set(y) {
    let label = this.yScale.invert(y).toLocaleString();
    let textWidth = getTextDimensions(label, this.fontSize)[0];
    let labelWidth = textWidth + this.labelMargin.left + this.labelMargin.right;

    this.cursor.attr("opacity", 1).attr("transform", `translate(0,${y})`);

    if (this.secondaryYScale !== null) {
      let secondaryLabel = this.secondaryYScale.invert(y).toLocaleString();
      let secondaryTextWidth = getTextDimensions(
        secondaryLabel,
        this.fontSize
      )[0];
      let secondaryLabelWidth =
        secondaryTextWidth + this.labelMargin.left + this.labelMargin.right;

      this.cursor
        .select(".cursor-label.secondary")
        .attr("width", secondaryLabelWidth)
        .attr(
          "transform",
          `translate(${this.width + 5}, ${this.labelHeight / 2})`
        );
      this.cursor
        .select("text.secondary")
        .attr(
          "transform",
          `translate(${this.width + 5 + this.labelMargin.left}, ${this.labelHeight / 2 - this.labelMargin.bottom
          })`
        )
        .text(secondaryLabel);
    }

    this.cursor
      .select(".cursor-label.primary")
      .attr("width", labelWidth)
      .attr(
        "transform",
        `translate(${-(5 + labelWidth)}, ${this.labelHeight / 2})`
      );
    this.cursor
      .select("text.primary")
      .attr(
        "transform",
        `translate(${-(5 + labelWidth)}, ${this.labelHeight / 2 - this.labelMargin.bottom
        })`
      )
      .text(label);
  }
}

class XCursor {
  constructor(config) {
    this.parent = config?.element ? config.element : document.body;
    this.fontSize = config?.fontSize ? config.fontSize : "0.8rem";
    this.height = config?.height ? config.height : 100;
    this.labelMargin = config?.labelMargin
      ? config.labelMargin
      : { top: 2, right: 5, bottom: 5, left: 5 };
    this.labelHeight = getTextDimensions("0,", this.fontSize)[1];
    this.xScale = config?.xScale ? config.xScale : null;
    this.hidden = true;

    if (this.xScale === null) {
      throw Error("scale cannot be null");
    }

    this.cursor = this.parent
      .append("g")
      .lower()
      .attr("class", "vertical-cursor cursor")
      .attr("opacity", this.hidden ? 0 : 1);

    this.cursor
      .append("line")
      .attr("class", "cursor-line")
      .attr("y1", 0)
      .attr("y2", this.height)
      .attr("stroke", "black");

    this.cursor
      .append("rect")
      .attr("class", "cursor-label")
      .attr("y", this.height + 5)
      .attr(
        "height",
        this.labelHeight + this.labelMargin.top + this.labelMargin.bottom
      )
      .attr("rx", 3)
      .attr("stroke", "black")
      .attr("fill", "white");

    this.cursor
      .append("text")
      .attr("x", 5)
      .attr("y", this.height + 5 + this.labelMargin.top + this.labelHeight)
      .attr("fill", "black")
      .attr("font-size", this.fontSize);
  }

  show() {
    this.hidden = false;
    this.cursor.attr("opacity", 1);
  }

  hide() {
    this.hidden = true;
    this.cursor.attr("opacity", 0);
  }

  set(x) {
    let verticalLabel = Math.floor(this.xScale.invert(x)).toLocaleString();
    let verticalLabelWidth = getTextDimensions(verticalLabel, this.fontSize)[0];

    this.cursor.attr("opacity", 1).attr("transform", `translate(${x}, 0)`);

    this.cursor
      .select(".cursor-label")
      .attr("x", -verticalLabelWidth / 2 - this.labelMargin.left)
      .attr(
        "width",
        verticalLabelWidth + this.labelMargin.left + this.labelMargin.right
      );
    this.cursor
      .select("text")
      .attr("x", -verticalLabelWidth / 2)
      .text(verticalLabel);
  }
}

class ChromosomePlot extends EventTarget {
  #data;
  #activeCaller;
  #fitToData;
  #cytobands;
  #plotArea;
  #lrArea;
  #bafArea;
  #ratios;
  #segments;
  #showAllData;
  #equalDistance;
  #highlightedRegion;
  #canvas;
  #ctx;

  constructor(config) {
    super();

    this.#showAllData = config?.showAllData ? config.showAllData : false;
    this.#equalDistance = config?.equalDistance ? config.equalDistance : false;
    this.#highlightedRegion = null;
    this.element = config?.element ? config.element : document.body;
    this.name = config?.name ? config.name : "";
    this.#data = config?.data;
    this.#activeCaller = config?.caller ? config.caller : 0;
    this.zoomRange = [0, this.length];
    this.minZoomRange = config?.minZoomRange ? config.minZoomRange : 20;
    this.#fitToData = config?.fitToData ? config.fitToData : false;
    this.baselineOffset = config?.baselineOffset ? config.baselineOffset : 0;
    this.simulatePurity = config?.simulatePurity
      ? config.simulatePurity
      : false;
    this.tc = config?.tc ? config.tc : 1;
    this.animationDuration = config?.animationDuration
      ? config.animationDuration
      : 500;
    this.height = config?.height ? config.height : 400;
    this.width = config?.width ? config.width : 800;
    this.margin = config?.margin
      ? config.margin
      : {
        top: 10,
        right: 60,
        bottom: 40,
        left: 60,
        between: 40,
      };

    this.cytobandHeight = this.#data.cytobands ? 10 : 0;
    if (this.cytobandHeight > 0) {
      this.margin.top += this.cytobandHeight + this.margin.top;
    }

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
    this.cnYScale = d3.scaleLog().base(2).range([this.plotHeight, 0]);
    this.bafYScale = d3
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
    this.cnYAxis = (g) => g.call(d3.axisRight(this.cnYScale).ticks(5));
    this.bafYAxis = (g) => g.call(d3.axisLeft(this.bafYScale).ticks(5));

    this.svg = d3
      .select(this.element)
      .attr("preserveAspectRatio", "xMinYMin meet")
      .attr("viewBox", [0, 0, this.width, this.height])
      .attr("style", "height: auto; height: intrinsic;");

    this.#drawAxes();

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

    if (this.#data.cytobands) {
      this.#cytobands = this.svg
        .append("g")
        .attr("clip-path", "url(#lr-area-clip)")
        .attr("transform", `translate(${this.margin.left},10)`)
        .attr("class", "cytobands-area");
    }

    this.#lrArea = this.#plotArea
      .append("g")
      .attr("id", "lr-plot")
      .attr("clip-path", "url(#lr-area-clip)");
    this.#bafArea = this.#plotArea
      .append("g")
      .attr("id", "baf-plot")
      .attr("clip-path", "url(#lr-area-clip)")
      .attr(
        "transform",
        `translate(0, ${this.plotHeight + this.margin.between})`
      );

    this.#ratios = this.#lrArea
      .append("g")
      .attr("class", "ratios")
      .attr("data-chromosome", this.data.chromosome)
      .attr("data-caller", this.#activeCaller);
    this.#segments = this.#lrArea
      .append("g")
      .attr("class", "segments")
      .attr("data-chromosome", this.data.chromosome)
      .attr("data-caller", this.#activeCaller);

    this.#setLabels();

    this.annotations = this.#plotArea
      .append("g")
      .attr("class", "annotation-container");

    this.initialiseMousetrap();

    // Setup Canvas for high-performance scatter rendering
    this.#setupCanvas();

    this.update();
  }

  set activeCaller(caller) {
    this.#activeCaller = caller;
    this.#ratios.attr("data-caller", caller);
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

  set showAllData(x) {
    this.#showAllData = x;
    this.update();
  }

  get showAllData() {
    return this.#showAllData;
  }

  get data() {
    return this.#data;
  }

  setData(data, start, end) {
    const prevChromosome = this.data.chromosome;

    if (data && data.chromosome !== prevChromosome) {
      this.#highlightedRegion = null;
      this.#data = data;
      if (!start && !end) {
        this.resetZoom();
      }
      this.#drawAxes();
    }

    if (start || end) {
      if (this.equalDistance) {
        start = this.getRatioIndex(start);
        end = this.getRatioIndex(end);
      }
      this.zoomTo(start, end);
    }

    this.#ratios.attr("data-chromosome", this.data.chromosome);
    this.#segments.attr("data-chromosome", this.data.chromosome);

    this.update();
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

  transformBAF(x) {
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

  get length() {
    if (this.equalDistance) {
      return this.#data.callers[this.#activeCaller].ratios.length;
    }
    return this.#data.length;
  }

  set equalDistance(equalDistance) {
    this.#equalDistance = equalDistance;
    this.resetZoom();
    this.update();
  }

  get equalDistance() {
    return this.#equalDistance;
  }

  getRatioIndex(pos) {
    const ratios = this.#data.callers[this.#activeCaller].ratios;
    let l = 0;
    let r = ratios.length - 1;
    let mid;
    while (l <= r) {
      mid = Math.floor((l + r) / 2);
      if (ratios[mid].start < pos) {
        l = mid + 1;
      } else if (ratios[mid].start > pos) {
        r = mid - 1;
      } else {
        return mid;
      }
    }
    return l < ratios.length ? l : ratios.length - 1;
  }

  #plotCytobands() {
    if (!this.#cytobands) {
      return;
    }

    const cytobandPolygon = ({ start, end, direction }) => {
      const s = this.equalDistance ? this.getRatioIndex(start) : start;
      const e = this.equalDistance ? this.getRatioIndex(end) : end;
      const bandWidthPx = this.xScale(e) - this.xScale(s);
      const arrowSize = Math.min(5, bandWidthPx);
      let points;
      switch (direction) {
        case "right":
          points = [
            [this.xScale(s), 0],
            [this.xScale(e) - arrowSize, 0],
            [this.xScale(e), this.cytobandHeight / 2],
            [this.xScale(e) - arrowSize, this.cytobandHeight],
            [this.xScale(s), this.cytobandHeight],
          ];
          break;
        case "left":
          points = [
            [this.xScale(e), 0],
            [this.xScale(s) + arrowSize, 0],
            [this.xScale(s), this.cytobandHeight / 2],
            [this.xScale(s) + arrowSize, this.cytobandHeight],
            [this.xScale(e), this.cytobandHeight],
          ];
          break;
        case "none":
          points = [
            [this.xScale(s), 0],
            [this.xScale(e), 0],
            [this.xScale(e), this.cytobandHeight],
            [this.xScale(s), this.cytobandHeight],
          ];
          break;
        default:
          throw new Error(`invalid cytoband direction: ${direction}`);
      }
      return points.join(" ");
    };

    const cytobandLabel = (d) => {
      const s = this.equalDistance ? this.getRatioIndex(d.start) : d.start;
      const e = this.equalDistance ? this.getRatioIndex(d.end) : d.end;
      const xDomain = this.xScale.domain();
      let visibleWidth, xPos;
      if (
        xDomain[0] > s &&
        xDomain[0] < e &&
        xDomain[1] > s &&
        xDomain[1] < e
      ) {
        // both ends outside
        visibleWidth = xDomain[1] - xDomain[0];
        xPos = xDomain[0] + visibleWidth / 2;
      } else if (xDomain[0] > s && xDomain[0] < e) {
        // hanging out on the left side
        visibleWidth = e - xDomain[0];
        xPos = xDomain[0] + visibleWidth / 2;
      } else if (xDomain[1] > s && xDomain[1] < e) {
        // hanging out on the right side
        visibleWidth = xDomain[1] - s;
        xPos = s + visibleWidth / 2;
      } else if (xDomain[0] > e || xDomain[1] < s) {
        // the band is out of view
        visibleWidth = 0;
        xPos = s;
      } else {
        // the whole band is visible
        visibleWidth = e - s;
        xPos = s + visibleWidth / 2;
      }

      const labelWidth = getTextDimensions(d.name, this.cytobandHeight)[0];

      if (this.xScale(xDomain[0] + visibleWidth) < labelWidth + 4) {
        return {
          label: "",
          x: xPos,
        };
      }

      return {
        label: d.name,
        x: xPos,
      };
    };

    const colorBrightness = (color) => {
      let colorType;
      if (color.indexOf("#") === 0) {
        colorType = "hex";
      } else if (color.indexOf("rgb") === 0) {
        colorType = "rgb";
      } else {
        throw new Error(`unknown color type, should be hex or rgb: ${color}`);
      }

      let m, r, g, b;

      if (colorType === "hex") {
        if (color.length === 7) {
          m = color.substr(1).match(/(\S{2})/g);
        } else if (color.length === 4) {
          m = color.substr(1).match(/(\S{1})/g);
        } else {
          throw new Error(`invalid hex color value: ${color}`);
        }

        if (m) {
          r = parseInt(m[0], 16);
          g = parseInt(m[1], 16);
          b = parseInt(m[2], 16);
        }
      }

      if (colorType === "rgb") {
        m = color.match(/(\d+){3}/g);
        if (m) {
          r = parseInt(m[0]);
          g = parseInt(m[1]);
          b = parseInt(m[2]);
        }
      }

      if (r !== undefined) {
        return (r * 299 + g * 587 + b * 114) / 1000;
      }
    };

    this.#cytobands
      .selectAll(".cytoband-group")
      .data(this.data.cytobands, (d) => [d.name, d.start, d.end, d.giemsa])
      .join(
        (enter) => {
          return enter
            .append("g")
            .classed("cytoband-group", true)
            .call((g) =>
              g
                .append("polygon")
                .attr("points", cytobandPolygon)
                .attr("fill", (d) => d.color)
                .classed("cytoband", true)
                .classed("centromere", (d) => d.giemsa === "acen")
                .on("mouseenter mousemove", (e, d) => {
                  this.showCytobandName(d.name, e);
                })
                .on("mouseout", () => {
                  this.showCytobandName(null);
                })
            )
            .call((g) =>
              g
                .append("text")
                .attr("x", (d) => this.xScale(cytobandLabel(d).x))
                .attr("y", this.cytobandHeight - 2)
                .attr("fill", (d) =>
                  colorBrightness(d.color) < 128 ? "#ffffff" : "#000000"
                )
                .attr("text-anchor", "middle")
                .classed("label", true)
                .style("font-size", `${this.cytobandHeight}px`)
                .style("pointer-events", "none")
                .text((d) => cytobandLabel(d).label)
            );
        },
        (update) => {
          update
            .call((g) =>
              g
                .select(".label")
                .transition()
                .duration(this.animationDuration)
                .attr("x", (d) => this.xScale(cytobandLabel(d).x))
                .text((d) => cytobandLabel(d).label)
            )
            .call((g) =>
              g
                .select(".cytoband")
                .transition()
                .duration(this.animationDuration)
                .attr("points", cytobandPolygon)
            );
        },
        (exit) => exit.remove()
      );
  }

  showCytobandName(name, event) {
    event?.preventDefault();
    d3.select("body")
      .selectAll(".cytoband-name")
      .data(name ? [name] : [], (d) => d)
      .join(
        (enter) => {
          let div = enter
            .append("div")
            .classed("cytoband-name", true)
            .style("position", "fixed")
            .attr("pointer-events", "none");
          div.append("p").text((d) => d);
          if (event?.clientX && event?.clientY) {
            div = div
              .style("top", `${event.clientY + 30}px`)
              .style("left", `${event.clientX}px`);
          }
          return div;
        },
        (update) =>
          update
            .style("top", `${event?.clientY + 20}px`)
            .style("left", `${event?.clientX}px`)
            .text((d) => d),
        (exit) => exit.remove()
      );
  }

  #plotRatios() {
    const self = this;

    let ratioData;
    if (this.equalDistance) {
      const ratios = this.#data.callers[this.#activeCaller].ratios;
      const [xMin, xMax] = this.xScale.domain();
      const startIdx = Math.max(0, Math.floor(xMin));
      const endIdx = Math.min(ratios.length, Math.ceil(xMax));

      ratioData = ratios.slice(startIdx, endIdx).map((p, i) => {
        let tp = { ...p };
        tp.log2 = self.transformLog2Ratio(tp.log2);
        tp.realStart = tp.start;
        tp.start = startIdx + i;
        tp.end = startIdx + i + 1;
        return tp;
      });
    } else {
      ratioData = this.#data.callers[this.#activeCaller].ratios
        .filter(
          (p) =>
            p.start >= this.xScale.domain()[0] &&
            p.start <= this.xScale.domain()[1]
        )
        .map((p) => {
          let tp = { ...p };
          tp.log2 = self.transformLog2Ratio(tp.log2);
          return tp;
        });
    }

    if (ratioData.length > MAX_POINTS && !this.#showAllData) {
      ratioData = slidingPixelWindow(
        ratioData,
        this.xScale,
        "start",
        "log2",
        0,
        3,
        true
      );
    }

    ratioData = ratioData.map((d) => {
      let td = { ...d };
      td.caller = self.activeCaller;
      return td;
    });

    // Draw ratios on Canvas (scatter only)
    this.#ctx.save();
    this.#ctx.translate(this.margin.left, this.margin.top);
    this.#ctx.fillStyle = "#333";
    this.#ctx.globalAlpha = 0.3;

    ratioData.forEach((d) => {
      if (d.mean === undefined) {
        const x = this.xScale(this.equalDistance ? (d.start + d.end) / 2 : (d.pos !== undefined ? d.pos : (d.start + d.end) / 2));
        const y = this.ratioYScale(d.log2);
        if (x >= 0 && x <= this.width - this.margin.left - this.margin.right) {
          this.#ctx.beginPath();
          this.#ctx.arc(x, y, 2, 0, 2 * Math.PI);
          this.#ctx.fill();
        }
      }
    });
    this.#ctx.restore();

    // Use SVG only for summarized data (means/rects)
    const svgData = ratioData.filter((d) => d.mean !== undefined);

    this.#ratios
      .selectAll(".data-point")
      .data(svgData, (d) => {
        const suffix = "summary";
        return `${d.caller}-${self.data.chromosome}-${d.start}-${d.end}-${suffix}`;
      })
      .join(
        (enter) => {
          // Summarised data
          let g = enter
            .append("g")
            .attr("class", "data-point")
            .attr("opacity", 0);

          g.append("rect")
            .attr("class", "variance-rect")
            .attr("x", (d) => this.xScale(d.start))
            .attr("y", (d) => this.ratioYScale(d.mean ? d.mean + d.sd : 0))
            .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            .attr("height", (d) =>
              isNaN(d.sd)
                ? 0
                : this.ratioYScale(this.ratioYScale.domain()[1] - 2 * d.sd)
            )
            .attr("fill", "#333")
            .attr("opacity", (d) => (isNaN(d.mean) ? 0 : 0.3));

          g.append("line")
            .attr("class", "mean")
            .attr("x1", (d) => this.xScale(d.start))
            .attr("x2", (d) => this.xScale(d.end))
            .attr("y1", (d) =>
              isNaN(d.mean) ? this.ratioYScale.range()[0] : this.ratioYScale(d.mean)
            )
            .attr("y2", (d) =>
              isNaN(d.mean) ? this.ratioYScale.range()[0] : this.ratioYScale(d.mean)
            )
            .attr("stroke", "#333")
            .attr("opacity", (d) => (isNaN(d.mean) ? 0 : 0.3));

          g.append("polygon")
            .attr("class", "outlier")
            .attr("points", (d) => {
              const start = self.xScale(d.start);
              const x0 = start + (self.xScale(d.end) - start) / 2;
              const x1 = x0 - 2;
              const x2 = x0 + 2;
              const y0 = self.ratioYScale.range()[0];
              const y1 = self.ratioYScale.range()[0] - 3;
              return `${x0},${y0},${x1},${y1},${x2},${y1}`;
            })
            .attr("fill", "red")
            .attr("opacity", (d) => (d.hasOutliers ? 1 : 0));

          return g.transition().duration(this.animationDuration).attr("opacity", 1);
        },
        (update) => {
          update
            .selectAll(".mean")
            .data((d) => [d])
            .transition()
            .duration(this.animationDuration)
            .attr("x1", (d) => this.xScale(d.start))
            .attr("x2", (d) => this.xScale(d.end))
            .attr("y1", (d) =>
              isNaN(d.mean) ? this.ratioYScale.range()[0] : this.ratioYScale(d.mean)
            )
            .attr("y2", (d) =>
              isNaN(d.mean) ? this.ratioYScale.range()[0] : this.ratioYScale(d.mean)
            )
            .attr("opacity", (d) => (isNaN(d.mean) ? 0 : 0.3));

          update
            .selectAll(".variance-rect")
            .data((d) => [d])
            .transition()
            .duration(this.animationDuration)
            .attr("x", (d) => this.xScale(d.start))
            .attr("y", (d) =>
              isNaN(d.mean)
                ? this.ratioYScale.range()[0]
                : this.ratioYScale(d.mean + d.sd)
            )
            .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            .attr("height", (d) =>
              isNaN(d.sd)
                ? 0
                : this.ratioYScale(this.ratioYScale.domain()[1] - 2 * d.sd)
            )
            .attr("opacity", (d) => (isNaN(d.mean) ? 0 : 0.3));

          update
            .selectAll(".outlier")
            .data((d) => [d])
            .transition()
            .duration(this.animationDuration)
            .attr("opacity", (d) => (d.hasOutliers ? 1 : 0));

          return update;
        },
        (exit) => {
          exit.transition().duration(this.animationDuration).attr("opacity", 0).remove();
        }
      );
  }

  #plotSegments() {
    const self = this;
    this.#segments
      .selectAll(".segment")
      .data(
        this.#data.callers[this.#activeCaller].segments.map((d) => {
          let ts = { ...d };
          ts.log2 = self.transformLog2Ratio(ts.log2);
          ts.caller = self.activeCaller;
          if (self.equalDistance) {
            ts.start = self.getRatioIndex(ts.start);
            ts.end = self.getRatioIndex(ts.end);
          }
          return ts;
        }),
        function (d) {
          return `${d.caller}-${self.data.chromosome}-${d.start}-${d.end}`;
        }
      )
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
            .attr("stroke-width", 2)
            .attr("stroke-opacity", 0)
            .call((enter) => enter.transition().attr("stroke-opacity", 1)),
        (update) =>
          update.attr("stroke-opacity", 1).call((update) =>
            update
              .transition()
              .duration(this.animationDuration)
              .attr("stroke-opacity", 1)
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

  #plotBAF() {
    let bafData;
    if (this.equalDistance) {
      bafData = this.#data.baf
        .map((p) => {
          let di = { ...p };
          di.start = this.getRatioIndex(di.pos);
          di.end = di.start + 1;
          di.baf = this.transformBAF(di.baf);
          return di;
        })
        .filter(
          (p) =>
            p.start >= this.xScale.domain()[0] && p.start <= this.xScale.domain()[1]
        );
    } else {
      bafData = this.#data.baf
        .filter(
          (p) =>
            p.pos >= this.xScale.domain()[0] && p.pos <= this.xScale.domain()[1]
        )
        .map((p) => {
          let di = { ...p };
          di.baf = this.transformBAF(di.baf);
          return di;
        });
    }

    if (bafData.length > MAX_POINTS && !this.#showAllData) {
      bafData = slidingPixelWindowBAF(bafData, this.xScale, this.equalDistance ? "start" : "pos");
    }

    // Draw BAF on Canvas (scatter only)
    this.#ctx.save();
    this.#ctx.translate(this.margin.left, this.margin.top + this.plotHeight + this.margin.between);
    this.#ctx.fillStyle = "#333";
    this.#ctx.globalAlpha = 0.3;

    bafData.forEach((d) => {
      if (d.mean === undefined) {
        const x = this.xScale(this.equalDistance ? (d.start + d.end) / 2 : (d.pos !== undefined ? d.pos : (d.start + d.end) / 2));
        const y = this.bafYScale(d.baf);
        if (x >= 0 && x <= this.width - this.margin.left - this.margin.right) {
          this.#ctx.beginPath();
          this.#ctx.arc(x, y, 2, 0, 2 * Math.PI);
          this.#ctx.fill();
        }
      }
    });
    this.#ctx.restore();

    // Use SVG only for summarized BAF data
    const svgData = bafData.filter((d) => d.mean !== undefined);

    this.#bafArea
      .selectAll(".data-point")
      .data(svgData, (d) => {
        return `${d.pos}-${d.start}-${d.end}:${d.mean < 0.5 ? "-" : "+"}`;
      })
      .join(
        (enter) => {
          let g = enter.append("g").attr("class", "data-point").attr("opacity", 0);

          g.append("rect")
            .attr("class", "variance-rect")
            .attr("x", (d) => this.xScale(d.start))
            .attr("y", (d) => this.bafYScale(d.mean + d.sd))
            .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            .attr("height", (d) =>
              this.bafYScale(this.bafYScale.domain()[1] - 2 * d.sd)
            )
            .attr("fill", "#333")
            .attr("opacity", 0.3);

          g.append("line")
            .attr("class", "mean")
            .attr("x1", (d) => this.xScale(d.start))
            .attr("x2", (d) => this.xScale(d.end))
            .attr("y1", (d) => this.bafYScale(d.mean))
            .attr("y2", (d) => this.bafYScale(d.mean))
            .attr("stroke", "#333")
            .attr("opacity", 0.5);

          return g.transition().duration(this.animationDuration).attr("opacity", 1);
        },
        (update) => {
          update
            .selectAll(".variance-rect")
            .data((d) => [d])
            .transition()
            .duration(this.animationDuration)
            .attr("x", (d) => this.xScale(d.start))
            .attr("y", (d) => this.bafYScale(d.mean + d.sd))
            .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            .attr("height", (d) =>
              this.bafYScale(this.bafYScale.domain()[1] - 2 * d.sd)
            );

          update
            .selectAll(".mean")
            .data((d) => [d])
            .transition()
            .duration(this.animationDuration)
            .attr("x1", (d) => this.xScale(d.start))
            .attr("x2", (d) => this.xScale(d.end))
            .attr("y1", (d) => this.bafYScale(d.mean))
            .attr("y2", (d) => this.bafYScale(d.mean));

          return update;
        },
        (exit) => exit.transition().duration(this.animationDuration).attr("opacity", 0).remove()
      );
  }

  #plotAnnotations() {
    const annotData = this.#data.annotations.map(d => {
      let ad = { ...d };
      if (this.equalDistance) {
        ad.start = this.getRatioIndex(d.start);
        ad.end = this.getRatioIndex(d.end);
      }
      return ad;
    });

    this.annotations
      .selectAll(".annotation")
      .data(annotData, (d) => [d.name, d.start, d.end])
      .join(
        (enter) => {
          return enter
            .append("g")
            .attr("class", "annotation")
            .attr("clip-path", "url(#annotation-clip)")
            .attr("opacity", 0)
            .call((enter) =>
              enter
                .append("rect")
                .attr("class", "annotation-marker")
                .attr("x", (d) => this.xScale(d.start))
                .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
                .attr(
                  "height",
                  this.height - this.margin.top - this.margin.bottom
                )
                .attr("stroke", "#000")
                .attr("stroke-width", 0.5)
                .attr("fill", "#333")
                .attr("fill-opacity", 0.05)
                .attr("pointer-events", "none")
            )
            .call((enter) =>
              enter
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
                    this.plotHeight +
                    this.margin.between / 2 -
                    labelHeight / 2 -
                    2
                  );
                })
                .attr(
                  "width",
                  (d) => getTextDimensions(d.name, "0.8rem")[0] + 10
                )
                .attr(
                  "height",
                  (d) => getTextDimensions(d.name, "0.8rem")[1] + 4
                )
                .attr("fill", "#EEE")
                .attr("rx", 4)
            )
            .call((enter) =>
              enter
                .append("text")
                .attr("class", "annotation-label")
                .text((d) => d.name)
                .attr("x", (d) => this.xScale(d.start + (d.end - d.start) / 2))
                .attr("y", this.plotHeight + this.margin.between / 2)
                .attr("text-anchor", "middle")
                .attr("dominant-baseline", "central")
            )
            .call((enter) =>
              enter
                .transition()
                .duration(this.animationDuration)
                .attr("opacity", 1)
            );
        },
        (update) =>
          update
            .call((update) =>
              update
                .selectAll(".annotation-label")
                .transition()
                .duration(this.animationDuration)
                .attr("x", (d) => this.xScale(d.start + (d.end - d.start) / 2))
            )
            .call((update) =>
              update
                .selectAll(".annotation-label-background")
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
            )
            .call((update) =>
              update
                .selectAll(".annotation-marker")
                .transition()
                .duration(this.animationDuration)
                .attr("x", (d) => this.xScale(d.start))
                .attr("width", (d) => this.xScale(d.end) - this.xScale(d.start))
            )
            .call((update) =>
              update
                .transition()
                .duration(this.animationDuration)
                .attr("opacity", 1)
            ),
        (exit) => {
          return exit
            .transition()
            .duration(this.animationDuration)
            .attr("opacity", 0)
            .remove();
        }
      );
  }

  #setupCanvas() {
    const container = this.element.parentNode;
    this.#canvas = document.createElement("canvas");
    this.#canvas.className = "plot-canvas print-visible";
    container.appendChild(this.#canvas);
    this.#ctx = this.#canvas.getContext("2d");

    // Handle high-DPI displays
    const dpr = window.devicePixelRatio || 1;
    this.#canvas.width = this.width * dpr;
    this.#canvas.height = this.height * dpr;
    this.#canvas.style.width = "100%";
    this.#canvas.style.height = "auto";
    this.#ctx.scale(dpr, dpr);
  }

  #clearCanvas() {
    this.#ctx.clearRect(0, 0, this.width, this.height);
  }

  #drawAxes() {
    this.svg.selectAll(".y-axis").remove();
    this.svg
      .insert("g", "#lr-area-clip")
      .attr("transform", `translate(${this.margin.left},${this.margin.top})`)
      .attr("class", "y-axis primary-y-axis ratio-y-axis")
      .transition()
      .duration(this.animationDuration)
      .call(this.ratioYAxis);
    this.svg
      .insert("g", "#lr-area-clip")
      .attr(
        "transform",
        `translate(${this.width - this.margin.right},${this.margin.top})`
      )
      .attr("class", "y-axis secondary-y-axis cn-y-axis")
      .transition()
      .duration(this.animationDuration)
      .call(this.cnYAxis);
    this.svg
      .insert("g", "#lr-area-clip")
      .attr(
        "transform",
        `translate(${this.margin.left},${this.margin.top + this.plotHeight + this.margin.between
        })`
      )
      .attr("class", "y-axis primary-y-axis baf-y-axis")
      .transition()
      .duration(this.animationDuration)
      .call(this.bafYAxis);

    this.svg.select(".x-axis").remove();
    this.svg
      .insert("g", "#lr-area-clip")
      .attr(
        "transform",
        `translate(${this.margin.left},${this.height - this.margin.bottom})`
      )
      .attr("class", "x-axis")
      .transition()
      .duration(this.animationDuration)
      .call(this.xAxis);
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
        `translate(${this.width},${this.margin.top + this.plotHeight / 2
        }) rotate(90)`
      )
      .attr("class", "y-label")
      .text("Copy number")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");

    this.svg
      .append("text")
      .attr(
        "transform",
        `translate(0,${this.height - this.margin.bottom - this.plotHeight / 2
        }) rotate(-90)`
      )
      .attr("class", "y-label")
      .text("BAF")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "text-before-edge");
  }

  initialiseMousetrap() {
    let isDragging = false;
    const mouseTrap = this.svg
      .append("g")
      .attr("transform", `translate(${this.margin.left},${this.margin.top})`);
    mouseTrap
      .append("g")
      .attr("id", "lr-mousetrap")
      .append("rect")
      .attr("class", "mousetrap")
      .attr("width", this.width - this.margin.left - this.margin.right)
      .attr("height", this.plotHeight)
      .attr("fill", "none")
      .attr("pointer-events", "all");
    mouseTrap
      .append("g")
      .attr("transform", `translate(0,${this.plotHeight})`)
      .attr("id", "zoom-mousetrap")
      .append("rect")
      .attr("class", "mousetrap")
      .attr("width", this.width - this.margin.left - this.margin.right)
      .attr("height", this.margin.between)
      .attr("fill", "none")
      .attr("pointer-events", "all");
    mouseTrap
      .append("g")
      .attr(
        "transform",
        `translate(0,${this.plotHeight + this.margin.between})`
      )
      .attr("id", "baf-mousetrap")
      .append("rect")
      .attr("class", "mousetrap")
      .attr("width", this.width - this.margin.left - this.margin.right)
      .attr("height", this.plotHeight)
      .attr("fill", "none")
      .attr("pointer-events", "all");

    mouseTrap.select("#lr-mousetrap").on("mouseenter mousemove", (e) => {
      if (isDragging) {
        return;
      }
      let pos = d3.pointer(e);
      ratioCursor.set(pos[1]);
      positionCursor.set(pos[0]);
    });

    mouseTrap.select("#zoom-mousetrap").on("mouseenter mousemove", (e) => {
      let pos = d3.pointer(e);
      positionCursor.set(pos[0]);
    });

    mouseTrap.select("#baf-mousetrap").on("mouseenter mousemove", (e) => {
      if (isDragging) {
        return;
      }
      let pos = d3.pointer(e);
      bafCursor.set(pos[1]);
      positionCursor.set(pos[0]);
    });

    // Where the zooming elements are drawn
    const zoomLayer = mouseTrap.append("g").lower().attr("id", "zoom-layer");

    // Handle zooming
    mouseTrap
      .selectAll(".mousetrap")
      .on("mouseleave", () => {
        positionCursor.hide();
        ratioCursor.hide();
        bafCursor.hide();
      })
      .call(
        d3
          .drag()
          .on("start", (e) => {
            isDragging = true;
            zoomLayer
              .append("rect")
              .attr("id", "zoom-region")
              .attr("x", e.x)
              .attr(
                "height",
                this.height - this.margin.bottom - this.margin.top
              )
              .attr("stroke-width", 0)
              .attr("fill-opacity", 0.1)
              .attr("pointer-events", "none");
            ratioCursor.hide();
            bafCursor.hide();
          })
          .on("drag", (e) => {
            const leftBound = Math.min(e.x, e.subject.x);
            let width = Math.abs(Math.max(0, e.x) - e.subject.x);

            if (leftBound + width > this.xScale.range()[1]) {
              width = this.xScale.range()[1] - leftBound;
            }

            const genomeWidth =
              this.xScale.invert(Math.max(e.x, e.subject.x)) -
              this.xScale.invert(leftBound);

            const zoomRegion = zoomLayer
              .select("#zoom-region")
              .attr("x", Math.max(0, Math.min(e.x, e.subject.x)))
              .attr("width", width);

            if (genomeWidth < this.minZoomRange) {
              zoomRegion.attr("fill", "red");
            } else {
              zoomRegion.attr("fill", "#333");
            }

            positionCursor.set(e.x);
          })
          .on("end", (e) => {
            zoomLayer.select("#zoom-region").remove();
            const xMin = Math.max(0, Math.min(e.x, e.subject.x));
            const xMax = Math.min(
              this.xScale.range()[1],
              Math.max(e.x, e.subject.x)
            );
            isDragging = false;
            if (xMax - xMin < 3) {
              return;
            }
            this.zoomTo(this.xScale.invert(xMin), this.xScale.invert(xMax));
            this.update();
          })
      )
      .on("click", (e) => {
        isDragging = false;
        const [xMin, xMax] = this.xScale.domain();
        if (xMax - xMin !== this.length) {
          this.resetZoom();
          this.update();
        }
        const pos = d3.pointer(e);
        positionCursor.set(pos[0]);
        if (e.target.parentElement.id === "lr-mousetrap") {
          ratioCursor.set(pos[1]);
        } else if (e.target.parentElement.id == "baf-mousetrap") {
          bafCursor.set(pos[1]);
        }
      });

    const positionCursor = new XCursor({
      element: mouseTrap,
      height: this.plotHeight * 2 + this.margin.between,
      xScale: this.xScale,
    });

    const ratioCursor = new YCursor({
      element: mouseTrap.select("#lr-mousetrap"),
      width: this.width - this.margin.left - this.margin.right,
      yScale: this.ratioYScale,
      secondaryYScale: this.cnYScale,
    });

    const bafCursor = new YCursor({
      element: mouseTrap.select("#baf-mousetrap"),
      width: this.width - this.margin.left - this.margin.right,
      yScale: this.bafYScale,
    });
  }

  #updateAxes() {
    const [staticYMin, staticYMax] = [-2, 2];

    const [xMin, xMax] = this.zoomRange;
    let yMin, yMax;

    this.xScale.domain([0, this.length]);

    if (xMin > 0 || xMax < this.length) {
      this.xScale.domain([Math.max(xMin, 0), Math.min(xMax, this.length)]);

      let yValues;
      if (this.equalDistance) {
        yValues = this.#data.callers[this.#activeCaller].ratios.slice(
          Math.max(xMin, 0),
          Math.min(xMax, this.length)
        );
      } else {
        yValues = this.#data.callers[this.#activeCaller].ratios.filter(
          (d) => d.start > xMin && d.start < xMax
        );
      }

      if (yValues.length === 0) {
        yMin = staticYMin;
        yMax = staticYMax;
      } else {
        yMin = yValues
          .map((d) => this.transformLog2Ratio(d.log2))
          .reduce((a, d) => (d < a ? d : a));
        yMax = yValues
          .map((d) => this.transformLog2Ratio(d.log2))
          .reduce((a, d) => (d > a ? d : a));
      }
    } else {
      [yMin, yMax] = d3.extent(
        this.#data.callers[this.#activeCaller].ratios,
        (d) => this.transformLog2Ratio(d.log2)
      );
      if (!yMin && !yMax) {
        yMin = staticYMin;
        yMax = staticYMax;
      }
    }

    if (this.fitToData) {
      this.dispatchEvent(
        new CustomEvent("zoom", {
          detail: { dataOutsideRange: false },
        })
      );
      const padding = (yMax - yMin) * 0.05;
      this.ratioYScale.domain([yMin - padding, yMax + padding]);
      this.cnYScale.domain([
        cnFromRatio(yMin - padding),
        cnFromRatio(yMax + padding),
      ]);
    } else {
      this.dispatchEvent(
        new CustomEvent("zoom", {
          detail: { dataOutsideRange: yMin < staticYMin || yMax > staticYMax },
        })
      );
      const padding = (staticYMax - staticYMin) * 0.05;
      this.ratioYScale.domain([staticYMin - padding, staticYMax + padding]);
      this.cnYScale.domain([
        cnFromRatio(staticYMin - padding),
        cnFromRatio(staticYMax + padding),
      ]);
    }

    this.svg
      .transition()
      .select(".x-axis")
      .duration(this.animationDuration)
      .call(this.xAxis);
    this.svg
      .transition()
      .select(".ratio-y-axis")
      .duration(this.animationDuration)
      .call(this.ratioYAxis);
    this.svg
      .transition()
      .select(".cn-y-axis")
      .duration(this.animationDuration)
      .call(this.cnYAxis);

    this.svg.selectAll(".gridline").remove();
    this.svg
      .selectAll(".primary-y-axis .tick")
      .lower()
      .append("line")
      .attr("class", (d) => {
        return d == 0 ? "gridline baseline" : "gridline";
      })
      .attr("x2", this.xScale.range()[1]);

    this.svg.select(".x-label").text(this.#data.label);
  }

  getZoomRange() {
    return this.zoomRange;
  }

  zoomTo(start, end) {
    if (end - start < this.minZoomRange) {
      this.dispatchEvent(new CustomEvent("max-zoom-reached", {}));
      return this;
    }
    this.zoomRange = [start, end];
    this.xScale.domain(this.zoomRange);
  }

  resetZoom() {
    this.zoomRange = [0, this.length];
    this.xScale.domain(this.zoomRange);
  }

  setBaselineOffset(dy) {
    this.baselineOffset = dy;
    this.update();
  }

  update() {
    // Clear highlights when chromosome changes or plot updates
    this.#plotArea.selectAll(".highlight-region").remove();
    this.#bafArea.selectAll(".highlight-region").remove();
    this.#clearCanvas();

    this.#updateAxes();
    this.#plotRatios();
    this.#plotSegments();
    this.#plotBAF();
    this.#plotAnnotations();
    this.#plotCytobands();

    if (this.#highlightedRegion) {
      this.#drawHighlight(
        this.#highlightedRegion.start,
        this.#highlightedRegion.end,
        this.#highlightedRegion.name
      );
    }
    return this;
  }

  highlightRegion(start, end, name) {
    this.#highlightedRegion = { start, end, name };
    this.#drawHighlight(start, end, name);
  }

  #drawHighlight(start, end, name) {
    // Clear existing highlights first
    this.#plotArea.selectAll(".highlight-region").remove();
    this.#bafArea.selectAll(".highlight-region").remove();

    // Determine coordinates based on view mode (equal distance vs genomic)
    let xStart, xEnd;
    let s = start;
    let e = end;
    if (this.equalDistance) {
      s = this.getRatioIndex(s);
      e = this.getRatioIndex(e);
    }
    xStart = this.xScale(s);
    xEnd = this.xScale(e);

    // If region is too small, make it visible (min 2px)
    if (xEnd - xStart < 2) {
      const mid = (xStart + xEnd) / 2;
      xStart = mid - 1;
      xEnd = mid + 1;
    }

    const drawHighlight = (parentGroup, height) => {
      const g = parentGroup.append("g")
        .attr("class", "highlight-region");

      g.append("rect")
        .attr("x", xStart)
        .attr("y", 0)
        .attr("width", xEnd - xStart)
        .attr("height", height)
        .attr("fill", "red")
        .attr("stroke", "red")
        .attr("stroke-width", 1)
        .attr("opacity", 0)
        .transition()
        .duration(300)
        .attr("opacity", 0.2);

      return g;
    };

    // Draw in Main Plot (LR)
    drawHighlight(this.#plotArea, this.plotHeight);

    // Draw in BAF Plot
    drawHighlight(this.#bafArea, this.plotHeight);

    // Add Label - centered between plots with red text
    if (name) {
      const labelGroup = this.#plotArea.append("g")
        .attr("class", "highlight-region");

      const fontSize = "12px";
      const [labelWidth, labelHeight] = getTextDimensions(name, fontSize);
      const centerX = (xStart + xEnd) / 2;
      const centerY = this.plotHeight + this.margin.between / 2;

      // Draw box border (background)
      labelGroup.append("rect")
        .attr("x", centerX - labelWidth / 2 - 4)
        .attr("y", centerY - labelHeight / 2 - 2)
        .attr("width", labelWidth + 8)
        .attr("height", labelHeight + 4)
        .attr("fill", "white")
        .attr("stroke", "red")
        .attr("stroke-width", 1)
        .attr("rx", 4)
        .attr("opacity", 0)
        .transition()
        .duration(300)
        .attr("opacity", 1);

      // Draw text
      labelGroup.append("text")
        .attr("x", centerX)
        .attr("y", centerY)
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "central")
        .attr("fill", "red")
        .attr("font-size", fontSize)
        .attr("font-weight", "bold")
        .attr("opacity", 0)
        .text(name)
        .transition()
        .duration(300)
        .attr("opacity", 1);
    }
  }
}
