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

d3.select("#reports-version")
  .call((p) =>
    p
      .select("a")
      .attr(
        "href",
        `https://github.com/hydra-genetics/reports/tree/v${VERSION}`
      )
  )
  .select("span")
  .text(VERSION);

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

function* generateWindowSlices(points, scale, posAttr, windowSize = 5) {
  let offset = scale.domain()[0];
  let currentWindow = [];
  for (let p of points) {
    if (p[posAttr] < offset) {
      continue;
    }

    if (offset > scale.domain()[1]) {
      break;
    }

    if (p[posAttr] < offset + windowSize) {
      currentWindow.push(p);
    } else {
      yield currentWindow;
      offset += windowSize;
      currentWindow = [p];
    }
  }
}

function slidingPixelWindowVAF(points, scale, pixelWindowSize = 5) {
  let windowSize = Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]);
  if (windowSize < 4) {
    let reducedPoints = points.filter(
      (p) => p.pos >= scale.domain()[0] && p.pos < scale.domain()[1]
    );
    return reducedPoints;
  }

  let reducedPoints = [];

  for (let window of generateWindowSlices(points, scale, "pos", windowSize)) {
    let hist = Array(5).fill(0);
    window.forEach((p) => {
      let bin = Math.floor(5 * p.vaf);
      hist[bin] += 1;
    });
    let maxBin = hist.indexOf(Math.max(...hist));
    if (maxBin !== 2) {
      // potentially bimodal
      let hiPoints = window.filter((p) => p.vaf >= 0.5);
      if (hiPoints.length > 0) {
        let hiSum = hiPoints.reduce((a, b) => a + b.vaf, 0);
        let hiMean = hiSum / hiPoints.length;
        let hiSD = Math.sqrt(
          hiPoints.reduce((a, b) => a + Math.pow(b.vaf - hiMean, 2), 0) /
            (hiPoints.length - 1)
        );
        reducedPoints.push({
          start: hiPoints[0].pos,
          end: hiPoints[hiPoints.length - 1].pos,
          mean: hiMean,
          sd: hiPoints.length > 1 ? hiSD : 0,
        });
      }

      let loPoints = window.filter((p) => p.vaf < 0.5);
      if (loPoints.length > 0) {
        let loSum = loPoints.reduce((a, b) => a + b.vaf, 0);
        let loMean = loSum / loPoints.length;
        let loSD = Math.sqrt(
          loPoints.reduce((a, b) => a + Math.pow(b.vaf - loMean, 2), 0) /
            (loPoints.length - 1)
        );

        reducedPoints.push({
          start: loPoints[0].pos,
          end: loPoints[loPoints.length - 1].pos,
          mean: loMean,
          sd: loPoints.length > 1 ? loSD : 0,
        });
      }
    } else {
      let winSum = window.reduce((a, b) => a + b.vaf, 0);
      let winMean = winSum / window.length;
      let winSD = Math.sqrt(
        window.reduce((a, b) => a + Math.pow(b.vaf - winMean, 2), 0) /
          (window.length - 1)
      );
      reducedPoints.push({
        start: window[0].pos,
        end: window[window.length - 1].pos,
        mean: winMean,
        sd: window.length > 1 ? winSD : 0,
      });
    }
  }

  return reducedPoints;
}

function slidingPixelWindow(
  points,
  scale,
  posAttr,
  valAttr,
  pixelWindowSize = 5
) {
  let windowSize = Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]);
  if (windowSize < 4) {
    let reducedPoints = points.filter(
      (p) => p[posAttr] >= scale.domain()[0] && p[posAttr] < scale.domain()[1]
    );
    return reducedPoints;
  }

  let reducedPoints = [];

  for (let window of generateWindowSlices(points, scale, posAttr, windowSize)) {
    let winSum = window.reduce((a, b) => a + b[valAttr], 0);
    let winMean = winSum / window.length;
    let winSD = Math.sqrt(
      window.reduce((a, b) => a + Math.pow(b[valAttr] - winMean, 2), 0) /
        (window.length - 1)
    );
    reducedPoints.push({
      start: window[0][posAttr],
      end: window[window.length - 1][posAttr],
      mean: winMean,
      sd: winSD,
    });
  }

  return reducedPoints;
}

const chromosomePlot = new ChromosomePlot({
  element: document.querySelector("#chromosome-view"),
  data: cnvData[0],
});

const genomePlot = new GenomePlot({
  element: document.querySelector("#genome-view"),
  data: cnvData,
});

const resultsTable = new ResultsTable(d3.select("#cnv-table"), {
  data: cnvData,
  filter: d3.select("#table-filter-toggle").node().checked,
});

const messageModal = document.querySelector("dialog");
messageModal.addEventListener("click", (e) => {
  const modalDimensions = messageModal.getBoundingClientRect();
  if (
    e.clientX < modalDimensions.left ||
    e.clientX > modalDimensions.right ||
    e.clientY < modalDimensions.top ||
    e.clientY > modalDimensions.bottom
  ) {
    messageModal.close();
  }
});

document.querySelector("dialog button.close").addEventListener("click", (e) => {
  e.currentTarget.parentNode.close();
});

function setModalMessage(msg, className) {
  let message = document.createElement("p");
  const messageText = document.createTextNode(msg);

  let icon = document.createElement("i");

  if (className === "error") {
    icon.className = "fa-solid fa-circle-exclamation";
  } else if (className === "warning") {
    icon.className = "fa-solid fa-triangle-exclamation";
  } else if (className === "info") {
    icon.className = "fa-solid fa-circle-info";
  }

  message.appendChild(icon);
  message.appendChild(messageText);

  messageModal.className = className ? className : "";
  messageModal.firstChild?.remove();
  messageModal.prepend(message);
}

chromosomePlot.addEventListener("zoom", (e) => {
  d3.selectAll(".data-range-warning").classed(
    "hidden",
    !e.detail.dataOutsideRange
  );
});

chromosomePlot.addEventListener("max-zoom-reached", () => {
  setModalMessage(
    "Trying to zoom in too far. " +
      `Current lower limit is ${chromosomePlot.minZoomRange} bp.`,
    "error"
  );
  messageModal.showModal();
});

genomePlot.addEventListener("chromosome-change", (e) => {
  chromosomePlot.setData(
    cnvData[e.detail.chromosome],
    e.detail.start,
    e.detail.end
  );
});

genomePlot.addEventListener("chromosome-zoom", (e) => {
  chromosomePlot.setData(null, e.detail.start, e.detail.end);
});

resultsTable.addEventListener("zoom-to-region", (e) => {
  genomePlot.selectChromosome(
    e.detail.chromosome,
    e.detail.start,
    e.detail.start + e.detail.length
  );
});

d3.select("#table-filter-toggle").on("change", (event) => {
  resultsTable.filter = event.target.checked;
});

d3.select("#chromosome-fit-to-data").on("change", (e) => {
  chromosomePlot.fitToData = e.target.checked;
});

d3.selectAll("input[name=dataset]").on("change", (e) => {
  chromosomePlot.activeCaller = parseInt(e.target.value);
  genomePlot.activeCaller = parseInt(e.target.value);
  resultsTable.activeCaller = parseInt(e.target.value);
});
