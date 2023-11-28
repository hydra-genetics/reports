const MAX_POINTS = 300;

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
      // Fast forward to the next window with data
      while (p[posAttr] >= offset + windowSize) {
        yield currentWindow;
        offset += windowSize;
        currentWindow = [];
      }
      currentWindow = [p];
    }
  }
}

function summariseWindow(points, windowStart, windowSize, valAttr) {
  let sum = points.reduce((a, b) => a + b[valAttr], 0);
  let mean = sum / points.length;
  let sd = Math.sqrt(
    points.reduce((a, b) => a + Math.pow(b[valAttr] - mean, 2), 0) /
      points.length
  );
  return {
    start: windowStart,
    end: windowStart + windowSize,
    mean: mean,
    sd: sd,
  };
}

function slidingPixelWindowVAF(
  points,
  scale,
  pixelWindowSize = 5,
  force = false
) {
  points = points
    .filter((p) => p.pos >= scale.domain()[0] && p.pos < scale.domain()[1])
    .map((p) =>
      p.vaf.map((v) => {
        return { pos: p.pos, vaf: v };
      })
    )
    .flat();
  let windowSize = Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]);
  let windowStart = Math.floor(scale.domain()[0]);

  if (!force && (windowSize < 4 || points.length <= MAX_POINTS)) {
    return points;
  }

  let reducedPoints = [];

  for (let window of generateWindowSlices(points, scale, "pos", windowSize)) {
    if (window.length === 0) {
      windowStart += windowSize;
      continue;
    }
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
        reducedPoints.push(
          summariseWindow(hiPoints, windowStart, windowSize, "vaf")
        );
      }

      let loPoints = window.filter((p) => p.vaf < 0.5);
      if (loPoints.length > 0) {
        reducedPoints.push(
          summariseWindow(loPoints, windowStart, windowSize, "vaf")
        );
      }
    } else {
      reducedPoints.push(
        summariseWindow(window, windowStart, windowSize, "vaf")
      );
    }
    windowStart += windowSize;
  }

  return reducedPoints;
}

function slidingPixelWindow(
  points,
  scale,
  posAttr,
  valAttr,
  pixelWindowSize = 5,
  force = false
) {
  points = points.filter(
    (p) => p[posAttr] >= scale.domain()[0] && p[posAttr] < scale.domain()[1]
  );
  let windowSize = Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]);
  let windowStart = Math.floor(scale.domain()[0]);
  if (!force && (windowSize < 4 || points.length <= MAX_POINTS)) {
    return points;
  }

  let reducedPoints = [];

  for (let window of generateWindowSlices(points, scale, posAttr, windowSize)) {
    if (window.length === 0) {
      windowStart += windowSize + 1;
      continue;
    }
    reducedPoints.push(
      summariseWindow(window, windowStart, windowSize, valAttr)
    );
    windowStart += windowSize + 1;
  }

  return reducedPoints;
}
