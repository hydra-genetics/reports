const MAX_POINTS = 300;
const MIN_COPY_NUMBER = 0.001;
const MIN_LOG2_RATIO = Math.log2(MIN_COPY_NUMBER / 2);

function* generateWindowSlices(points, scale, posAttr, windowSize = 5) {
  let offset = scale.domain()[0];
  let currentWindow = [];
  for (let p of points) {
    if (p[posAttr] < offset) continue;
    if (offset > scale.domain()[1]) break;

    if (p[posAttr] < offset + windowSize) {
      currentWindow.push(p);
    } else {
      while (p[posAttr] >= offset + windowSize) {
        yield { window: currentWindow, start: offset };
        offset += windowSize;
        currentWindow = [];
      }
      currentWindow = [p];
    }
  }
  if (currentWindow.length > 0) {
    yield { window: currentWindow, start: offset };
  }
}

function summariseWindow(points, windowStart, windowSize, valAttr, minValue = undefined, offset = 0) {
  if (!points || points.length === 0) return null;
  let hasOutliers = false;
  let filtered = points;
  if (minValue !== undefined) {
    hasOutliers = points.some((x) => (x[valAttr] ?? x.mean ?? 0) + offset <= minValue);
    filtered = points.filter((x) => (x[valAttr] ?? x.mean ?? 0) + offset > minValue);
  }
  if (filtered.length === 0) return null;

  let sum = filtered.reduce((a, b) => a + (b[valAttr] ?? b.mean ?? 0), 0);
  let mean = sum / filtered.length;

  // Calculate pooled variance: Between-group variance + Within-group variance
  // Between-group: Variance of the means
  let betweenVar = filtered.reduce((a, b) => a + Math.pow((b[valAttr] ?? b.mean ?? 0) - mean, 2), 0) / filtered.length;
  
  // Within-group: Average of the variances (sd^2) of the input points (if available)
  // If input is raw, sd is undefined/0, so withinVar contributes 0.
  let withinVar = filtered.reduce((a, b) => a + Math.pow(b.sd ?? 0, 2), 0) / filtered.length;

  let sd = Math.sqrt(betweenVar + withinVar);
  
  return {
    start: windowStart,
    end: windowStart + windowSize,
    pos: windowStart + windowSize / 2,
    baf: mean,
    mean: mean,
    sd: sd,
    hasOutliers: hasOutliers
  };
}

function slidingPixelWindowBAF(points, scale, posAttr = "pos", pixelWindowSize = 5, force = false) {
  const [d0, d1] = scale.domain();
  points = points.filter(p => (p.end ?? p[posAttr]) >= d0 && (p.start ?? p[posAttr]) < d1);
  if (points.length === 0) return [];

  // Sort by position to ensure correct binning
  points.sort((a, b) => (a[posAttr] ?? a.start) - (b[posAttr] ?? b.start));

  if (points[0]?.baf instanceof Array) {
    points = points.flatMap(p => p.baf.map(v => ({ ...p, baf: v })));
  }
  
  let windowSize = Math.max(1, Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]));
  if (!force && (windowSize < 4 || points.length <= MAX_POINTS)) return points;

  let reducedPoints = [];
  for (let slice of generateWindowSlices(points, scale, posAttr, windowSize)) {
    let { window, start } = slice;
    if (window.length === 0) continue;

    let hist = Array(5).fill(0);
    window.forEach(p => {
      let b = p.baf ?? p.mean;
      if (b !== undefined) hist[Math.min(4, Math.floor(b * 5))] += 1;
    });

    if (hist.indexOf(Math.max(...hist)) !== 2) {
      let hi = window.filter(p => (p.baf ?? p.mean) >= 0.5);
      let lo = window.filter(p => (p.baf ?? p.mean) < 0.5);
      let sHi = summariseWindow(hi, start, windowSize, "baf");
      let sLo = summariseWindow(lo, start, windowSize, "baf");
      if (sHi) reducedPoints.push(sHi);
      if (sLo) reducedPoints.push(sLo);
    } else {
      let s = summariseWindow(window, start, windowSize, "baf");
      if (s) reducedPoints.push(s);
    }
  }
  return reducedPoints;
}



function slidingPixelWindow(
  points,
  scale,
  posAttr,
  valAttr,
  offset = 0,
  pixelWindowSize = 5,
  force = false
) {
  const [d0, d1] = scale.domain();
  points = points.filter(p => (p.end ?? p[posAttr]) >= d0 && (p.start ?? p[posAttr]) < d1);
  if (points.length === 0) return [];

  let windowSize = Math.max(1, Math.ceil(scale.invert(pixelWindowSize) - scale.domain()[0]));
  if (!force && (windowSize < 4 || points.length <= MAX_POINTS)) {
    return points;
  }

  let reducedPoints = [];

  for (let slice of generateWindowSlices(points, scale, posAttr, windowSize)) {
    let window = slice.window;
    let windowStart = slice.start;
    if (window.length === 0) {
      continue;
    }
    let s = summariseWindow(
      window,
      windowStart,
      windowSize,
      valAttr,
      MIN_LOG2_RATIO,
      offset
    );
    if (s) reducedPoints.push(s);
  }

  return reducedPoints;
}

