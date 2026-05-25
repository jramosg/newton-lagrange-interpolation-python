const functions = {
  runge: (x) => 1 / (1 + 25 * x * x),
  sin: (x) => Math.sin(x),
  gaussian: (x) => Math.exp(-20 * x * x),
};

const controls = {
  functionSelect: document.getElementById("functionSelect"),
  nodeSelect: document.getElementById("nodeSelect"),
  countInput: document.getElementById("countInput"),
  pointInput: document.getElementById("pointInput"),
  demoResult: document.getElementById("demoResult"),
};

function equispacedNodes(n) {
  return Array.from({ length: n }, (_, i) => -1 + (2 * i) / (n - 1));
}

function chebyshevNodes(n) {
  return Array.from({ length: n }, (_, i) => {
    const k = i + 1;
    return Math.cos(((2 * k - 1) * Math.PI) / (2 * n));
  }).sort((a, b) => a - b);
}

function dividedDifferences(xs, ys) {
  const n = xs.length;
  const table = Array.from({ length: n }, () => Array(n).fill(0));
  ys.forEach((y, i) => {
    table[i][0] = y;
  });
  for (let order = 1; order < n; order += 1) {
    for (let row = 0; row < n - order; row += 1) {
      const numerator = table[row + 1][order - 1] - table[row][order - 1];
      const denominator = xs[row + order] - xs[row];
      table[row][order] = numerator / denominator;
    }
  }
  return table[0];
}

function evaluateNewton(coefficients, xs, x) {
  let result = coefficients[coefficients.length - 1];
  for (let i = coefficients.length - 2; i >= 0; i -= 1) {
    result = result * (x - xs[i]) + coefficients[i];
  }
  return result;
}

function getSeries(settings) {
  const fn = functions[settings.functionName];
  const xs = settings.nodeType === "chebyshev"
    ? chebyshevNodes(settings.n)
    : equispacedNodes(settings.n);
  const ys = xs.map(fn);
  const coefficients = dividedDifferences(xs, ys);
  const real = [];
  const interp = [];
  for (let i = 0; i < 360; i += 1) {
    const x = -1 + (2 * i) / 359;
    real.push([x, fn(x)]);
    interp.push([x, evaluateNewton(coefficients, xs, x)]);
  }
  return { xs, ys, coefficients, real, interp };
}

function drawChart(canvas, settings, options = {}) {
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  const padding = { top: 24, right: 28, bottom: 38, left: 48 };
  const series = getSeries(settings);
  const yValues = [
    ...series.real.map((p) => p[1]),
    ...series.interp.map((p) => p[1]),
    ...series.ys,
  ];
  const yMin = Math.min(-1.2, ...yValues);
  const yMax = Math.max(1.2, ...yValues);
  const xToPx = (x) =>
    padding.left +
    ((x + 1) / 2) * (width - padding.left - padding.right);
  const yToPx = (y) =>
    padding.top +
    ((yMax - y) / (yMax - yMin)) *
      (height - padding.top - padding.bottom);

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  ctx.strokeStyle = "#d7dde5";
  ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i += 1) {
    const x = padding.left + i * (width - padding.left - padding.right) / 4;
    ctx.beginPath();
    ctx.moveTo(x, padding.top);
    ctx.lineTo(x, height - padding.bottom);
    ctx.stroke();
  }
  for (let i = 0; i <= 4; i += 1) {
    const y = padding.top + i * (height - padding.top - padding.bottom) / 4;
    ctx.beginPath();
    ctx.moveTo(padding.left, y);
    ctx.lineTo(width - padding.right, y);
    ctx.stroke();
  }

  drawLine(ctx, series.real, xToPx, yToPx, "#17202a", 2);
  drawLine(ctx, series.interp, xToPx, yToPx, "#1e5eff", 2.5);

  ctx.fillStyle = "#c0392b";
  series.xs.forEach((x, i) => {
    ctx.beginPath();
    ctx.arc(xToPx(x), yToPx(series.ys[i]), 4.2, 0, Math.PI * 2);
    ctx.fill();
  });

  if (options.evaluationPoint !== undefined) {
    const x = options.evaluationPoint;
    const y = evaluateNewton(series.coefficients, series.xs, x);
    ctx.fillStyle = "#087f5b";
    ctx.beginPath();
    ctx.arc(xToPx(x), yToPx(y), 5.5, 0, Math.PI * 2);
    ctx.fill();
  }

  ctx.fillStyle = "#5d6673";
  ctx.font = "13px system-ui, sans-serif";
  ctx.fillText("-1", padding.left - 6, height - 14);
  ctx.fillText("1", width - padding.right - 4, height - 14);
  ctx.fillText("x", width / 2, height - 14);
}

function drawLine(ctx, points, xToPx, yToPx, color, lineWidth) {
  ctx.strokeStyle = color;
  ctx.lineWidth = lineWidth;
  ctx.beginPath();
  points.forEach(([x, y], index) => {
    if (index === 0) {
      ctx.moveTo(xToPx(x), yToPx(y));
    } else {
      ctx.lineTo(xToPx(x), yToPx(y));
    }
  });
  ctx.stroke();
}

function currentSettings() {
  return {
    functionName: controls.functionSelect.value,
    nodeType: controls.nodeSelect.value,
    n: Math.max(3, Math.min(25, Number(controls.countInput.value) || 11)),
  };
}

function updateDemo() {
  const settings = currentSettings();
  const x = Math.max(-1, Math.min(1, Number(controls.pointInput.value)));
  const series = getSeries(settings);
  const y = evaluateNewton(series.coefficients, series.xs, x);
  controls.demoResult.textContent = `P(${x.toFixed(2)}) = ${y.toFixed(8)}`;
  drawChart(document.getElementById("demoChart"), settings, {
    evaluationPoint: x,
  });
}

["change", "input"].forEach((eventName) => {
  Object.values(controls).forEach((control) => {
    if (control instanceof HTMLElement && control.tagName !== "DIV") {
      control.addEventListener(eventName, updateDemo);
    }
  });
});

drawChart(document.getElementById("heroChart"), {
  functionName: "runge",
  nodeType: "equispaced",
  n: 11,
});
updateDemo();
