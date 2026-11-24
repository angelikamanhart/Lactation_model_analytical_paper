import React, { useState, useMemo, useRef, useEffect, useLayoutEffect, useCallback } from "react";
import { createRoot } from "react-dom/client";


/* ------------------------------------------------------------------ tokens */
const C = {
  paper: "#E7ECE9",
  panel: "#FBFCFB",
  ink: "#13211E",
  ink2: "#53645F",
  ink3: "#87968F",
  rule: "#C3CFCA",
  grid: "#D9E2DD",
  invol: "#F0F2F2", // involution wins
  lact: "#C6E8AF", // bistable, lactation state
  empt: "#F6DBAA", // bistable, emptied state
  unstTT: "#42A860", // unstable target-met fixed point, target-met / target-met attractor
  unstET: "#EDC28D", // unstable target-met fixed point, emptying / target-met attractor
  sep: "#281F8C", // separatrix (blue in the manuscript)
  flip: "#992923", // stability-loss boundary (red in Fig. 8)
  traj: "#0E6B6B", // trajectory
  mark: "#D9930F", // selection
  opt: "#6F4C9B", // numerical M_DAY optimum trace
};
const MONO = "'IBM Plex Mono', ui-monospace, SFMono-Regular, Menlo, monospace";
const SANS = "'IBM Plex Sans', ui-sans-serif, system-ui, sans-serif";
const SERIF = "'IBM Plex Serif', Georgia, serif";

/* Active stimulation law (notation and parameters as in the manuscript). */
function ResponseFormula({P, mode="amount"}) {
  const input = mode === "amount" ? "M" : "τ";
  const suffix = mode === "amount" ? "M" : "τ";
  return (
    <div style={{fontFamily:SERIF,fontSize:12,color:C.ink2,lineHeight:1.55,
      background:C.paper,border:`1px solid ${C.rule}`,padding:"6px 8px",marginBottom:12}}>
      <span style={{fontFamily:SANS,fontSize:10,color:C.ink3,display:"block",marginBottom:2}}>
        {P.resp === "sig" ? "Hill response" : "Michaelis–Menten response"}
      </span>
      <i>a</i>({input}) = <i>a</i><sub>max</sub> · {input}{P.resp === "sig" ? <sup>2</sup> : null}
      {" / ("}{input}{P.resp === "sig" ? <sup>2</sup> : null}{" + μ"}<sub>{suffix}</sub>{P.resp === "sig" ? <sup>2</sup> : null}{")"}
    </div>
  );
}

/* -------------------------------------------------------------------- math */
/* Stimulation-efficiency maximum: Michaelis–Menten versus Hill. */
function amountRmax(P) {
  return P.amax / (P.resp === "sig" ? 2 * P.muM : P.muM);
}

function model(P, T, Mt) {
  const g = Math.exp(-P.gam * T);
  const d = Math.exp(-P.dlt * T);
  const G =
    Math.abs(P.dlt - P.gam) < 1e-9
      ? T * Math.exp(-P.gam * T)
      : (Math.exp(-P.gam * T) - Math.exp(-P.dlt * T)) / (P.dlt - P.gam);
  const Gamma = (1 - g) / G;
  const kappa = amountRmax(P) / Gamma;
  const c = P.muM;
  const sig = P.resp === "sig";
  const a = sig
    ? (M) => (P.amax * M * M) / (M * M + c * c)
    : (M) => (P.amax * M) / (M + P.muM);

  let psep = null, pemp = null, Msep = null, Memp = null;
  if (kappa > 1) {
    if (sig) {
      const disc = Math.sqrt(Math.max((P.amax * G) ** 2 - (2 * P.muM * (1 - g)) ** 2, 0));
      psep = (P.amax * G - disc) / (2 * G * (1 - g));
      pemp = (P.amax * G + disc) / (2 * G * (1 - g));
      Msep = psep * G;
      Memp = pemp * G;
    } else {
      pemp = P.amax / (1 - g) - P.muM / G;
      Memp = pemp * G;
    }
  }

  const lacExists = Mt > 0 && a(Mt) / Mt > Gamma;
  const plac = lacExists ? a(Mt) / (1 - g) : null;
  const mlac = lacExists ? (plac * G - Mt) / (1 - d) : null;
  const emptied = pemp !== null && Mt >= Memp;

  const involStable = sig || kappa <= 1;
  const bistable = sig && (lacExists || emptied);
  const sInvariant = !lacExists; // {m=0} forward invariant  <=>  a(Mt) <= Gamma*Mt

  const step = (m, p) => {
    const x = m * d + p * G;
    return [x > Mt ? x - Mt : 0, p * g + a(x < Mt ? x : Mt)];
  };

  return { P, T, Mt, g, d, G, Gamma, kappa, a, sig, psep, pemp, Msep, Memp,
    lacExists, plac, mlac, emptied, involStable,
    bistable, sInvariant, step };
}

/* separatrix branches: arc inside E, then affine pre-images outside E */
function separatrix(M) {
  if (!M.sig || M.psep == null || M.Mt < M.Msep) return [];
  const { g, d, G, Mt, a, psep, P } = M;
  const xInt = psep >= P.amax ? Infinity : P.muM * Math.sqrt(psep / (P.amax - psep));
  const xEnd = Math.min(Mt, xInt);
  const x0 = psep * G;
  if (!(xEnd > x0)) return [];
  const b1 = [];
  const n = 220;
  for (let i = 0; i <= n; i++) {
    const x = x0 + ((xEnd - x0) * i) / n;
    const p = (psep - a(x)) / g;
    const m = (x - (G / g) * (psep - a(x))) / d;
    b1.push([Math.max(m, 0), Math.max(p, 0)]);
  }
  const branches = [b1];
  const aM = a(Mt);
  let cur = b1;
  for (let k = 0; k < 14; k++) {
    const nxt = [];
    for (const [m1, p1] of cur) {
      if (p1 >= aM) {
        const p = (p1 - aM) / g;
        nxt.push([(m1 + Mt - p * G) / d, p]);
      }
    }
    if (nxt.length < 2) break;
    branches.push(nxt);
    cur = nxt;
  }
  return branches;
}

/* which attractor: 0 = involution, 1 = the other one */
function classify(m, p, M) {
  const { d, G, Mt, a, g, psep, sInvariant } = M;
  for (let i = 0; i < 400; i++) {
    const x = m * d + p * G;
    const nm = x > Mt ? x - Mt : 0;
    p = p * g + a(x < Mt ? x : Mt);
    m = nm;
    if (!isFinite(m) || m > 1e14) return 1;
    if (m === 0 && sInvariant) return p < psep ? 0 : 1; // exact: orbit is on S
  }
  return p < psep ? 0 : 1;
}

/* Critical wait interval for amount-dependent stimulation: Γ(Tcrit) = Rmax. */
function criticalWait(P) {
  const rmax = amountRmax(P);
  if (rmax <= P.gam) return null;
  let lo = 1e-7, hi = 1;
  while (hi < 1e5 && waitPars(P, hi).Gamma < rmax) hi *= 2;
  if (hi >= 1e5) return null;
  for (let i = 0; i < 65; i++) {
    const mid = (lo + hi) / 2;
    if (waitPars(P, mid).Gamma < rmax) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}
function markCriticalWait(ctx, w, h, P, Tmax, X) {
  const tcrit = criticalWait(P);
  if (tcrit == null) return;
  ctx.save();
  ctx.font = `10px ${MONO}`;
  if (tcrit <= Tmax) {
    const xx = X(tcrit);
    ctx.strokeStyle = C.ink2; ctx.lineWidth = 1.15; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(xx, PAD.t); ctx.lineTo(xx, h - PAD.b); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = C.ink2;
    ctx.textAlign = xx > w * 0.7 ? "right" : "left";
    ctx.fillText(`T_crit = ${fmt3(tcrit)} h`, xx + (xx > w * 0.7 ? -5 : 5), h - PAD.b - 7);
  } else {
    ctx.fillStyle = C.ink3; ctx.textAlign = "right";
    ctx.fillText(`T_crit = ${fmt3(tcrit)} h (outside range)`, w - PAD.r - 4, h - PAD.b - 8);
  }
  ctx.restore();
}

/* --------------------------------------------------------------- debounce */
/* Heavy grids (basins of attraction, region and output maps, numerical
   bifurcation diagrams) are recomputed only once the inputs have been still
   for a moment, so that sliders stay responsive while being dragged. The
   value passed in must keep its identity between renders (state or useMemo). */
function useDebounced(value, delay = 180) {
  const [v, setV] = useState(value);
  useEffect(() => {
    if (Object.is(v, value)) return;
    const id = setTimeout(() => setV(value), delay);
    return () => clearTimeout(id);
  }, [value, delay]);
  return v;
}

/* ------------------------------------------------------------------ canvas */
function useCanvas(height, draw, onPick) {
  const wrap = useRef(null);
  const cvs = useRef(null);
  const [w, setW] = useState(420);
  // Measure before paint on initial mount and again whenever a collapsed
  // panel becomes visible. ResizeObserver keeps the canvas exactly as wide
  // as its containing panel after browser/layout changes.
  useLayoutEffect(() => {
    const el = wrap.current;
    if (!el) return;
    const measure = () => {
      const width = el.getBoundingClientRect().width;
      if (width > 0) setW(Math.max(260, Math.floor(width)));
    };
    measure();
    const ro = new ResizeObserver(measure);
    ro.observe(el);
    window.addEventListener("resize", measure);
    return () => { ro.disconnect(); window.removeEventListener("resize", measure); };
  }, []);
  useEffect(() => {
    const el = cvs.current;
    if (!el) return;
    const dpr = window.devicePixelRatio || 1;
    el.width = w * dpr;
    el.height = height * dpr;
    el.style.width = w + "px";
    el.style.height = height + "px";
    const ctx = el.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, height);
    draw(ctx, w, height);
  }, [w, height, draw]);
  const click = (e) => {
    if (!onPick) return;
    const r = cvs.current.getBoundingClientRect();
    onPick(e.clientX - r.left, e.clientY - r.top, w, height);
  };
  return { wrap, cvs, w, click };
}

const PAD = { l: 48, r: 14, t: 14, b: 34 };

function frame(ctx, w, h, xlab, ylab, xmax, ymax, xmin = 0, ymin = 0) {
  const X = (v) => PAD.l + ((v - xmin) / (xmax - xmin)) * (w - PAD.l - PAD.r);
  const Y = (v) => h - PAD.b - ((v - ymin) / (ymax - ymin)) * (h - PAD.t - PAD.b);
  ctx.strokeStyle = C.grid;
  ctx.lineWidth = 1;
  ctx.fillStyle = C.ink3;
  ctx.font = `10px ${MONO}`;
  const fx = tickFormatter((xmax - xmin) / 4);
  const fy = tickFormatter((ymax - ymin) / 4);
  for (let i = 0; i <= 4; i++) {
    const xv = xmin + ((xmax - xmin) * i) / 4;
    const yv = ymin + ((ymax - ymin) * i) / 4;
    ctx.beginPath(); ctx.moveTo(X(xv), PAD.t); ctx.lineTo(X(xv), h - PAD.b); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(PAD.l, Y(yv)); ctx.lineTo(w - PAD.r, Y(yv)); ctx.stroke();
    ctx.textAlign = "center";
    ctx.fillText(fx(xv), X(xv), h - PAD.b + 14);
    ctx.textAlign = "right";
    ctx.fillText(fy(yv), PAD.l - 6, Y(yv) + 3);
  }
  ctx.strokeStyle = C.rule;
  ctx.strokeRect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
  ctx.fillStyle = C.ink2;
  ctx.font = `11px ${SANS}`;
  ctx.textAlign = "center";
  ctx.fillText(xlab, PAD.l + (w - PAD.l - PAD.r) / 2, h - 4);
  ctx.save();
  ctx.translate(11, PAD.t + (h - PAD.t - PAD.b) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(ylab, 0, 0);
  ctx.restore();
  return { X, Y };
}
/* Axis ticks: the number of decimals follows the tick spacing, so that
   neighbouring ticks never print the same label (e.g. on a 0–0.2 axis). */
function tickFormatter(step) {
  if (!(step > 0) || !isFinite(step)) return (v) => fmt(v);
  let dec = Math.max(0, Math.ceil(-Math.log10(step) - 1e-9));
  const scaled = step * Math.pow(10, dec);
  if (Math.abs(scaled - Math.round(scaled)) > 1e-6 * scaled) dec += 1;
  dec = Math.min(dec, 4);
  return (v) => (Math.abs(v) < 0.5 * Math.pow(10, -dec) ? "0" : v.toFixed(dec));
}
/* General values: integers above 10, two significant digits below. */
const fmt = (v) => {
  if (!isFinite(v)) return String(v);
  const a = Math.abs(v);
  if (a >= 10) return v.toFixed(0);
  if (a === 0) return "0";
  return String(parseFloat(v.toPrecision(2)));
};
const fmt3 = (v) =>
  v == null ? "—" : !isFinite(v) ? "∞" : Math.abs(v) >= 1000 ? v.toExponential(2) : v.toFixed(3);


/* ------------------------------------------------------- daily output + export */
// Rainbow/jet colormap used for daily-output heatmaps in Figure 11 of the paper.
// Low positive outputs are dark blue, high outputs are red.
const HEAT_STOPS = [
  [0.000, "#000080"],
  [0.125, "#0000FF"],
  [0.375, "#00FFFF"],
  [0.625, "#FFFF00"],
  [0.875, "#FF0000"],
  [1.000, "#800000"],
];
const clamp01 = (x) => Math.max(0, Math.min(1, x));
function mixHex(a, b, t) {
  const A = hex(a), B = hex(b), q = clamp01(t);
  const v = A.map((x, i) => Math.round(x + (B[i] - x) * q));
  return `rgb(${v[0]},${v[1]},${v[2]})`;
}
function heatNorm(v, vmax, scale = "linear") {
  if (!(vmax > 0)) return 0;
  const q = Math.max(0, Math.min(v, vmax));
  if (scale === "log") return clamp01(Math.log1p(q) / Math.log1p(vmax));
  return clamp01(q / vmax);
}
function heatValueAtNorm(u, vmax, scale = "linear") {
  const q = clamp01(u);
  if (!(vmax > 0)) return 0;
  return scale === "log" ? Math.expm1(q * Math.log1p(vmax)) : q * vmax;
}
function heatColor(v, vmax, scale = "linear") {
  if (!(v > 0)) return C.invol;  // zero output: grey, as in the paper
  const u = heatNorm(v, vmax, scale);
  for (let i = 0; i < HEAT_STOPS.length - 1; i++) {
    const [u0, c0] = HEAT_STOPS[i], [u1, c1] = HEAT_STOPS[i + 1];
    if (u <= u1) return mixHex(c0, c1, (u - u0) / (u1 - u0));
  }
  return HEAT_STOPS[HEAT_STOPS.length - 1][1];
}
function targetOutputInfo(M, P) {
  if (M.lacExists && M.mlac != null) {
    const x = M.mlac + M.Mt;
    const tau = x > 0 ? M.Mt / ((1 / P.tauEmp) * x) : 0;
    return { out: 24 * M.Mt / (M.T + tau), Mstar: M.Mt, tau, state: "target-met" };
  }
  if (M.emptied && M.Memp != null) {
    const tau = P.tauEmp;
    return { out: 24 * M.Memp / (M.T + tau), Mstar: M.Memp, tau, state: "emptied" };
  }
  return { out: 0, Mstar: 0, tau: null, state: "involution" };
}
function durationOutputInfo(M) {
  if (M.pemp == null) return { out: 0, Mstar: 0, state: "involution" };
  const Mstar = M.Lam * M.pemp;
  return { out: 24 * Mstar / (M.T + M.tau), Mstar, state: "positive" };
}

/* Numerical argmax trace of M_DAY in schedule space. For each T-column of
   an output grid, choose the smallest schedule parameter attaining the
   numerical maximum. This selects the lower edge of a flat optimum plateau. */
function outputMaximumTrace(grid) {
  if (!grid) return [];
  const pts = [];
  for (let i = 0; i < grid.NX; i++) {
    let vmax = -Infinity;
    for (let j = 0; j < grid.NY; j++) {
      const v = grid.values[j * grid.NX + i];
      if (isFinite(v)) vmax = Math.max(vmax, v);
    }
    if (!(vmax > 1e-10)) { pts.push([null, null]); continue; }

    const tol = 1e-7 * Math.max(1, Math.abs(vmax));
    let jbest = 0;
    while (jbest < grid.NY - 1 && grid.values[jbest * grid.NX + i] < vmax - tol) jbest++;

    let jf = jbest;
    // Quadratic sub-grid refinement for an isolated interior maximum.
    if (jbest > 0 && jbest < grid.NY - 1) {
      const vm = grid.values[(jbest - 1) * grid.NX + i];
      const v0 = grid.values[jbest * grid.NX + i];
      const vp = grid.values[(jbest + 1) * grid.NX + i];
      const den = vm - 2 * v0 + vp;
      if (vm < v0 - tol && vp < v0 - tol && Math.abs(den) > 1e-14) {
        const q = Math.max(-0.5, Math.min(0.5, 0.5 * (vm - vp) / den));
        jf = jbest + q;
      }
    }

    const x = grid.xmin + (i / (grid.NX - 1)) * (grid.xmax - grid.xmin);
    const y = grid.ymin + (jf / (grid.NY - 1)) * (grid.ymax - grid.ymin);
    pts.push([x, y]);
  }
  return pts;
}

function drawOutputMaximumTrace(ctx, pts, X, Y) {
  if (!pts || !pts.length) return;
  ctx.setLineDash([]);
  ctx.beginPath();
  let pen = false;
  for (const q of pts) {
    if (q[0] == null || q[1] == null || !isFinite(q[0]) || !isFinite(q[1])) { pen = false; continue; }
    const px = X(q[0]), py = Y(q[1]);
    if (!pen) { ctx.moveTo(px, py); pen = true; } else ctx.lineTo(px, py);
  }
  // Pale underlay keeps the optimum trace visible on dark parts of the heatmap.
  ctx.strokeStyle = C.panel;
  ctx.lineWidth = 5.0;
  ctx.stroke();
  ctx.strokeStyle = C.opt;
  ctx.lineWidth = 2.4;
  ctx.stroke();
}
function drawHeatField(ctx, w, h, grid, vmax, contourCount, X, Y, colorScale = "linear") {
  const off = document.createElement("canvas");
  off.width = grid.NX; off.height = grid.NY;
  const oc = off.getContext("2d");
  const img = oc.createImageData(grid.NX, grid.NY);
  for (let j = 0; j < grid.NY; j++) {
    for (let i = 0; i < grid.NX; i++) {
      const v = grid.values[j * grid.NX + i];
      const rgb = hexRgb(heatColor(v, vmax, colorScale));
      const py = grid.NY - 1 - j;
      const o = (py * grid.NX + i) * 4;
      img.data[o] = rgb[0]; img.data[o + 1] = rgb[1]; img.data[o + 2] = rgb[2]; img.data[o + 3] = 255;
    }
  }
  oc.putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(off, PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);

  if (contourCount > 0 && vmax > 0) {
    ctx.save();
    ctx.strokeStyle = "rgba(19,33,30,0.72)";
    ctx.lineWidth = 0.9;
    for (let k = 1; k <= contourCount; k++) {
      const level = heatValueAtNorm(k / (contourCount + 1), vmax, colorScale);
      ctx.beginPath();
      for (let j = 0; j < grid.NY - 1; j++) {
        for (let i = 0; i < grid.NX - 1; i++) {
          const v00 = grid.values[j * grid.NX + i];
          const v10 = grid.values[j * grid.NX + i + 1];
          const v11 = grid.values[(j + 1) * grid.NX + i + 1];
          const v01 = grid.values[(j + 1) * grid.NX + i];
          const pts = [];
          const hit = (va, vb, xa, ya, xb, yb) => {
            if ((va < level && vb >= level) || (vb < level && va >= level)) {
              const q = (level - va) / (vb - va);
              pts.push([xa + q * (xb - xa), ya + q * (yb - ya)]);
            }
          };
          hit(v00, v10, i, j, i + 1, j);
          hit(v10, v11, i + 1, j, i + 1, j + 1);
          hit(v11, v01, i + 1, j + 1, i, j + 1);
          hit(v01, v00, i, j + 1, i, j);
          const seg = (a, b) => {
            const x1 = grid.xmin + (a[0] / (grid.NX - 1)) * (grid.xmax - grid.xmin);
            const y1 = grid.ymin + (a[1] / (grid.NY - 1)) * (grid.ymax - grid.ymin);
            const x2 = grid.xmin + (b[0] / (grid.NX - 1)) * (grid.xmax - grid.xmin);
            const y2 = grid.ymin + (b[1] / (grid.NY - 1)) * (grid.ymax - grid.ymin);
            ctx.moveTo(X(x1), Y(y1)); ctx.lineTo(X(x2), Y(y2));
          };
          if (pts.length === 2) seg(pts[0], pts[1]);
          else if (pts.length === 4) { seg(pts[0], pts[1]); seg(pts[2], pts[3]); }
        }
      }
      ctx.stroke();
    }
    ctx.restore();
  }

  // compact legend for the heatmap
  const lx = w - PAD.r - 110, ly = PAD.t + 8, lw = 94, lh = 8;
  for (let i = 0; i < lw; i++) {
    const lv = heatValueAtNorm(i / (lw - 1), vmax, colorScale);
    ctx.fillStyle = heatColor(Math.max(lv, 1e-12), vmax, colorScale);
    ctx.fillRect(lx + i, ly, 1.2, lh);
  }
  ctx.strokeStyle = C.rule; ctx.strokeRect(lx, ly, lw, lh);
  ctx.fillStyle = C.ink2; ctx.font = `9px ${MONO}`; ctx.textAlign = "left";
  ctx.fillText("0", lx, ly + 19);
  ctx.textAlign = "right"; ctx.fillText(fmt(vmax), lx + lw, ly + 19);
  ctx.textAlign = "center"; ctx.fillText(`Daily output (ml/day)${colorScale === "log" ? " · log" : ""}`, lx + lw / 2, ly + 32);
}
function hexRgb(c) {
  if (c[0] === "#") return hex(c);
  const m = c.match(/rgb\((\d+),(\d+),(\d+)\)/);
  return m ? [parseInt(m[1]), parseInt(m[2]), parseInt(m[3])] : [0, 0, 0];
}
/* Four orbit panels: values at each iteration (n=0 is the prescribed initial state).
   M_n and τ_n belong to removal session n, so they start at n=1. */
function drawOrbitPanels(ctx, w, h, traj) {
  ctx.fillStyle = C.panel; ctx.fillRect(0, 0, w, h);
  if (!traj) {
    ctx.fillStyle = C.ink3; ctx.textAlign = "center";
    ctx.font = `12px ${SANS}`;
    ctx.fillText("Select an initial condition in the phase plane.", w / 2, h / 2);
    return;
  }
  const nmax = Math.max(1, traj.length - 1);
  const left = 53, right = w - 13, top = 9, bottom = h - 29;
  const rowHeight = (bottom - top) / 4;
  const X = n => left + (n / nmax) * (right - left);
  const rows = [
    { index: 0, label: "m_n (ml)", color: C.ink },
    { index: 1, label: "p_n (ml/h)", color: C.traj },
    { index: 3, label: "M_n (ml)", color: "#B66B15" },
    { index: 2, label: "τ_n (h)", color: C.sep }
  ];
  rows.forEach((r, j) => {
    const yTop = top + j * rowHeight + 9;
    const yBottom = top + (j + 1) * rowHeight - 15;
    const available = traj.map(z => z[r.index]).filter(v => v != null && Number.isFinite(v));
    const vmax = Math.max(1e-6, ...available);
    const Y = v => yBottom - Math.max(0, v) / vmax * (yBottom - yTop);
    ctx.strokeStyle = C.grid; ctx.lineWidth = 1;
    ctx.beginPath();ctx.moveTo(left, yBottom);ctx.lineTo(right, yBottom);ctx.stroke();
    for (let k = 1; k <= 3; k++) {
      const xx = left + (k / 4) * (right-left);
      ctx.beginPath(); ctx.moveTo(xx,yTop);ctx.lineTo(xx,yBottom);ctx.stroke();
    }
    ctx.fillStyle = r.color; ctx.textAlign = "left";ctx.font = `10px ${SANS}`;
    ctx.fillText(r.label, left + 3, yTop - 3);
    ctx.fillStyle = C.ink3;ctx.textAlign = "right";ctx.font = `9px ${MONO}`;
    ctx.fillText(fmt(vmax), left - 5, yTop + 5);
    ctx.fillText("0", left - 5, yBottom + 3);
    ctx.strokeStyle = r.color;ctx.lineWidth = 1.6;ctx.beginPath();
    let started=false;
    traj.forEach((v,i) => {
      const value=v[r.index];
      if (value==null || !Number.isFinite(value)) {started=false;return;}
      if (!started) {ctx.moveTo(X(i),Y(value));started=true;}
      else ctx.lineTo(X(i),Y(value));
    });
    ctx.stroke();
    ctx.fillStyle=r.color;
    traj.forEach((v,i)=>{
      const value=v[r.index];if(value==null||!Number.isFinite(value))return;
      ctx.beginPath();ctx.arc(X(i),Y(value),2.35,0,2*Math.PI);ctx.fill();
    });
  });
  ctx.fillStyle=C.ink2;ctx.textAlign="center";ctx.font=`10px ${MONO}`;
  for(let j=0;j<=4;j++){
    const k=Math.round(nmax*j/4);ctx.fillText(String(k),X(k),h-17);
  }
  ctx.font=`11px ${SANS}`;ctx.fillText("iteration n",(left+right)/2,h-3);
}

function FixedPoints({points}) {
  return <div className="w-[245px] shrink-0 rounded-sm p-3" style={{background:C.panel,border:`1px solid ${C.rule}`}}>
    <div style={{fontFamily:SANS,fontSize:12,fontWeight:600,color:C.ink,marginBottom:9}}>Fixed points</div>
    <table style={{width:"100%",fontFamily:MONO,fontSize:11,color:C.ink,borderCollapse:"collapse",fontVariantNumeric:"tabular-nums"}}>
      <thead><tr>{["", "m (ml)", "p (ml/h)", "M (ml)"].map(t => <th key={t} style={{fontWeight:500,color:C.ink2,textAlign:"right",paddingBottom:7,fontSize:9.5}}>{t}</th>)}</tr></thead>
      <tbody>{points.map((q,i)=><tr key={i}>
        <td style={{padding:"5px 1px",borderTop:`1px solid ${C.grid}`,fontWeight:600}}>{q.label}</td>
        {[q.m,q.p,q.M].map((v,k)=><td key={k} style={{padding:"5px 0",textAlign:"right",borderTop:`1px solid ${C.grid}`}}>{fmt3(v)}</td>)}
      </tr>)}</tbody>
    </table>
  </div>;
}

/* ------------------------------------------------------------------- panel */
function Panel({ title, hint, children }) {
  return (
    <div
      className="flex-1 min-w-[300px] rounded-sm overflow-hidden"
      style={{ background: C.panel, border: `1px solid ${C.rule}` }}
    >
      <div className="flex items-center justify-between px-3 py-2 gap-2" style={{ borderBottom: `1px solid ${C.rule}` }}>
        <span style={{ fontFamily: SANS, fontSize: 12.5, fontWeight: 600, color: C.ink }}>{title}</span>
        <div className="flex items-center gap-2">
          <span style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink3 }}>{hint}</span>
        </div>
      </div>
      {children}
    </div>
  );
}

function Row({ k, v, strong }) {
  return (
    <div className="flex items-baseline justify-between gap-2 py-[3px]">
      <span style={{ fontFamily: MONO, fontSize: 11, color: C.ink3 }}>{k}</span>
      <span style={{ fontFamily: MONO, fontSize: 11.5, color: strong ? C.sep : C.ink, fontWeight: strong ? 600 : 400 }}>{v}</span>
    </div>
  );
}

function Slider({ label, sym, value, set, min, max, step }) {
  return (
    <label className="block mb-2.5">
      <div className="flex items-baseline justify-between">
        <span style={{ fontFamily: SANS, fontSize: 11, color: C.ink2 }}>
          {label} <span style={{ fontFamily: MONO, color: C.ink3 }}>{sym}</span>
        </span>
        <span style={{ fontFamily: MONO, fontSize: 11, color: C.ink }}>{label === "emptying time" && typeof value === "number" ? value.toFixed(2) : value}</span>
      </div>
      <input
        type="range" min={min} max={max} step={step} value={value}
        onChange={(e) => set(parseFloat(e.target.value))}
        className="w-full mt-1 cursor-pointer"
        style={{ accentColor: C.traj, height: 16 }}
      />
    </label>
  );
}

/* -------------------------------------------------------------------- main */
/* Collapsible bifurcation area; the control remains visible when the plot is hidden. */
function CollapsibleBifurcation({ title, children }) {
  const [open, setOpen] = useState(true);
  return (
    <div className="w-full" style={{width:"100%",minWidth:0}}>
      <div className="flex items-center justify-end mb-1">
        <button type="button" onClick={() => setOpen((v) => !v)}
          aria-expanded={open} className="px-3 py-1 rounded-sm"
          style={{fontFamily:SANS,fontSize:11,color:C.ink2,background:C.panel,border:`1px solid ${C.rule}`}}>
          {open ? "Hide" : "Show"} {title}
        </button>
      </div>
      {/* Keep the canvas in the DOM while collapsed: its ResizeObserver must
          remain attached so reopening restores the full plot width. */}
      <div style={{display:open?"block":"none",width:"100%"}}>{children}</div>
    </div>
  );
}

/* Analytical bifurcation diagrams for amount-dependent stimulation.
   The theorems establish fixed-point convergence: no non-constant attractors occur here.
   Curves show stable (solid) and unstable (dashed) fixed points. */
function AmountBifurcation({ P, kind, T, Tmax, param, paramMax, setT, setParam }) {
  const [vary, setVary] = useState("T");
  const [pAxisMax, setPAxisMax] = useState(250);
  const isTarget = kind === "target";
  const isTauDuration = kind === "tau-duration";
  const pname = isTarget ? "M_tar" : "τ";
  const punit = isTarget ? "ml" : "h";
  const pmin = isTarget ? 0.01 : 0.005;

  // The fixed bifurcation slice IS the selected schedule: no independent slider state.
  const dat = useMemo(() => {
    const N = 650;
    const xmin = vary === "T" ? 0.15 : pmin;
    const xmax = vary === "T" ? Math.max(0.2, Tmax) : Math.max(pmin * 2, paramMax);
    const rows = [];
    for (let i = 0; i < N; i++) {
      const x = xmin + i * (xmax - xmin) / (N - 1);
      const t = vary === "T" ? x : Math.max(0.15, Math.min(T, Tmax));
      const q = vary === "T" ? Math.max(pmin, Math.min(param, paramMax)) : x;
      const mm = isTauDuration ? modelDurTauStim(P,t,q) : isTarget ? model(P,t,q) : modelDur(P,t,q);
      const item = {
        x, zero: isTauDuration ? null : {p:0, stable:mm.involStable},
        target:null, emptied:null, saddle:null,
      };
      if (isTauDuration) {
        // The duration-dependent set-duration model is an affine contraction:
        // its unique positive fixed point is globally attracting for τ>0.
        // Color only records whether the removal session empties the breast.
        item[q >= P.tauEmp ? "emptied" : "target"] = {p:mm.pstar,stable:true};
      } else if (isTarget) {
        // A branch must be a FIXED POINT of the selected removal protocol, not
        // just a root of R(M) = Gamma. These are the admissibility conditions.
        if (mm.lacExists) item.target = {p:mm.plac,stable:true};
        if (mm.emptied) item.emptied = {p:mm.pemp,stable:true};
        if (mm.sig && mm.psep != null && q >= mm.Msep) item.saddle = {p:mm.psep,stable:false};
      } else {
        if (mm.pemp != null) item.target = {p:mm.pemp,stable:true};
        if (mm.sig && mm.psep != null) item.saddle = {p:mm.psep,stable:false};
      }
      rows.push(item);
    }
    return {rows, xmin, xmax};
  },[P,kind,isTauDuration,vary,T,param,Tmax,paramMax]);

  const draw = useCallback((ctx,w,h) => {
    const xl = vary === "T" ? "wait interval  T  (h)" : `${pname}  (${punit})`;
    const {X,Y} = frame(ctx,w,h,xl,"fixed-point production rate  p  (ml/h)",dat.xmax,pAxisMax,dat.xmin,0);
    ctx.save();ctx.beginPath();ctx.rect(PAD.l,PAD.t,w-PAD.l-PAD.r,h-PAD.t-PAD.b);ctx.clip();
    const colors={zero:C.ink2, target:"#309C50", emptied:"#D28B2C", saddle:C.sep};
    for (const key of ["zero","target","emptied","saddle"]) {
      const col=colors[key];
      for (const stable of [true,false]) {
        ctx.strokeStyle=col;ctx.lineWidth=stable?2:1.7;
        ctx.setLineDash(stable?[]:[5,4]);ctx.beginPath();let pen=false;
        for (const row of dat.rows) {
          const z=row[key];
          if (!z||z.stable!==stable||!isFinite(z.p)) {pen=false;continue;}
          const xx=X(row.x), yy=Y(z.p);
          if (!pen) {ctx.moveTo(xx,yy);pen=true;}else ctx.lineTo(xx,yy);
        }
        ctx.stroke();ctx.setLineDash([]);
      }
    }
    if (vary === "T") {
      const crit=criticalWait(P);
      if (crit != null && crit>=dat.xmin && crit<=dat.xmax) {
        const xx=X(crit);
        ctx.strokeStyle=C.ink3;ctx.lineWidth=1;ctx.setLineDash([3,4]);
        ctx.beginPath();ctx.moveTo(xx,PAD.t);ctx.lineTo(xx,h-PAD.b);ctx.stroke();ctx.setLineDash([]);
        ctx.font=`10px ${MONO}`;ctx.textAlign="left";ctx.fillStyle=C.ink2;
        ctx.fillText("T_crit",Math.min(xx+5,w-PAD.r-43),h-PAD.b-7);
      }
    }
    const sel=vary==="T"?T:param;
    if(sel>=dat.xmin&&sel<=dat.xmax){
      ctx.strokeStyle=C.mark;ctx.lineWidth=1;ctx.setLineDash([3,4]);
      ctx.beginPath();ctx.moveTo(X(sel),PAD.t);ctx.lineTo(X(sel),h-PAD.b);ctx.stroke();ctx.setLineDash([]);
    }
    ctx.restore();
  },[P,vary,T,param,pname,punit,dat,pAxisMax]);

  const canvas=useCanvas(285,draw,(px,py,w,h)=>{
    const frac=Math.max(0,Math.min(1,(px-PAD.l)/(w-PAD.l-PAD.r)));
    const val=dat.xmin+frac*(dat.xmax-dat.xmin);
    if(vary==="T") setT(+val.toFixed(3));
    else setParam(+val.toFixed(3));
  });
  const fixedValue = vary === "T" ? param : T;
  return (
    <CollapsibleBifurcation title="bifurcation diagram">
      <Panel title="Analytical bifurcation diagram" hint="fixed-point branches; click to select schedule">
        <div className="flex flex-wrap gap-3 items-center px-3 pt-2 pb-1" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}>
          <div className="flex gap-1">
            {[["T","vary wait interval T"],["parameter",`vary ${pname}`]].map(([key,label])=>(
              <button key={key} onClick={()=>setVary(key)} className="px-2 py-1 rounded-sm"
                style={{background:vary===key?C.ink:"transparent",color:vary===key?C.panel:C.ink2,
                  border:`1px solid ${vary===key?C.ink:C.rule}`}}>{label}</button>
            ))}
          </div>
          <div className="flex-1 min-w-[160px]">
            <Slider label={vary==="T"?`fixed ${pname}`:"fixed T"} sym={vary==="T"?punit:"h"}
              value={fixedValue} set={vary==="T"?setParam:setT}
              min={vary==="T"?pmin:0.15} max={vary==="T"?paramMax:Tmax}
              step={vary==="T"?(isTarget?1:0.005):0.05}/>
          </div>
          <div className="flex-1 min-w-[150px]">
            <Slider label="p-axis maximum" sym="ml/h" value={pAxisMax}
              set={setPAxisMax} min={0.5} max={5000} step={0.5}/>
          </div>
        </div>
        <div ref={canvas.wrap} style={{width:"100%"}}><canvas ref={canvas.cvs} onClick={canvas.click} tabIndex={0}
          style={{cursor:"crosshair",display:"block"}} /></div>
        <div className="flex flex-wrap gap-4 px-3 py-2" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,borderTop:`1px solid ${C.rule}`}}>
          <span><span style={{color:"#309C50"}}>━━</span> {isTauDuration ? "p* (partial emptying)" : isTarget ? "p_tar (stable target-met)" : "p_lac (stable)"}</span>
          {(isTarget || isTauDuration) && <span><span style={{color:"#D28B2C"}}>━━</span> {isTauDuration ? "p* (complete emptying)" : "p_emp (stable emptying)"}</span>}
          {!isTauDuration && P.resp==="sig" && <span><span style={{color:C.sep}}>┄┄</span> p_sep (unstable)</span>}
          {!isTauDuration && <span><span style={{color:C.ink2}}>━━</span> zero-production state</span>}
          <span>{isTauDuration
            ? "Unique attracting fixed point; branch color indicates partial versus complete emptying."
            : isTarget ? "Only admissible branches are shown; p_emp requires M_tar ≥ M_emp(T)."
            : "Only admissible positive fixed-point branches are shown."}</span>
        </div>
      </Panel>
    </CollapsibleBifurcation>
  );
}

function TargetView({ P, setP, T, setT }) {
  const [Mt, setMt] = useState(170);
  const [Tmax, setTmax] = useState(12);
  const [Mmax, setMmax] = useState(600);
  const [zoom, setZoom] = useState(1);
  const [fixedAxes, setFixedAxes] = useState(false);
  const [mapMode, setMapMode] = useState("regions");
  const [outMax, setOutMax] = useState(1200);
  const [outContours, setOutContours] = useState(6);
  const [heatScale, setHeatScale] = useState("linear");
  const [ic, setIc] = useState(null);
  const [orbitShow, setOrbitShow] = useState(50);
  const set = (k) => (v) => setP((q) => ({ ...q, [k]: v }));

  const M = useMemo(() => model(P, T, Mt), [P, T, Mt]);
  const seps = useMemo(() => separatrix(M), [M]);

  /* ---- regime map curves */
  const curves = useMemo(() => {
    const sep = [], emp = [];
    const brk = (arr) => {
      if (arr.length && arr[arr.length - 1][0] !== null) arr.push([null, null]);
    };
    const n = 260;
    for (let i = 0; i <= n; i++) {
      const t = 0.15 + ((Tmax - 0.15) * i) / n;
      const q = model(P, t, Mt);
      if (q.kappa > 1) {
        if (q.sig) sep.push([t, q.Msep]); else brk(sep);
        emp.push([t, q.Memp]);
      } else { brk(sep); brk(emp); }
    }
    return { sep, emp };
  }, [P, Tmax, Mt]);

  const gridIn = useDebounced(useMemo(() => ({ P, Tmax, Mmax }), [P, Tmax, Mmax]));
  const outGrid = useMemo(() => {
    const { P, Tmax, Mmax } = gridIn;
    const NX = 150, NY = 145;
    const values = new Float64Array(NX * NY);
    let observedMax = 0;
    for (let j = 0; j < NY; j++) {
      const mt = (j / (NY - 1)) * Mmax;
      for (let i = 0; i < NX; i++) {
        const t = Math.max(0.05, (i / (NX - 1)) * Tmax);
        const info = targetOutputInfo(model(P, t, Math.max(mt, 1e-8)), P);
        values[j * NX + i] = info.out;
        observedMax = Math.max(observedMax, info.out);
      }
    }
    return { NX, NY, values, xmin: 0, xmax: Tmax, ymin: 0, ymax: Mmax, observedMax };
  }, [gridIn]);
  const optTrace = useMemo(() => outputMaximumTrace(outGrid), [outGrid]);

  /* ---- phase plane ranges.  In "fixed" mode the frame is built from the whole
         M_tar axis, so it does not rescale as M_tar is swept. */
  const view = useMemo(() => {
    const collect = (mm) => {
      const cand = [mm.Mt / mm.d];
      if (mm.mlac != null) cand.push(mm.mlac);
      for (const b of separatrix(mm)) for (const [m] of b) cand.push(m);
      return {
        m: Math.max(...cand.filter(isFinite), 0),
        p: Math.max(mm.pemp || 0, mm.plac || 0, mm.psep || 0, 0),
      };
    };
    let mBest = 0, pBest = 0;
    if (fixedAxes) {
      for (let i = 1; i <= 40; i++) {
        const c = collect(model(P, T, (Mmax * i) / 40));
        mBest = Math.max(mBest, c.m);
        pBest = Math.max(pBest, c.p);
      }
    } else {
      const c = collect(M);
      mBest = c.m;
      pBest = c.p;
    }
    const mMax = Math.max(mBest, 1) * 1.25 * zoom;
    const pMax = Math.max(pBest, (P.amax / (1 - M.g)) * 0.35) * 1.25 * zoom;
    return { mMax, pMax };
  }, [M, P, T, Mmax, zoom, fixedAxes]);

  /* ---- basin raster */
  const basinIn = useDebounced(useMemo(() => ({ M, view }), [M, view]));
  const basins = useMemo(() => {
    const { M, view } = basinIn;
    const NX = 150, NY = 138;
    const img = new Uint8ClampedArray(NX * NY * 4);
    const uniform = !M.sig || M.kappa <= 1;
    const solid = !M.sig && M.kappa > 1 ? 1 : 0;
    const col = [hex(C.invol), hex(M.lacExists ? C.lact : C.empt)];
    for (let j = 0; j < NY; j++) {
      const p = ((NY - 1 - j) / (NY - 1)) * view.pMax;
      for (let i = 0; i < NX; i++) {
        const m = (i / (NX - 1)) * view.mMax;
        const k = uniform ? solid : classify(m, p, M);
        const o = (j * NX + i) * 4;
        img[o] = col[k][0]; img[o + 1] = col[k][1]; img[o + 2] = col[k][2]; img[o + 3] = 255;
      }
    }
    return { NX, NY, img };
  }, [basinIn]);

  const traj = useMemo(() => {
    if (!ic) return null;
    let [m, p] = ic;
    const pts = [[m, p, null, null, null]];
    for (let i = 0; i < orbitShow; i++) {
      const x = M.d*m + M.G*p;
      const amount = Math.min(M.Mt, x);
      const duration = x > 0 ? P.tauEmp * amount/x : P.tauEmp;
      [m, p] = M.step(m, p);
      if (!Number.isFinite(m) || !Number.isFinite(p)) break;
      pts.push([m, p, duration, amount, x > M.Mt]);
    }
    return pts;
  }, [ic, M, P.tauEmp, orbitShow]);

  /* ---- draw: schedule space */
  const drawMap = useCallback(
    (ctx, w, h) => {
      const { X, Y } = frame(ctx, w, h, "wait interval  T  (h)", "target  M_tar  (ml)", Tmax, Mmax);
      ctx.save();
      ctx.beginPath();
      ctx.rect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      ctx.clip();

      if (mapMode === "output") {
        drawHeatField(ctx, w, h, outGrid, outMax, outContours, X, Y, heatScale);
      } else {
        ctx.fillStyle = C.invol;
        ctx.fillRect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
        const n = 240;
        for (let i = 0; i < n; i++) {
          const t0 = 0.15 + ((Tmax - 0.15) * i) / n;
          const t1 = 0.15 + ((Tmax - 0.15) * (i + 1)) / n;
          const md = model(P, (t0 + t1) / 2, Mt);
          if (md.kappa <= 1) continue;
          const x0 = X(t0), x1 = X(t1);
          if (md.sig) {
            ctx.fillStyle = C.lact;
            ctx.fillRect(x0, Y(md.Memp), x1 - x0 + 0.6, Y(md.Msep) - Y(md.Memp));
            ctx.fillStyle = C.empt;
            ctx.fillRect(x0, PAD.t, x1 - x0 + 0.6, Y(md.Memp) - PAD.t);
          } else {
            ctx.fillStyle = C.lact;
            ctx.fillRect(x0, Y(md.Memp), x1 - x0 + 0.6, h - PAD.b - Y(md.Memp));
            ctx.fillStyle = C.empt;
            ctx.fillRect(x0, PAD.t, x1 - x0 + 0.6, Y(md.Memp) - PAD.t);
          }
        }
      }

      const line = (pts, color, dash, width) => {
        ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash);
        ctx.beginPath();
        let pen = false;
        for (const q of pts) {
          if (q[0] === null) { pen = false; continue; }
          const px = X(q[0]), py = Y(q[1]);
          if (!pen) { ctx.moveTo(px, py); pen = true; } else ctx.lineTo(px, py);
        }
        ctx.stroke(); ctx.setLineDash([]);
      };
      line(curves.emp, mapMode === "output" ? C.ink2 : C.ink, [], 1.4);
      if (M.sig) line(curves.sep, mapMode === "output" ? C.ink2 : C.ink, [], 1.4);
      if (mapMode === "output") drawOutputMaximumTrace(ctx, optTrace, X, Y);
      markCriticalWait(ctx, w, h, P, Tmax, X);

      const cx = X(T), cy = Y(Mt);
      ctx.strokeStyle = C.mark; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(PAD.l, cy); ctx.lineTo(w - PAD.r, cy);
      ctx.moveTo(cx, PAD.t); ctx.lineTo(cx, h - PAD.b); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = C.mark; ctx.strokeStyle = C.ink; ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.arc(cx, cy, 4.5, 0, 7); ctx.fill(); ctx.stroke();
      ctx.restore();

      if (mapMode === "regions") {
        ctx.font = `10px ${MONO}`; ctx.textAlign = "left";
        const lab = (pts, txt, color) => {
          const q = pts.filter((z) => z[0] !== null && z[1] < Mmax * 0.97 && z[1] > Mmax * 0.03);
          if (!q.length) return;
          const z = q[Math.floor(q.length * 0.62)];
          ctx.fillStyle = color;
          ctx.fillText(txt, X(z[0]) + 4, Y(z[1]) - 4);
        };
        lab(curves.emp, "M_emp", C.ink);
        if (M.sig) lab(curves.sep, "M_sep", C.ink);
      }
    },
    [P, T, Mt, Tmax, Mmax, curves, M.sig, mapMode, outGrid, optTrace, outMax, outContours, heatScale]
  );

  /* ---- draw: phase plane */
  const drawPhase = useCallback(
    (ctx, w, h) => {
      const { X, Y } = frame(ctx, w, h, "milk after removal  m  (ml)", "production rate  p  (ml/h)", view.mMax, view.pMax);
      ctx.save();
      ctx.beginPath();
      ctx.rect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      ctx.clip();
      // basins
      const off = document.createElement("canvas");
      off.width = basins.NX; off.height = basins.NY;
      off.getContext("2d").putImageData(new ImageData(basins.img, basins.NX, basins.NY), 0, 0);
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(off, PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      // region E boundary:  m d + p G = Mt
      ctx.strokeStyle = C.ink2; ctx.lineWidth = 1.2; ctx.setLineDash([5, 3]);
      ctx.beginPath(); ctx.moveTo(X(M.Mt / M.d), Y(0)); ctx.lineTo(X(0), Y(M.Mt / M.G)); ctx.stroke();
      ctx.setLineDash([]);
      // separatrix branches
      seps.forEach((b, k) => {
        ctx.strokeStyle = C.sep; ctx.lineWidth = k === 0 ? 2.2 : 1.6;
        ctx.setLineDash(k === 0 ? [] : [6, 3]);
        ctx.beginPath();
        b.forEach(([m, p], i) => (i ? ctx.lineTo(X(m), Y(p)) : ctx.moveTo(X(m), Y(p))));
        ctx.stroke(); ctx.setLineDash([]);
      });
      // trajectory
      if (traj) {
        ctx.strokeStyle = C.traj; ctx.lineWidth = 1.2; ctx.globalAlpha = 0.65;
        ctx.beginPath();
        traj.forEach(([m, p], i) => (i ? ctx.lineTo(X(m), Y(p)) : ctx.moveTo(X(m), Y(p))));
        ctx.stroke(); ctx.globalAlpha = 1;
        traj.slice(0, 60).forEach(([m, p], i) => {
          ctx.fillStyle = C.traj; ctx.globalAlpha = Math.max(0.15, 1 - i / 28);
          ctx.beginPath(); ctx.arc(X(m), Y(p), 2.2, 0, 7); ctx.fill();
        });
        ctx.globalAlpha = 1;
        ctx.strokeStyle = C.mark; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(X(traj[0][0]), Y(traj[0][1]), 5, 0, 7); ctx.stroke();
      }
      // steady states
      const pt = (m, p, stable, label) => {
        const x = X(m), y = Y(p);
        ctx.beginPath(); ctx.arc(x, y, 5, 0, 7);
        ctx.fillStyle = stable ? C.ink : C.panel;
        ctx.strokeStyle = C.ink; ctx.lineWidth = 1.6;
        ctx.fill(); ctx.stroke();
        ctx.fillStyle = C.ink; ctx.font = `10px ${MONO}`; ctx.textAlign = "left";
        ctx.fillText(label, x + 8, y - 6);
      };
      pt(0, 0, M.involStable, "z₀");
      if (M.sig && M.psep != null && M.Mt >= M.Msep) pt(0, M.psep, false, "z₁");
      if (M.emptied) pt(0, M.pemp, true, "z₂");
      if (M.lacExists) pt(M.mlac, M.plac, true, "z₂");
      ctx.restore();
    },
    [M, seps, traj, view, basins]
  );

  /* ---- draw: time series */
  const drawSeries = useCallback((ctx,w,h) => drawOrbitPanels(ctx,w,h,traj), [traj]);

  const mapC = useCanvas(300, drawMap, (px, py, w, h) => {
    const t = ((px - PAD.l) / (w - PAD.l - PAD.r)) * Tmax;
    const mv = ((h - PAD.b - py) / (h - PAD.t - PAD.b)) * Mmax;
    setT(Math.min(Tmax, Math.max(0.15, +t.toFixed(2))));
    setMt(Math.max(1, +mv.toFixed(1)));
    setIc(null);
  });
  const phaseC = useCanvas(300, drawPhase, (px, py, w, h) => {
    const m = ((px - PAD.l) / (w - PAD.l - PAD.r)) * view.mMax;
    const p = ((h - PAD.b - py) / (h - PAD.t - PAD.b)) * view.pMax;
    if (m >= 0 && p >= 0) setIc([m, p]);
  });
  const seriesC = useCanvas(380, drawSeries, null);


  const other = M.lacExists ? "→ lactation state" : M.emptied ? "→ emptied state" : "→ involution";
  const fate = !traj
    ? null
    : !M.bistable
    ? M.involStable && !M.lacExists && !M.emptied
      ? "→ involution"
      : other
    : classify(traj[0][0], traj[0][1], M) === 0
    ? "→ involution"
    : other;

  return (
    <div>
      <header className="mb-3">
        <h1 style={{ fontFamily: SERIF, fontSize: 16, color: C.ink, fontWeight: 600, letterSpacing: "-0.01em" }}>
          Amount-dependent stimulation, set-target removal
        </h1>
        <p style={{ fontFamily: SANS, fontSize: 12, color: C.ink2, maxWidth: 720, marginTop: 2 }}>
          <span style={{ fontFamily: MONO }}>m′ = max(x − M_tar, 0), p′ = gp + a(min(M_tar, x)), x = dm + Gp.</span>
        </p>
      </header>

      <div className="flex flex-wrap gap-3 items-start">
        {/* sidebar */}
        <div className="w-[236px] shrink-0 rounded-sm p-3" style={{ background: C.panel, border: `1px solid ${C.rule}` }}>
          <div className="flex gap-1 mb-3">
            {[["con", "Michaelis–Menten"], ["sig", "Hill"]].map(([k, l]) => (
              <button key={k} onClick={() => set("resp")(k)}
                className="flex-1 py-1 rounded-sm transition-colors"
                style={{
                  fontFamily: SANS, fontSize: 11,
                  background: P.resp === k ? C.ink : "transparent",
                  color: P.resp === k ? C.panel : C.ink2,
                  border: `1px solid ${P.resp === k ? C.ink : C.rule}`,
                }}>
                {l}
              </button>
            ))}
          </div>
          <ResponseFormula P={P} mode="amount" />
          <Slider label="milk decay" sym="δ" value={P.dlt} set={set("dlt")} min={0.005} max={0.4} step={0.005} />
          <Slider label="production decay" sym="γ" value={P.gam} set={set("gam")} min={0.005} max={0.4} step={0.005} />
          <Slider label="max stimulation" sym="a_max" value={P.amax} set={set("amax")} min={2} max={80} step={0.5} />
          <Slider label="half-response" sym={<span>μ<sub>M</sub></span>} value={P.muM} set={set("muM")} min={10} max={600} step={5} />
          <Slider label="emptying time" sym={<span>τ<sub>emp</sub></span>} value={P.tauEmp} set={set("tauEmp")} min={0.08} max={5} step={0.01} />
          <div className="my-2" style={{ borderTop: `1px solid ${C.rule}` }} />
          <Slider label="T" sym="h" value={T} set={(v) => { setT(v); setIc(null); }} min={0.15} max={Tmax} step={0.05} />
          <Slider label="M_tar" sym="ml" value={Mt} set={(v) => { setMt(v); setIc(null); }} min={1} max={Mmax} step={1} />
          <div className="my-2" style={{ borderTop: `1px solid ${C.rule}` }} />
          <Slider label="T axis" sym="max" value={Tmax} set={setTmax} min={2} max={30} step={0.5} />
          <Slider label="M axis" sym="max" value={Mmax} set={setMmax} min={50} max={2000} step={10} />
          <Slider label="orbit iterations" sym="N" value={orbitShow} set={setOrbitShow} min={10} max={300} step={10} />
          <div className="my-2" style={{ borderTop: `1px solid ${C.rule}` }} />
          <div className="flex gap-1 mb-2">
            {[['regions','Regions View'],['output','Daily Output View']].map(([k,l]) => (
              <button key={k} onClick={() => setMapMode(k)} className="flex-1 py-1 rounded-sm"
                style={{ fontFamily:SANS, fontSize:11, background:mapMode===k?C.ink:"transparent", color:mapMode===k?C.panel:C.ink2, border:`1px solid ${mapMode===k?C.ink:C.rule}` }}>{l}</button>
            ))}
          </div>
          {mapMode === "output" && <>
            <Slider label="heatmap max" sym="ml/day" value={outMax} set={setOutMax} min={100} max={5000} step={50} />
            <Slider label="contour lines" sym="n" value={outContours} set={setOutContours} min={0} max={12} step={1} />
            <div className="flex gap-1 mb-2">{[["linear","linear color"],["log","log color"]].map(([k,l]) => (
              <button key={k} onClick={() => setHeatScale(k)} className="flex-1 py-1 rounded-sm"
                style={{fontFamily:SANS,fontSize:10.5,background:heatScale===k?C.traj:"transparent",color:heatScale===k?C.panel:C.ink2,border:`1px solid ${heatScale===k?C.traj:C.rule}`}}>{l}</button>
            ))}</div>
            <button onClick={() => setOutMax(Math.max(100, Math.ceil(outGrid.observedMax / 50) * 50))} className="w-full py-1 rounded-sm mb-2"
              style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,border:`1px solid ${C.rule}`}}>auto max ({fmt(outGrid.observedMax)} ml/day)</button>
            <div style={{fontFamily:SANS,fontSize:9.5,color:C.ink3,lineHeight:1.35,marginBottom:8}}>daily milk output at the positive steady state; target-session time uses linear emptying.</div>
          </>}
          <button onClick={() => setFixedAxes((v) => !v)}
            className="w-full py-1 rounded-sm transition-colors"
            style={{
              fontFamily: SANS, fontSize: 11,
              background: fixedAxes ? C.ink : "transparent",
              color: fixedAxes ? C.panel : C.ink2,
              border: `1px solid ${fixedAxes ? C.ink : C.rule}`,
            }}>
            {fixedAxes ? "phase axes: fixed" : "phase axes: auto"}
          </button>
        </div>

        {/* plots */}
        <div className="flex-1 min-w-[300px] flex flex-col gap-3">
          <div className="flex flex-wrap gap-3">
            <Panel title="Schedule space (T, M_tar)" hint="click to select">
              <div ref={mapC.wrap}>
                <canvas ref={mapC.cvs} onClick={mapC.click} tabIndex={0} style={{ cursor: "crosshair", display: "block" }} />
              </div>
              <div className="flex flex-wrap gap-3 px-3 py-2" style={{ borderTop: `1px solid ${C.rule}` }}>
                {mapMode === "regions" ? <>
                  {[[C.invol, "involution only"],
                    [C.lact, M.sig ? "bistable: lactation (target met) / involution" : "lactation (target met)"],
                    [C.empt, M.sig ? "bistable: lactation (breast emptied) / involution" : "lactation (breast emptied)"]].map(([c, l]) => (
                    <span key={l} className="flex items-center gap-1.5" style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink2 }}>
                      <i style={{ width: 11, height: 11, background: c, border: `1px solid ${C.rule}`, display: "inline-block" }} />
                      {l}
                    </span>
                  ))}
                </> : <><span style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}>Heatmap: daily milk output at the positive steady state; {heatScale} color scale and correspondingly spaced contour levels.</span><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:16,height:0,borderTop:`2.4px solid ${C.opt}`,display:"inline-block"}}/>numerical M_DAY argmax (lowest M_tar on a plateau; plotted range)</span></>}
              </div>
            </Panel>

            <Panel title="Phase plane (m, p)" hint="click to drop an initial condition">
              <div ref={phaseC.wrap}>
                <canvas ref={phaseC.cvs} onClick={phaseC.click} tabIndex={0} style={{ cursor: "crosshair", display: "block" }} />
              </div>
              <div className="px-3 pt-2"><Slider label="phase zoom" sym="×" value={zoom} set={setZoom} min={0.4} max={3} step={0.1} /></div>
              <div className="flex flex-wrap gap-3 px-3 py-2" style={{ borderTop: `1px solid ${C.rule}` }}>
                <span className="flex items-center gap-1.5" style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink2 }}>
                  <i style={{ width: 14, height: 0, borderTop: `2px solid ${C.sep}`, display: "inline-block" }} />
                  separatrix, inner arc
                </span>
                <span className="flex items-center gap-1.5" style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink2 }}>
                  <i style={{ width: 14, height: 0, borderTop: `1.5px dashed ${C.sep}`, display: "inline-block" }} />
                  outer branches
                </span>
                <span className="flex items-center gap-1.5" style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink2 }}>
                  <i style={{ width: 14, height: 0, borderTop: `1.5px dashed ${C.ink2}`, display: "inline-block" }} />
                  emptying boundary
                </span>
              </div>
            </Panel>
          </div>

          <div className="flex flex-wrap gap-3">
            <AmountBifurcation P={P} kind="target" T={T} Tmax={Tmax}
              param={Mt} paramMax={Mmax} setT={setT} setParam={setMt} />
          </div>

          <div className="flex flex-wrap gap-3">
            <Panel title="Orbit: mₙ, pₙ, Mₙ, τₙ" hint={fate || "no initial condition"}>
              <div ref={seriesC.wrap}>
                <canvas ref={seriesC.cvs} style={{ display: "block" }} />
              </div>
            </Panel>

            <FixedPoints points={[
              {label:"z₀",m:0,p:0,M:0},
              ...(M.sig && M.psep != null && M.Mt >= M.Msep ? [{label:"z₁",m:0,p:M.psep,M:M.Msep}] : []),
              ...(M.lacExists ? [{label:"z₂",m:M.mlac,p:M.plac,M:M.Mt}]
                : M.emptied ? [{label:"z₂",m:0,p:M.pemp,M:M.Memp}] : [])
            ]} />
          </div>
        </div>
      </div>
    </div>
  );
}


/* ================================================================ DURATION */
/* Fixed removal duration:  F(m,p) = ( x·h , p·g + a(x(1-h)) ),  x = m·d + p·G.
   Smooth, single-branch, and invertible; no region E.  The scalar problem is
   the fixed-target one with G replaced by the effective gain Λ.            */

function modelDur(P, T, tau) {
  const g = Math.exp(-P.gam * T);
  const d = Math.exp(-P.dlt * T);
  const G =
    Math.abs(P.dlt - P.gam) < 1e-9
      ? T * Math.exp(-P.gam * T)
      : (Math.exp(-P.gam * T) - Math.exp(-P.dlt * T)) / (P.dlt - P.gam);
  const h = Math.max(1 - (1 / P.tauEmp) * tau, 0);
  const Lam = (G * (1 - h)) / (1 - h * d);
  const Gamma = (1 - g) / G;
  const kappa = amountRmax(P) / Gamma;              // fixed-target efficiency ratio
  const kapT = kappa * (1 - h) / (1 - h * d);   // κ_τ = R_max / Γ_τ
  const c = P.muM;
  const sig = P.resp === "sig";
  const a = sig
    ? (M) => (P.amax * M * M) / (M * M + c * c)
    : (M) => (P.amax * M) / (M + P.muM);
  const ap = sig
    ? (M) => (2 * P.amax * c * c * M) / Math.pow(M * M + c * c, 2)
    : (M) => (P.amax * P.muM) / Math.pow(M + P.muM, 2);

  let psep = null, pemp = null;
  if (kapT > 1) {
    if (sig) {
      const D = Math.sqrt(Math.pow(P.amax * Lam, 2) - Math.pow(2 * P.muM * (1 - g), 2));
      psep = (P.amax * Lam - D) / (2 * Lam * (1 - g));
      pemp = (P.amax * Lam + D) / (2 * Lam * (1 - g));
    } else {
      pemp = P.amax / (1 - g) - P.muM / Lam;
    }
  }
  const mOf = (p) => (p * G * h) / (1 - d * h);
  const step = (m, p) => {
    const x = m * d + p * G;
    return [x * h, p * g + a(x * (1 - h))];
  };
  const inv = (m, p) => {
    if (h <= 1e-12) return [NaN, NaN];
    const x = m / h;
    const pp = (p - a(x * (1 - h))) / g;
    return [(x - pp * G) / d, pp];
  };
  const jac = (p) => {
    const M = Lam * p;
    return [[d * h, G * h], [ap(M) * d * (1 - h), g + ap(M) * G * (1 - h)]];
  };
  return { P, T, tau, g, d, G, h, Lam, Gamma, kappa, kapT, a, ap, sig,
    psep, pemp, mOf, step, inv, jac,
    involStable: sig || kapT <= 1 };
}

/* minimum removal duration τ_min(T) (κ_τ = 1), in closed form */
function tauCritical(P, T) {
  const g = Math.exp(-P.gam * T);
  const d = Math.exp(-P.dlt * T);
  const G =
    Math.abs(P.dlt - P.gam) < 1e-9
      ? T * Math.exp(-P.gam * T)
      : (Math.exp(-P.gam * T) - Math.exp(-P.dlt * T)) / (P.dlt - P.gam);
  const Gamma = (1 - g) / G;
  const kappa = amountRmax(P) / Gamma;
  if (!(kappa > 1)) return null;
  const hc = (kappa - 1) / (kappa - d);
  if (hc >= 1) return 0;
  return (1 - hc) / (1 / P.tauEmp);
}

/* stable manifold of the interior saddle, by backward iteration */
function separatrixDur(M) {
  if (!M.sig || M.psep == null) return [];
  if (M.h <= 1e-10) return separatrixDurEmptied(M);
  const z1 = [M.mOf(M.psep), M.psep];
  const J = M.jac(M.psep);
  const tr = J[0][0] + J[1][1];
  const dt = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  const D = Math.sqrt(Math.max(tr * tr - 4 * dt, 0));
  const l2 = (tr - D) / 2;
  let v = [J[0][1], l2 - J[0][0]];
  const nv = Math.hypot(v[0], v[1]) || 1;
  v = [v[0] / nv, v[1] / nv];
  const eps = 1e-5 * (Math.hypot(z1[0], z1[1]) + 1);
  let pts = [];
  for (const s of [1, -1]) {
    let cur = [];
    for (let i = 0; i <= 60; i++) {
      const t = eps * Math.pow(1 / Math.max(l2, 1e-9), i / 60);
      cur.push([z1[0] + s * t * v[0], z1[1] + s * t * v[1]]);
    }
    for (let k = 0; k < 40; k++) {
      pts = pts.concat(cur.filter((q) => q[0] >= 0 && q[1] >= 0));
      const nx = [];
      for (const [m, p] of cur) {
        const q = M.inv(m, p);
        if (q[0] >= 0 && q[1] >= 0 && q[0] < 1e7 && q[1] < 1e7) nx.push(q);
      }
      if (nx.length < 2) break;
      cur = nx;
    }
  }
  pts.push(z1);
  pts.sort((A, B) => A[0] - B[0]);
  return pts;
}

/* For tau >= tau_emp (h = 0) the map is F(m,p) = (0, pg + a(x)), x = dm + Gp, so
   every point lands on the p-axis after one session and
   Sigma = { (m,p) : pg + a(dm + Gp) = p_sep }.  Parametrised by x in
   [x0, x1] with x0 = G p_sep (the saddle z1 = (0, p_sep)) and
   x1 = a^{-1}(p_sep) (where Sigma meets the m-axis):
   p = (p_sep - a(x)) / g,   m = (x - G p) / d.                              */
function separatrixDurEmptied(M) {
  const { g, d, G, a, psep, P } = M;
  const x0 = psep * G;
  // Hill inverse; if p_sep >= a_max, Sigma never reaches the m-axis: draw a finite piece.
  const x1 = psep < P.amax ? P.muM * Math.sqrt(psep / (P.amax - psep)) : 8 * x0;
  if (!(x1 > x0)) return [];
  const pts = [];
  const n = 240;
  for (let i = 0; i <= n; i++) {
    const x = x0 + ((x1 - x0) * i) / n;
    const p = Math.max((psep - a(x)) / g, 0);
    const m = Math.max((x - G * p) / d, 0);
    pts.push([m, p]);
  }
  return pts;
}

function classifyDur(m, p, M) {
  for (let i = 0; i < 500; i++) {
    [m, p] = M.step(m, p);
    if (!isFinite(m) || m > 1e14) return 1;
  }
  return p < M.psep ? 0 : 1;
}

function DurationView({ P, setP, T, setT }) {
  const [tau, setTau] = useState(0.2);
  const [Tmax, setTmax] = useState(12);
  const [TauMax, setTauMax] = useState(1.5);
  const [zoom, setZoom] = useState(1);
  const [mapMode, setMapMode] = useState("regions");
  const [outMax, setOutMax] = useState(1200);
  const [outContours, setOutContours] = useState(6);
  const [heatScale, setHeatScale] = useState("linear");
  const [ic, setIc] = useState(null);
  const [orbitShow, setOrbitShow] = useState(50);
  const set = (k) => (v) => setP((q) => ({ ...q, [k]: v }));

  const M = useMemo(() => modelDur(P, T, tau), [P, T, tau]);
  const sep = useMemo(() => separatrixDur(M), [M]);
  const gridIn = useDebounced(useMemo(() => ({ P, Tmax, TauMax }), [P, Tmax, TauMax]));
  const outGrid = useMemo(() => {
    const { P, Tmax, TauMax } = gridIn;
    const NX = 150, NY = 145;
    const values = new Float64Array(NX * NY);
    let observedMax = 0;
    for (let j = 0; j < NY; j++) {
      const tv = (j / (NY - 1)) * TauMax;
      for (let i = 0; i < NX; i++) {
        const t = Math.max(0.05, (i / (NX - 1)) * Tmax);
        const info = durationOutputInfo(modelDur(P, t, tv));
        values[j * NX + i] = info.out;
        observedMax = Math.max(observedMax, info.out);
      }
    }
    return { NX, NY, values, xmin: 0, xmax: Tmax, ymin: 0, ymax: TauMax, observedMax };
  }, [gridIn]);
  const optTrace = useMemo(() => outputMaximumTrace(outGrid), [outGrid]);

  const states = useMemo(() => {
    const out = [{ z: [0, 0], stable: M.involStable, lab: "z₀" }];
    if (M.psep != null) out.push({ z: [M.mOf(M.psep), M.psep], stable: false, lab: "z₁" });
    if (M.pemp != null) out.push({ z: [M.mOf(M.pemp), M.pemp], stable: true, lab: "z₂" });
    return out;
  }, [M]);

  const view = useMemo(() => {
    const cand = [1, ...states.map((s) => s.z[0]), ...sep.map((q) => q[0])];
    const pc = [1, ...states.map((s) => s.z[1]), ...sep.map((q) => q[1])];
    return {
      mMax: Math.max(...cand.filter(isFinite)) * 1.25 * zoom,
      pMax: Math.max(...pc.filter(isFinite)) * 1.25 * zoom,
    };
  }, [states, sep, zoom]);

  const basinIn = useDebounced(useMemo(() => ({ M, view }), [M, view]));
  const basins = useMemo(() => {
    const { M, view } = basinIn;
    const NX = 140, NY = 130;
    const img = new Uint8ClampedArray(NX * NY * 4);
    const uniform = !M.sig || M.psep == null;
    const solid = !M.sig && M.kapT > 1 ? 1 : 0;
    const col = [hex(C.invol), hex(C.lact)];
    for (let j = 0; j < NY; j++) {
      const p = ((NY - 1 - j) / (NY - 1)) * view.pMax;
      for (let i = 0; i < NX; i++) {
        const m = (i / (NX - 1)) * view.mMax;
        const k = uniform ? solid : classifyDur(m, p, M);
        const o = (j * NX + i) * 4;
        img[o] = col[k][0]; img[o + 1] = col[k][1]; img[o + 2] = col[k][2]; img[o + 3] = 255;
      }
    }
    return { NX, NY, img };
  }, [basinIn]);

  const traj = useMemo(() => {
    if (!ic) return null;
    let [m, p] = ic;
    const pts = [[m, p, null, null, null]];
    for (let i = 0; i < orbitShow; i++) {
      const x = M.d*m + M.G*p;
      const amount = x*(1-M.h);
      [m,p] = M.step(m,p);
      if (!Number.isFinite(m) || !Number.isFinite(p)) break;
      pts.push([m,p,M.tau,amount,false]);
    }
    return pts;
  }, [ic, M, orbitShow]);

  /* ---- draw: schedule space (T, tau) */
  const drawMap = useCallback(
    (ctx, w, h) => {
      const { X, Y } = frame(ctx, w, h, "wait interval  T  (h)", "removal duration  τ  (h)", Tmax, TauMax);
      ctx.save();
      ctx.beginPath();
      ctx.rect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      ctx.clip();

      if (mapMode === "output") {
        drawHeatField(ctx, w, h, outGrid, outMax, outContours, X, Y, heatScale);
      } else {
        ctx.fillStyle = C.invol;
        ctx.fillRect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      }

      const n = 300;
      const curve = [];
      for (let i = 0; i <= n; i++) {
        const t = 0.05 + ((Tmax - 0.05) * i) / n;
        const tc = tauCritical(P, t);
        if (tc === null) { curve.push([null, null]); continue; }
        curve.push([t, tc]);
        if (mapMode === "regions") {
          ctx.fillStyle = C.lact;
          const x0 = X(t), x1 = X(t + (Tmax - 0.05) / n) + 0.6;
          ctx.fillRect(x0, PAD.t, x1 - x0, Y(Math.max(tc, 0)) - PAD.t);
          // Complete emptying (τ >= τ_emp) is orange in Fig. 7 of the paper.
          if (P.tauEmp < TauMax) {
            ctx.fillStyle = C.empt;
            ctx.fillRect(x0, PAD.t, x1 - x0, Y(P.tauEmp) - PAD.t);
          }
        }
      }
      ctx.strokeStyle = mapMode === "output" ? C.ink2 : C.ink; ctx.lineWidth = 1.4;
      ctx.beginPath();
      let pen = false;
      for (const q of curve) {
        if (q[0] === null) { pen = false; continue; }
        const px = X(q[0]), py = Y(q[1]);
        if (!pen) { ctx.moveTo(px, py); pen = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();

      {
        const tsat = P.tauEmp;
        if (tsat <= TauMax) {
          ctx.strokeStyle = C.ink3; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
          ctx.beginPath(); ctx.moveTo(PAD.l, Y(tsat)); ctx.lineTo(w - PAD.r, Y(tsat)); ctx.stroke();
          ctx.setLineDash([]);
        }
      }
      if (mapMode === "output") drawOutputMaximumTrace(ctx, optTrace, X, Y);
      markCriticalWait(ctx, w, h, P, Tmax, X);
      const cx = X(T), cy = Y(tau);
      ctx.strokeStyle = C.mark; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(PAD.l, cy); ctx.lineTo(w - PAD.r, cy);
      ctx.moveTo(cx, PAD.t); ctx.lineTo(cx, h - PAD.b); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = C.mark; ctx.strokeStyle = C.ink; ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.arc(cx, cy, 4.5, 0, 7); ctx.fill(); ctx.stroke();
      ctx.restore();

      ctx.fillStyle = C.ink; ctx.font = `10px ${MONO}`; ctx.textAlign = "left";
      const lab = curve.filter((q) => q[0] !== null && q[1] < TauMax * 0.9);
      if (lab.length) {
        const z = lab[Math.floor(lab.length * 0.55)];
        ctx.fillText("τ_min(T)", X(z[0]) + 4, Y(z[1]) - 5);
      }
      if (P.tauEmp <= TauMax) {
        ctx.fillStyle = C.ink3;
        ctx.fillText("τ = τ_emp", X(Tmax * 0.74), Y(P.tauEmp) - 5);
      }
    },
    [P, T, tau, Tmax, TauMax, mapMode, outGrid, optTrace, outMax, outContours, heatScale]
  );

  /* ---- draw: phase plane */
  const drawPhase = useCallback(
    (ctx, w, h) => {
      const { X, Y } = frame(ctx, w, h, "milk after removal  m  (ml)", "production rate  p  (ml h⁻¹)", view.mMax, view.pMax);
      ctx.save();
      ctx.beginPath();
      ctx.rect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      ctx.clip();
      const off = document.createElement("canvas");
      off.width = basins.NX; off.height = basins.NY;
      off.getContext("2d").putImageData(new ImageData(basins.img, basins.NX, basins.NY), 0, 0);
      ctx.drawImage(off, PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      if (sep.length > 1) {
        ctx.strokeStyle = C.sep; ctx.lineWidth = 2.2;
        ctx.beginPath();
        sep.forEach(([m, p], i) => (i ? ctx.lineTo(X(m), Y(p)) : ctx.moveTo(X(m), Y(p))));
        ctx.stroke();
      }
      if (traj) {
        ctx.strokeStyle = C.traj; ctx.lineWidth = 1.2; ctx.globalAlpha = 0.65;
        ctx.beginPath();
        traj.forEach(([m, p], i) => (i ? ctx.lineTo(X(m), Y(p)) : ctx.moveTo(X(m), Y(p))));
        ctx.stroke(); ctx.globalAlpha = 1;
        traj.slice(0, 60).forEach(([m, p], i) => {
          ctx.fillStyle = C.traj; ctx.globalAlpha = Math.max(0.15, 1 - i / 28);
          ctx.beginPath(); ctx.arc(X(m), Y(p), 2.2, 0, 7); ctx.fill();
        });
        ctx.globalAlpha = 1;
        ctx.strokeStyle = C.mark; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(X(traj[0][0]), Y(traj[0][1]), 5, 0, 7); ctx.stroke();
      }
      for (const s of states) {
        ctx.beginPath(); ctx.arc(X(s.z[0]), Y(s.z[1]), 5, 0, 7);
        ctx.fillStyle = s.stable ? C.ink : C.panel;
        ctx.strokeStyle = C.ink; ctx.lineWidth = 1.6;
        ctx.fill(); ctx.stroke();
        ctx.fillStyle = C.ink; ctx.font = `10px ${MONO}`; ctx.textAlign = "left";
        ctx.fillText(s.lab, X(s.z[0]) + 8, Y(s.z[1]) - 6);
      }
      ctx.restore();
    },
    [M, sep, traj, view, basins, states]
  );

  const drawSeries = useCallback((ctx,w,h) => drawOrbitPanels(ctx,w,h,traj), [traj]);

  const mapC = useCanvas(300, drawMap, (px, py, w, h) => {
    setT(Math.max(0.15, +(((px - PAD.l) / (w - PAD.l - PAD.r)) * Tmax).toFixed(2)));
    setTau(Math.max(0.01, +(((h - PAD.b - py) / (h - PAD.t - PAD.b)) * TauMax).toFixed(3)));
    setIc(null);
  });
  const phaseC = useCanvas(300, drawPhase, (px, py, w, h) => {
    const m = ((px - PAD.l) / (w - PAD.l - PAD.r)) * view.mMax;
    const p = ((h - PAD.b - py) / (h - PAD.t - PAD.b)) * view.pMax;
    if (m >= 0 && p >= 0) setIc([m, p]);
  });
  const seriesC = useCanvas(380, drawSeries, null);

  const fate = !traj ? null
    : !M.sig || M.kapT <= 1
      ? (M.kapT > 1 ? "→ lactation state" : "→ involution")
      : classifyDur(traj[0][0], traj[0][1], M) === 0 ? "→ involution" : "→ lactation state";

  return (
    <div>
      <header className="mb-3">
        <h1 style={{ fontFamily: SERIF, fontSize: 16, color: C.ink, fontWeight: 600, letterSpacing: "-0.01em" }}>
          Amount-dependent stimulation, set-duration removal
        </h1>
        <p style={{ fontFamily: SANS, fontSize: 12, color: C.ink2, maxWidth: 760, marginTop: 2 }}>
          <span style={{ fontFamily: MONO }}>m′ = hx, p′ = gp + a((1−h)x), x = dm + Gp, h = max(1−τ/τ_emp, 0).</span>
        </p>
      </header>

      <div className="flex flex-wrap gap-3 items-start">
        <div className="w-[236px] shrink-0 rounded-sm p-3" style={{ background: C.panel, border: `1px solid ${C.rule}` }}>
          <div className="flex gap-1 mb-3">
            {[["con", "Michaelis–Menten"], ["sig", "Hill"]].map(([k, l]) => (
              <button key={k} onClick={() => set("resp")(k)}
                className="flex-1 py-1 rounded-sm transition-colors"
                style={{
                  fontFamily: SANS, fontSize: 11,
                  background: P.resp === k ? C.ink : "transparent",
                  color: P.resp === k ? C.panel : C.ink2,
                  border: `1px solid ${P.resp === k ? C.ink : C.rule}`,
                }}>
                {l}
              </button>
            ))}
          </div>
          <ResponseFormula P={P} mode="amount" />
          <Slider label="milk decay" sym="δ" value={P.dlt} set={set("dlt")} min={0.005} max={0.4} step={0.005} />
          <Slider label="production decay" sym="γ" value={P.gam} set={set("gam")} min={0.005} max={0.4} step={0.005} />
          <Slider label="max stimulation" sym="a_max" value={P.amax} set={set("amax")} min={2} max={80} step={0.5} />
          <Slider label="half-response" sym={<span>μ<sub>M</sub></span>} value={P.muM} set={set("muM")} min={10} max={600} step={5} />
          <Slider label="emptying time" sym={<span>τ<sub>emp</sub></span>} value={P.tauEmp} set={set("tauEmp")} min={0.08} max={5} step={0.01} />
          <div className="my-2" style={{ borderTop: `1px solid ${C.rule}` }} />
          <Slider label="T" sym="h" value={T} set={(v) => { setT(v); setIc(null); }} min={0.15} max={Tmax} step={0.05} />
          <Slider label="τ" sym="h" value={tau} set={(v) => { setTau(v); setIc(null); }} min={0.01} max={TauMax} step={0.01} />
          <div className="my-2" style={{ borderTop: `1px solid ${C.rule}` }} />
          <Slider label="T axis" sym="max" value={Tmax} set={setTmax} min={2} max={30} step={0.5} />
          <Slider label="τ axis" sym="max" value={TauMax} set={setTauMax} min={0.2} max={4} step={0.1} />
          <Slider label="orbit iterations" sym="N" value={orbitShow} set={setOrbitShow} min={10} max={300} step={10} />
          <div className="my-2" style={{ borderTop: `1px solid ${C.rule}` }} />
          <div className="flex gap-1 mb-2">
            {[['regions','Regions View'],['output','Daily Output View']].map(([k,l]) => (
              <button key={k} onClick={() => setMapMode(k)} className="flex-1 py-1 rounded-sm"
                style={{ fontFamily:SANS, fontSize:11, background:mapMode===k?C.ink:"transparent", color:mapMode===k?C.panel:C.ink2, border:`1px solid ${mapMode===k?C.ink:C.rule}` }}>{l}</button>
            ))}
          </div>
          {mapMode === "output" && <>
            <Slider label="heatmap max" sym="ml/day" value={outMax} set={setOutMax} min={100} max={5000} step={50} />
            <Slider label="contour lines" sym="n" value={outContours} set={setOutContours} min={0} max={12} step={1} />
            <div className="flex gap-1 mb-2">{[["linear","linear color"],["log","log color"]].map(([k,l]) => (
              <button key={k} onClick={() => setHeatScale(k)} className="flex-1 py-1 rounded-sm"
                style={{fontFamily:SANS,fontSize:10.5,background:heatScale===k?C.traj:"transparent",color:heatScale===k?C.panel:C.ink2,border:`1px solid ${heatScale===k?C.traj:C.rule}`}}>{l}</button>
            ))}</div>
            <button onClick={() => setOutMax(Math.max(100, Math.ceil(outGrid.observedMax / 50) * 50))} className="w-full py-1 rounded-sm mb-2"
              style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,border:`1px solid ${C.rule}`}}>auto max ({fmt(outGrid.observedMax)} ml/day)</button>
          </>}
        </div>

        <div className="flex-1 min-w-[300px] flex flex-col gap-3">
          <div className="flex flex-wrap gap-3">
            <Panel title="Schedule space (T, τ)" hint="click to select">
              <div ref={mapC.wrap}>
                <canvas ref={mapC.cvs} onClick={mapC.click} tabIndex={0} style={{ cursor: "crosshair", display: "block" }} />
              </div>
              <div className="flex flex-wrap gap-3 px-3 py-2" style={{ borderTop: `1px solid ${C.rule}` }}>
                {mapMode === "regions" ? [[C.invol, "involution only"],
                  [C.lact, M.sig ? "bistable: lactation / involution" : "lactation"],
                  [C.empt, M.sig ? "bistable: lactation (breast emptied) / involution" : "lactation (breast emptied)"]].map(([c, l]) => (
                  <span key={l} className="flex items-center gap-1.5" style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink2 }}>
                    <i style={{ width: 11, height: 11, background: c, border: `1px solid ${C.rule}`, display: "inline-block" }} />
                    {l}
                  </span>
                )) : <><span style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}>Heatmap: daily milk output at the positive steady state; {heatScale} color scale and correspondingly spaced contour levels.</span><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:16,height:0,borderTop:`2.4px solid ${C.opt}`,display:"inline-block"}}/>numerical M_DAY argmax (within plotted τ range)</span></>}
              </div>
            </Panel>

            <Panel title="Phase plane (m, p)" hint="click to drop an initial condition">
              <div ref={phaseC.wrap}>
                <canvas ref={phaseC.cvs} onClick={phaseC.click} tabIndex={0} style={{ cursor: "crosshair", display: "block" }} />
              </div>
              <div className="px-3 pt-2"><Slider label="phase zoom" sym="×" value={zoom} set={setZoom} min={0.4} max={3} step={0.1} /></div>
              <div className="flex flex-wrap gap-3 px-3 py-2" style={{ borderTop: `1px solid ${C.rule}` }}>
                <span className="flex items-center gap-1.5" style={{ fontFamily: SANS, fontSize: 10.5, color: C.ink2 }}>
                  <i style={{ width: 14, height: 0, borderTop: `2px solid ${C.sep}`, display: "inline-block" }} />
                  separatrix Σ = W^s(z₁)
                </span>
              </div>
            </Panel>
          </div>

          <div className="flex flex-wrap gap-3">
            <AmountBifurcation P={P} kind="duration" T={T} Tmax={Tmax}
              param={tau} paramMax={TauMax} setT={setT} setParam={setTau} />
          </div>

          <div className="flex flex-wrap gap-3">
            <Panel title="Orbit: mₙ, pₙ, Mₙ, τₙ" hint={fate || "no initial condition"}>
              <div ref={seriesC.wrap}>
                <canvas ref={seriesC.cvs} style={{ display: "block" }} />
              </div>
            </Panel>

            <FixedPoints points={[
              {label:"z₀",m:0,p:0,M:0},
              ...(M.psep != null ? [{label:"z₁",m:M.mOf(M.psep),p:M.psep,M:M.Lam*M.psep}] : []),
              ...(M.pemp != null ? [{label:"z₂",m:M.mOf(M.pemp),p:M.pemp,M:M.Lam*M.pemp}] : [])
            ]} />
          </div>
        </div>
      </div>
    </div>
  );
}

/* ===================================================== DURATION-DEPENDENT STIMULATION */
function waitPars(P, T) {
  const g = Math.exp(-P.gam * T);
  const d = Math.exp(-P.dlt * T);
  const G = Math.abs(P.dlt - P.gam) < 1e-9
    ? T * Math.exp(-P.gam * T)
    : (Math.exp(-P.gam * T) - Math.exp(-P.dlt * T)) / (P.dlt - P.gam);
  return { g, d, G, Gamma: (1 - g) / G };
}

function tauStimFns(P) {
  const mu = Math.max(P.muTau == null ? 0.2 : P.muTau, 1e-8);
  const sig = P.resp === "sig";
  const c = mu;
  const a = sig
    ? (t) => (P.amax * t * t) / (t * t + c * c)
    : (t) => (P.amax * t) / (t + mu);
  const ap = sig
    ? (t) => (2 * P.amax * c * c * t) / Math.pow(t * t + c * c, 2)
    : (t) => (P.amax * mu) / Math.pow(t + mu, 2);
  return { mu, sig, a, ap };
}

function modelDurTauStim(P, T, tau) {
  const { g, d, G, Gamma } = waitPars(P, T);
  const { a, ap, sig, mu } = tauStimFns(P);
  const h = Math.max(1 - (1 / P.tauEmp) * tau, 0);
  const r = (1 - h) / (1 - d * h);
  const aval = a(tau);
  const pstar = aval / (1 - g);
  const mstar = (pstar * G * h) / (1 - d * h);
  const Mstar = G * r * pstar;
  const out = 24 * Mstar / (T + tau);
  const step = (m, p) => {
    const x = d * m + G * p;
    return [x * h, g * p + aval];
  };
  return { P, T, tau, g, d, G, Gamma, h, r, a, ap, sig, mu,
    aval, pstar, mstar, Mstar, out, step };
}

function modelTargetTauStim(P, T, Mt) {
  const { g, d, G, Gamma } = waitPars(P, T);
  const { a, ap, sig, mu } = tauStimFns(P);
  const tauE = P.tauEmp;
  const aE = a(tauE);
  const pemp = aE / (1 - g);
  const Memp = aE / Gamma;

  const step = (m, p) => {
    const x = d * m + G * p;
    let dur, mn, hitTarget;
    if (x > Mt && x > 1e-14) {
      dur = Math.min(Mt / ((1 / P.tauEmp) * x), tauE);
      mn = x - Mt;
      hitTarget = true;
    } else {
      dur = tauE;
      mn = 0;
      hitTarget = false;
    }
    return [mn, g * p + a(dur), dur, x, hitTarget];
  };

  if (!(Mt > 0) || Mt >= Memp * (1 - 1e-10)) {
    const out = 24 * Memp / (T + tauE);
    return { P, T, Mt, g, d, G, Gamma, a, ap, sig, mu, tauE, aE, Memp,
      type: "emptied", tau: tauE, h: 0, r: 1, mstar: 0, pstar: pemp,
      stable: true, flipMargin: null, out, step };
  }

  const f = (tt) => {
    const h = 1 - (1 / P.tauEmp) * tt;
    const r = (1 - h) / (1 - d * h);
    return (r * a(tt)) / Gamma - Mt;
  };
  let lo = 1e-10, hi = tauE;
  for (let k = 0; k < 60; k++) {
    const mid = (lo + hi) / 2;
    if (f(mid) < 0) lo = mid; else hi = mid;
  }
  const tau = (lo + hi) / 2;
  const h = Math.max(1 - (1 / P.tauEmp) * tau, 0);
  const r = (1 - h) / (1 - d * h);
  const aval = a(tau);
  const pstar = aval / (1 - g);
  const mstar = (pstar * G * h) / (1 - d * h);
  const xstar = mstar + Mt;
  const cfb = (tau * ap(tau)) / Math.max(xstar, 1e-14);
  const tr = d + g - cfb * G;
  const det = d * g;
  const disc = tr * tr - 4 * det;
  let eig1, eig2;
  if (disc >= 0) {
    const D = Math.sqrt(disc);
    eig1 = (tr + D) / 2;
    eig2 = (tr - D) / 2;
  } else {
    const im = Math.sqrt(-disc) / 2;
    eig1 = { re: tr / 2, im };
    eig2 = { re: tr / 2, im: -im };
  }
  const absEig = (z) => typeof z === "number" ? Math.abs(z) : Math.hypot(z.re, z.im);
  const rho = Math.max(absEig(eig1), absEig(eig2));
  const flipMargin = (1 + d) * (1 + g) - cfb * G;
  const stable = rho < 1 - 1e-7;
  const out = 24 * Mt / (T + tau);
  return { P, T, Mt, g, d, G, Gamma, a, ap, sig, mu, tauE, aE, Memp,
    type: "target-met", tau, h, r, aval, mstar, pstar, stable, flipMargin, out, step };
}

function targetTauStimAttractorOutput(mm, burn = 420, sample = 140) {
  // Stable fixed points and the emptied state have exact daily output.
  if (mm.type === "emptied" || mm.stable) {
    return { out: mm.out, state: mm.type === "emptied" ? "emptied" : "stable fixed point", numerical: false };
  }

  // Beyond the flip, estimate long-term output on the attracting orbit reached
  // from a small perturbation of the (unstable) target-met fixed point.
  let m = Math.max(0, (mm.mstar || 0) + 1e-3 * Math.max(1, Math.abs(mm.mstar || 0)));
  let p = Math.max(0, (mm.pstar || 0) - 1e-3 * Math.max(1, Math.abs(mm.pstar || 0)));
  let milk = 0, time = 0, targetOnly = true;
  const total = burn + sample;
  for (let k = 0; k < total; k++) {
    const q = mm.step(m, p);
    if (!isFinite(q[0]) || !isFinite(q[1]) || q[0] > 1e14 || q[1] > 1e14) {
      return { out: NaN, state: "numerical failure", numerical: true };
    }
    if (k >= burn) {
      const removed = Math.max(0, q[3] - q[0]);
      milk += removed;
      time += mm.T + q[2];
      if (!q[4]) targetOnly = false;
    }
    m = q[0]; p = q[1];
  }
  return {
    out: time > 0 ? 24 * milk / time : 0,
    state: targetOnly ? "non-steady target-met" : "non-steady incl. emptying",
    numerical: true
  };
}

function detectTwoCycle(traj) {
  if (!traj || traj.length < 16) return false;
  const q = traj.slice(-10);
  let d2 = 0, d1 = 0, scale = 1;
  for (let i = 2; i < q.length; i++) {
    const a = q[i], b = q[i - 2], c = q[i - 1];
    d2 = Math.max(d2, Math.hypot(a[0] - b[0], a[1] - b[1]));
    d1 = Math.max(d1, Math.hypot(a[0] - c[0], a[1] - c[1]));
    scale = Math.max(scale, Math.hypot(a[0], a[1]));
  }
  return d2 < 2e-5 * scale && d1 > 2e-3 * scale;
}

function classifyTargetTauStimAttractor(mm, total=240, keep=12) {
  if (mm.type === "emptied") return "emptied";
  if (mm.stable) return "stable";
  let m = Math.max(0, (mm.mstar || 0) + 1e-3 * Math.max(1, Math.abs(mm.mstar || 0)));
  let p = Math.max(0, (mm.pstar || 0) - 1e-3 * Math.max(1, Math.abs(mm.pstar || 0)));
  const flags = [];
  for (let k = 0; k < total; k++) {
    const q = mm.step(m, p);
    m = q[0]; p = q[1];
    if (!isFinite(m) || !isFinite(p) || m > 1e14 || p > 1e14) break;
    if (k >= total - keep) flags.push(!!q[4]);
  }
  return flags.some(v => !v) ? "mixed" : "target-only";
}

function flipTargetAtT(P, T) {
  if (P.resp !== "sig") return null;
  const { g, d, G, Gamma } = waitPars(P, T);
  const { a, ap } = tauStimFns(P);
  const tauE = P.tauEmp;
  const margin = (tau) => {
    const h = Math.max(1 - (1 / P.tauEmp) * tau, 0);
    const aval = a(tau);
    if (!(aval > 0)) return 1;
    const E = tau * ap(tau) / aval;
    return (1 + d) * (1 + g) - (1 - g) * (1 - d * h) * E;
  };
  const n = 160;
  let t0 = 1e-6 * tauE, f0 = margin(t0);
  for (let i = 1; i <= n; i++) {
    const t1 = (i / n) * tauE;
    const f1 = margin(t1);
    if (f0 <= 0 && f1 > 0) {
      let lo = t0, hi = t1;
      for (let k = 0; k < 45; k++) {
        const mid = (lo + hi) / 2;
        if (margin(mid) <= 0) lo = mid; else hi = mid;
      }
      const tau = (lo + hi) / 2;
      const h = 1 - (1 / P.tauEmp) * tau;
      const r = (1 - h) / (1 - d * h);
      return { tau, Mt: r * a(tau) / Gamma };
    }
    t0 = t1; f0 = f1;
  }
  return null;
}

function TauStimDurationView({ P, setP, T, setT }) {
  const [tau, setTau] = useState(0.2);
  const [Tmax, setTmax] = useState(12);
  const [TauMax, setTauMax] = useState(2.0);
  const [zoom, setZoom] = useState(1);
  const [outMax, setOutMax] = useState(1200);
  const [outContours, setOutContours] = useState(6);
  const [heatScale, setHeatScale] = useState("linear");
  const [mapMode, setMapMode] = useState("regions");
  const [ic, setIc] = useState(null);
  const [orbitShow, setOrbitShow] = useState(50);
  const set = (k) => (v) => setP((q) => ({ ...q, [k]: v }));
  const M = useMemo(() => modelDurTauStim(P, T, tau), [P, T, tau]);

  const gridIn = useDebounced(useMemo(() => ({ P, Tmax, TauMax }), [P, Tmax, TauMax]));
  const outGrid = useMemo(() => {
    const { P, Tmax, TauMax } = gridIn;
    const NX = 145, NY = 140;
    const values = new Float64Array(NX * NY);
    let observedMax = 0;
    for (let j = 0; j < NY; j++) {
      const tv = Math.max(1e-4, (j / (NY - 1)) * TauMax);
      for (let i = 0; i < NX; i++) {
        const t = Math.max(0.05, (i / (NX - 1)) * Tmax);
        const mm = modelDurTauStim(P, t, tv);
        values[j * NX + i] = mm.out;
        observedMax = Math.max(observedMax, mm.out);
      }
    }
    return { NX, NY, values, xmin: 0, xmax: Tmax, ymin: 0, ymax: TauMax, observedMax };
  }, [gridIn]);
  const optTrace = useMemo(() => outputMaximumTrace(outGrid), [outGrid]);

  const view = useMemo(() => ({
    mMax: Math.max(1, M.mstar * 2.2) * zoom,
    pMax: Math.max(1, M.pstar * 2.2) * zoom,
  }), [M, zoom]);

  const traj = useMemo(() => {
    if (!ic) return null;
    let [m, p] = ic;
    const pts = [[m, p, null, null, null]];
    for (let i = 0; i < orbitShow; i++) {
      const x = M.d*m + M.G*p;
      const amount = x*(1-M.h);
      [m,p] = M.step(m,p);
      if (!Number.isFinite(m) || !Number.isFinite(p)) break;
      pts.push([m,p,M.tau,amount,false]);
    }
    return pts;
  }, [ic, M, orbitShow]);

  const drawMap = useCallback((ctx, w, h) => {
    const { X, Y } = frame(ctx, w, h, "wait interval  T  (h)", "removal duration  τ  (h)", Tmax, TauMax);
    ctx.save();
    ctx.beginPath(); ctx.rect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b); ctx.clip();
    if (mapMode === "output") {
      drawHeatField(ctx, w, h, outGrid, outMax, outContours, X, Y, heatScale);
    } else {
      // Duration-dependent stimulation always has a positive, attracting equilibrium.
      // Only the removal protocol changes at complete emptying.
      ctx.fillStyle = C.lact;
      ctx.fillRect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
      if (P.tauEmp < TauMax) {
        ctx.fillStyle = C.empt;
        ctx.fillRect(PAD.l, PAD.t, w - PAD.l - PAD.r, Y(P.tauEmp) - PAD.t);
      }
    }
    if (P.tauEmp <= TauMax) {
      ctx.strokeStyle = C.ink3; ctx.lineWidth = 1.2; ctx.setLineDash([5, 3]);
      ctx.beginPath(); ctx.moveTo(PAD.l, Y(P.tauEmp)); ctx.lineTo(w - PAD.r, Y(P.tauEmp)); ctx.stroke(); ctx.setLineDash([]);
    }
    if (mapMode === "output") drawOutputMaximumTrace(ctx, optTrace, X, Y);
    const cx = X(T), cy = Y(tau);
    ctx.strokeStyle = C.mark; ctx.lineWidth = 1; ctx.setLineDash([2,3]);
    ctx.beginPath(); ctx.moveTo(PAD.l, cy); ctx.lineTo(w - PAD.r, cy); ctx.moveTo(cx, PAD.t); ctx.lineTo(cx, h - PAD.b); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = C.mark; ctx.strokeStyle = C.ink; ctx.beginPath(); ctx.arc(cx, cy, 4.5, 0, 7); ctx.fill(); ctx.stroke();
    ctx.restore();
    if (P.tauEmp <= TauMax) {
      ctx.fillStyle = C.ink3; ctx.font = `10px ${MONO}`; ctx.textAlign = "left";
      ctx.fillText("τ = τ_emp", X(Tmax * 0.72), Y(P.tauEmp) - 5);
    }
  }, [P, T, tau, Tmax, TauMax, mapMode, outGrid, optTrace, outMax, outContours, heatScale]);

  const drawPhase = useCallback((ctx, w, h) => {
    const { X, Y } = frame(ctx, w, h, "milk after removal  m  (ml)", "production rate  p  (ml h⁻¹)", view.mMax, view.pMax);
    ctx.fillStyle = C.lact; ctx.fillRect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b);
    if (traj) {
      ctx.strokeStyle = C.traj; ctx.lineWidth = 1.2; ctx.globalAlpha = .7; ctx.beginPath();
      traj.forEach(([m,p],i) => i ? ctx.lineTo(X(m),Y(p)) : ctx.moveTo(X(m),Y(p))); ctx.stroke(); ctx.globalAlpha = 1;
      traj.slice(0,50).forEach(([m,p],i) => { ctx.fillStyle=C.traj;ctx.globalAlpha=Math.max(.15,1-i/25);ctx.beginPath();ctx.arc(X(m),Y(p),2.1,0,7);ctx.fill(); });
      ctx.globalAlpha=1;
    }
    ctx.beginPath(); ctx.arc(X(M.mstar), Y(M.pstar), 5, 0, 7); ctx.fillStyle=C.ink; ctx.strokeStyle=C.ink; ctx.fill(); ctx.stroke();
    ctx.fillStyle=C.ink;ctx.font=`10px ${MONO}`;ctx.textAlign="left";ctx.fillText("z₂",X(M.mstar)+8,Y(M.pstar)-6);
  }, [M, traj, view]);

  const drawSeries=useCallback((ctx,w,h)=>drawOrbitPanels(ctx,w,h,traj),[traj]);

  const mapC=useCanvas(300,drawMap,(px,py,w,h)=>{setT(Math.max(.15,+(((px-PAD.l)/(w-PAD.l-PAD.r))*Tmax).toFixed(2)));setTau(Math.max(.01,+(((h-PAD.b-py)/(h-PAD.t-PAD.b))*TauMax).toFixed(3)));setIc(null);});
  const phaseC=useCanvas(300,drawPhase,(px,py,w,h)=>{const m=((px-PAD.l)/(w-PAD.l-PAD.r))*view.mMax,p=((h-PAD.b-py)/(h-PAD.t-PAD.b))*view.pMax;if(m>=0&&p>=0)setIc([m,p]);});
  const seriesC=useCanvas(380,drawSeries,null);

  return <div>
    <header className="mb-3"><h1 style={{fontFamily:SERIF,fontSize:16,color:C.ink,fontWeight:600}}>Duration-dependent stimulation, set-duration removal</h1>
      <p style={{fontFamily:SANS,fontSize:12,color:C.ink2,maxWidth:800,marginTop:2}}><span style={{fontFamily:MONO}}>m′ = hx, p′ = gp + a(τ), x = dm + Gp, h = max(1−τ/τ_emp, 0).</span></p></header>
    <div className="flex flex-wrap gap-3 items-start">
      <div className="w-[236px] shrink-0 rounded-sm p-3" style={{background:C.panel,border:`1px solid ${C.rule}`}}>
        <div className="flex gap-1 mb-3">{[["con","Michaelis–Menten"],["sig","Hill"]].map(([k,l])=><button key={k} onClick={()=>set("resp")(k)} className="flex-1 py-1 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,background:P.resp===k?C.ink:"transparent",color:P.resp===k?C.panel:C.ink2,border:`1px solid ${P.resp===k?C.ink:C.rule}`}}>{l}</button>)}</div>
        <ResponseFormula P={P} mode="duration" />
        <Slider label="milk decay" sym="δ" value={P.dlt} set={set("dlt")} min={.005} max={.4} step={.005}/>
        <Slider label="production decay" sym="γ" value={P.gam} set={set("gam")} min={.005} max={.4} step={.005}/>
        <Slider label="max stimulation" sym="a_max" value={P.amax} set={set("amax")} min={2} max={80} step={.5}/>
        <Slider label="duration half-response" sym={<span>μ<sub>τ</sub></span>} value={P.muTau} set={set("muTau")} min={.02} max={4} step={.02}/>
        <Slider label="emptying time" sym={<span>τ<sub>emp</sub></span>} value={P.tauEmp} set={set("tauEmp")} min={0.08} max={5} step={0.01} />
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/>
        <Slider label="T" sym="h" value={T} set={v=>{setT(v);setIc(null);}} min={.15} max={Tmax} step={.05}/>
        <Slider label="τ" sym="h" value={tau} set={v=>{setTau(v);setIc(null);}} min={.01} max={TauMax} step={.01}/>
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/>
        <Slider label="T axis" sym="max" value={Tmax} set={setTmax} min={2} max={30} step={.5}/>
        <Slider label="τ axis" sym="max" value={TauMax} set={setTauMax} min={.2} max={5} step={.1}/>
          <Slider label="orbit iterations" sym="N" value={orbitShow} set={setOrbitShow} min={10} max={300} step={10} />
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/>
        <div className="flex gap-1 mb-2">{[["regions","Regions View"],["output","Daily Output View"]].map(([k,l])=><button key={k} onClick={()=>setMapMode(k)} className="flex-1 py-1 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,background:mapMode===k?C.ink:"transparent",color:mapMode===k?C.panel:C.ink2,border:`1px solid ${mapMode===k?C.ink:C.rule}`}}>{l}</button>)}</div>
        {mapMode === "output" && <>
          <Slider label="heatmap max" sym="ml/day" value={outMax} set={setOutMax} min={100} max={5000} step={50}/>
          <Slider label="contour lines" sym="n" value={outContours} set={setOutContours} min={0} max={12} step={1}/>
          <div className="flex gap-1 mb-2">{[["linear","linear color"],["log","log color"]].map(([k,l])=><button key={k} onClick={()=>setHeatScale(k)} className="flex-1 py-1 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,background:heatScale===k?C.traj:"transparent",color:heatScale===k?C.panel:C.ink2,border:`1px solid ${heatScale===k?C.traj:C.rule}`}}>{l}</button>)}</div>
          <button onClick={()=>setOutMax(Math.max(100,Math.ceil(outGrid.observedMax/50)*50))} className="w-full py-1 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,border:`1px solid ${C.rule}`}}>auto max ({fmt(outGrid.observedMax)} ml/day)</button>
        </>}
      </div>
      <div className="flex-1 min-w-[300px] flex flex-col gap-3">
        <div className="flex flex-wrap gap-3"><Panel title="Schedule space (T, τ)" hint={mapMode === "regions" ? "complete / incomplete emptying" : "daily output (ml/day)"}><div ref={mapC.wrap}><canvas ref={mapC.cvs} onClick={mapC.click} tabIndex={0} style={{cursor:"crosshair",display:"block"}}/></div><div className="flex flex-wrap gap-3 px-3 py-2" style={{borderTop:`1px solid ${C.rule}`}}>{mapMode === "regions" ? <><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:11,height:11,background:C.lact,border:`1px solid ${C.rule}`}}/>partial emptying</span><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:11,height:11,background:C.empt,border:`1px solid ${C.rule}`}}/>complete emptying</span></> : <span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:16,height:0,borderTop:`2.4px solid ${C.opt}`,display:"inline-block"}}/>numerical M_DAY argmax (within plotted τ range)</span>}<span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:14,height:0,borderTop:`1.2px dashed ${C.ink3}`,display:"inline-block"}}/>τ = τ_emp</span></div></Panel>
        <Panel title="Phase plane (m, p)" hint="click to drop an initial condition"><div ref={phaseC.wrap}><canvas ref={phaseC.cvs} onClick={phaseC.click} tabIndex={0} style={{cursor:"crosshair",display:"block"}}/></div><div className="px-3 pt-2"><Slider label="phase zoom" sym="×" value={zoom} set={setZoom} min={0.4} max={4} step={0.1}/></div></Panel></div>
        <div className="flex flex-wrap gap-3"><AmountBifurcation P={P} kind="tau-duration" T={T} Tmax={Tmax} param={tau} paramMax={TauMax} setT={setT} setParam={setTau}/></div>
        <div className="flex flex-wrap gap-3"><Panel title="Orbit: mₙ, pₙ, Mₙ, τₙ" hint="→ globally attracting fixed point"><div ref={seriesC.wrap}><canvas ref={seriesC.cvs} style={{display:"block"}}/></div></Panel>
          <FixedPoints points={[{label:"z₂",m:M.mstar,p:M.pstar,M:M.Mstar}]} />
        </div>
      </div>
    </div>
  </div>;
}

function TauStimTargetView({ P, setP, T, setT }) {
  const [Mt, setMt] = useState(40);
  const [Tmax, setTmax] = useState(14);
  const [Mmax, setMmax] = useState(300);
  const [zoom, setZoom] = useState(1);
  const [ic, setIc] = useState(null);
  const [bifMode, setBifMode] = useState("T"); // "T": fix T and vary M_tar; "M": fix M_tar and vary T
  const [bifT, setBifT] = useState(T);
  const [bifM, setBifM] = useState(Mt);
  const [bifIter, setBifIter] = useState(800);
  const [bifKeep, setBifKeep] = useState(30);
  const [bifSamples, setBifSamples] = useState(280);
  const [pBifMax, setPBifMax] = useState(200);
  const [orbitShow, setOrbitShow] = useState(80);
  const [mapMode, setMapMode] = useState("regions");
  const [outMax, setOutMax] = useState(1200);
  const [outContours, setOutContours] = useState(6);
  const [heatScale, setHeatScale] = useState("log");
  const set = (k) => (v) => setP((q) => ({ ...q, [k]: v }));
  const M = useMemo(() => modelTargetTauStim(P, T, Mt), [P, T, Mt]);

  // Keep the bifurcation slice synchronized with the selected schedule.
  useEffect(() => { setBifT(T); }, [T]);
  useEffect(() => { setBifM(Mt); }, [Mt]);

  useEffect(() => {
    if (bifT > Tmax) setBifT(Tmax);
    if (bifT < 0.15) setBifT(0.15);
  }, [Tmax]);
  useEffect(() => {
    if (bifM > Mmax) setBifM(Mmax);
    if (bifM < 0.01) setBifM(0.01);
  }, [Mmax]);

  const gridIn = useDebounced(useMemo(() => ({ P, Tmax, Mmax }), [P, Tmax, Mmax]));
  const regionGrid = useMemo(() => {
    if (mapMode !== "regions") return null;
    const { P, Tmax, Mmax } = gridIn;
    const NX=118,NY=110,img=new Uint8ClampedArray(NX*NY*4);
    const cols={stable:hex(C.lact),targetOnly:hex(C.unstTT),mixed:hex(C.unstET),emptied:hex(C.empt)};
    for(let j=0;j<NY;j++){
      const mt=Math.max(1e-6,(j/(NY-1))*Mmax);
      for(let i=0;i<NX;i++){
        const t=Math.max(.05,(i/(NX-1))*Tmax);const mm=modelTargetTauStim(P,t,mt);
        let c;
        if(mm.type==="emptied") c=cols.emptied;
        else if(mm.stable) c=cols.stable;
        else {
          const cls=classifyTargetTauStimAttractor(mm);
          // Fig. 8: mixed target-met / emptying non-constant orbits are checkered.
          c=cls==="mixed"
            ? ((Math.floor(i/4)+Math.floor(j/4))%2 ? cols.targetOnly : cols.emptied)
            : cols.targetOnly;
        }
        const o=((NY-1-j)*NX+i)*4;
        img[o]=c[0];img[o+1]=c[1];img[o+2]=c[2];img[o+3]=255;
      }
    }
    return {NX,NY,img};
  },[gridIn,mapMode]);

  const outGrid = useMemo(() => {
    if (mapMode !== "output") return null;
    const { P, Tmax, Mmax } = gridIn;
    const NX = 112, NY = 104;
    const values = new Float64Array(NX * NY);
    let observedMax = 0;
    for (let j = 0; j < NY; j++) {
      const mt = Math.max(1e-6, (j / (NY - 1)) * Mmax);
      for (let i = 0; i < NX; i++) {
        const t = Math.max(0.05, (i / (NX - 1)) * Tmax);
        const mm = modelTargetTauStim(P, t, mt);
        const info = targetTauStimAttractorOutput(mm);
        const v = isFinite(info.out) ? Math.max(0, info.out) : 0;
        values[j * NX + i] = v;
        observedMax = Math.max(observedMax, v);
      }
    }
    return { NX, NY, values, xmin: 0, xmax: Tmax, ymin: 0, ymax: Mmax, observedMax };
  }, [gridIn, mapMode]);
  const optTrace = useMemo(() => outputMaximumTrace(outGrid), [outGrid]);

  const curves=useMemo(()=>{
    const emp=[],flip=[];const n=260;
    for(let i=0;i<=n;i++){
      const t=.05+((Tmax-.05)*i)/n;const mm=modelTargetTauStim(P,t,1e12);emp.push([t,mm.Memp]);
      const ff=flipTargetAtT(P,t);flip.push(ff?[t,ff.Mt]:[null,null]);
    }
    return {emp,flip};
  },[P,Tmax]);

  const bifIn=useDebounced(useMemo(()=>({P,bifMode,bifT,bifM,bifIter,bifKeep,bifSamples,Tmax,Mmax}),
    [P,bifMode,bifT,bifM,bifIter,bifKeep,bifSamples,Tmax,Mmax]));
  const bifData=useMemo(()=>{
    const {P,bifMode,bifT,bifM,bifIter,bifKeep,bifSamples,Tmax,Mmax}=bifIn;
    const N=bifSamples;
    const total=Math.max(100,Math.round(bifIter));
    const keep=Math.max(1,Math.min(Math.round(bifKeep),total));
    // Retain the full K iterations numerically but plot at most ~320 per parameter
    // to keep the 720 x 2000 Fig. 8(c) demonstration responsive in a browser.
    const plotStride = Math.max(1, Math.ceil(keep / 320));
    const points=[];
    const fixed=[];
    const flips=[];
    let ymax=0;
    const xmin=0.01; // model requires T>0; Figure 8(c) runs essentially from 0 to 30 h
    const xmax=bifMode==="T"?Math.max(0.02,Mmax):Math.max(0.16,Tmax);
    let prev=null;
    for(let i=0;i<N;i++){
      const par=xmin+((xmax-xmin)*i)/(N-1);
      const t=bifMode==="T"?Math.min(Tmax,Math.max(.15,bifT)):par;
      const mt=bifMode==="T"?par:Math.min(Mmax,Math.max(.01,bifM));
      const mm=modelTargetTauStim(P,t,mt);
      fixed.push({x:par,y:mm.pstar,stable:mm.stable,type:mm.type,margin:mm.flipMargin});
      ymax=Math.max(ymax,mm.pstar||0);

      // Start each parameter value close to its fixed point, but not exactly on it.
      // This makes the diagram display the attractor reached from a small perturbation
      // of the fixed point rather than relying on continuation from the previous sample.
      let m=Math.max(0,(mm.mstar||0)+1e-3*Math.max(1,Math.abs(mm.mstar||0)));
      let pp=Math.max(0,(mm.pstar||0)-1e-3*Math.max(1,Math.abs(mm.pstar||0)));
      const ys=[];
      let ok=true;
      for(let k=0;k<total;k++){
        const q=mm.step(m,pp);m=q[0];pp=q[1];
        if(!isFinite(m)||!isFinite(pp)||m>1e14||pp>1e14){ok=false;break;}
        if(k>=total-keep){
          if ((k-(total-keep))%plotStride===0 || k===total-1) ys.push({y:pp,hitTarget:!!q[4]});
          ymax=Math.max(ymax,pp);
        }
      }
      if(ok) points.push({x:par,ys});

      // Mark analytically predicted lambda=-1 crossings along the selected slice.
      if(mm.type==="target-met" && mm.flipMargin!=null){
        if(prev && prev.type==="target-met" && prev.margin!=null && prev.margin*mm.flipMargin<0){
          const q=Math.abs(prev.margin)/(Math.abs(prev.margin)+Math.abs(mm.flipMargin));
          flips.push(prev.x+q*(par-prev.x));
        }
        prev={x:par,margin:mm.flipMargin,type:mm.type};
      } else {
        prev=null;
      }
    }
    return {points,fixed,flips,xmin,xmax,ymax:Math.max(1,ymax*1.08),total,keep};
  },[bifIn]);

  const traj=useMemo(()=>{
    if(!ic)return null;let [m,p]=ic;const pts=[[m,p,null,null,null]];
    for(let i=0;i<orbitShow;i++){
      const q=M.step(m,p);
      const removed=Math.max(0,q[3]-q[0]);
      const hitTarget=!!q[4];
      m=q[0];p=q[1];pts.push([m,p,q[2],removed,hitTarget]);
      if(!isFinite(m)||!isFinite(p)||m>1e12||p>1e12)break;
    }
    return pts;
  },[ic,M,orbitShow]);
  const has2=useMemo(()=>detectTwoCycle(traj),[traj]);

  const view=useMemo(()=>{
    let mBest=Math.max(1,M.mstar||0),pBest=Math.max(1,M.pstar||0);
    if(traj){for(const q of traj.slice(-120)){if(isFinite(q[0]))mBest=Math.max(mBest,q[0]);if(isFinite(q[1]))pBest=Math.max(pBest,q[1]);}}
    return {mMax:mBest*1.35*zoom,pMax:pBest*1.35*zoom};
  },[M,traj,zoom]);

  const drawMap=useCallback((ctx,w,h)=>{
    const {X,Y}=frame(ctx,w,h,"wait interval  T  (h)","target  M_tar  (ml)",Tmax,Mmax);
    ctx.save();ctx.beginPath();ctx.rect(PAD.l,PAD.t,w-PAD.l-PAD.r,h-PAD.t-PAD.b);ctx.clip();
    if (mapMode === "output" && outGrid) {
      drawHeatField(ctx,w,h,outGrid,outMax,outContours,X,Y,heatScale);
    } else if (regionGrid) {
      const off=document.createElement("canvas");off.width=regionGrid.NX;off.height=regionGrid.NY;off.getContext("2d").putImageData(new ImageData(regionGrid.img,regionGrid.NX,regionGrid.NY),0,0);ctx.drawImage(off,PAD.l,PAD.t,w-PAD.l-PAD.r,h-PAD.t-PAD.b);
    }
    const line=(pts,color,dash,width)=>{ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);ctx.beginPath();let pen=false;for(const q of pts){if(q[0]===null){pen=false;continue;}const xx=X(q[0]),yy=Y(q[1]);if(!pen){ctx.moveTo(xx,yy);pen=true;}else ctx.lineTo(xx,yy);}ctx.stroke();ctx.setLineDash([]);};
    line(curves.emp,mapMode === "output" ? C.ink2 : C.ink,[5,3],1.3);line(curves.flip,C.flip,[],1.6);
    if (mapMode === "output") drawOutputMaximumTrace(ctx, optTrace, X, Y);
    // Selected slice for the numerical bifurcation diagram.
    ctx.strokeStyle=C.traj;ctx.lineWidth=1.7;ctx.setLineDash([7,4]);ctx.beginPath();
    if(bifMode==="T"){
      const sx=X(Math.min(Tmax,Math.max(.15,bifT)));ctx.moveTo(sx,PAD.t);ctx.lineTo(sx,h-PAD.b);
    }else{
      const sy=Y(Math.min(Mmax,Math.max(.01,bifM)));ctx.moveTo(PAD.l,sy);ctx.lineTo(w-PAD.r,sy);
    }
    ctx.stroke();ctx.setLineDash([]);
    const cx=X(T),cy=Y(Mt);ctx.strokeStyle=C.mark;ctx.lineWidth=1;ctx.setLineDash([2,3]);ctx.beginPath();ctx.moveTo(PAD.l,cy);ctx.lineTo(w-PAD.r,cy);ctx.moveTo(cx,PAD.t);ctx.lineTo(cx,h-PAD.b);ctx.stroke();ctx.setLineDash([]);ctx.fillStyle=C.mark;ctx.strokeStyle=C.ink;ctx.beginPath();ctx.arc(cx,cy,4.5,0,7);ctx.fill();ctx.stroke();ctx.restore();
    if (mapMode === "regions") {
      ctx.font=`10px ${MONO}`;ctx.textAlign="left";ctx.fillStyle=C.ink;const e=curves.emp.filter(q=>q[0]!==null&&q[1]<.9*Mmax);if(e.length){const z=e[Math.floor(.55*e.length)];ctx.fillText("emptying boundary",X(z[0])+4,Y(z[1])-5);}const f=curves.flip.filter(q=>q[0]!==null&&q[1]<.9*Mmax);if(f.length){const z=f[Math.floor(.55*f.length)];ctx.fillStyle=C.flip;ctx.fillText("λ = −1",X(z[0])+4,Y(z[1])-5);}
    }
  },[P,T,Mt,Tmax,Mmax,regionGrid,curves,bifMode,bifT,bifM,mapMode,outGrid,optTrace,outMax,outContours,heatScale]);

  const drawPhase=useCallback((ctx,w,h)=>{
    const {X,Y}=frame(ctx,w,h,"milk after removal  m  (ml)","production rate  p  (ml h⁻¹)",view.mMax,view.pMax);
    ctx.fillStyle=M.type==="emptied"?C.empt:(M.stable?C.lact:C.unstTT);ctx.fillRect(PAD.l,PAD.t,w-PAD.l-PAD.r,h-PAD.t-PAD.b);
    // boundary x = M_tar
    ctx.strokeStyle=C.ink2;ctx.lineWidth=1.1;ctx.setLineDash([5,3]);ctx.beginPath();ctx.moveTo(X(Mt/M.d),Y(0));ctx.lineTo(X(0),Y(Mt/M.G));ctx.stroke();ctx.setLineDash([]);
    if(traj){ctx.strokeStyle=C.traj;ctx.lineWidth=1.15;ctx.globalAlpha=.65;ctx.beginPath();traj.forEach(([m,p],i)=>i?ctx.lineTo(X(m),Y(p)):ctx.moveTo(X(m),Y(p)));ctx.stroke();ctx.globalAlpha=1;traj.slice(0,70).forEach(([m,p,,,,hitTarget],i)=>{ctx.fillStyle=i===0?C.ink3:(hitTarget?C.traj:C.mark);ctx.globalAlpha=Math.max(.18,1-i/34);ctx.beginPath();ctx.arc(X(m),Y(p),2.1,0,7);ctx.fill();});ctx.globalAlpha=1;
      if(has2&&traj.length>3){for(const q of traj.slice(-2)){ctx.beginPath();ctx.arc(X(q[0]),Y(q[1]),6,0,7);ctx.fillStyle=C.panel;ctx.strokeStyle=C.sep;ctx.lineWidth=2;ctx.fill();ctx.stroke();}}
    }
    ctx.beginPath();ctx.arc(X(M.mstar),Y(M.pstar),5,0,7);ctx.fillStyle=M.stable?C.ink:C.panel;ctx.strokeStyle=C.ink;ctx.lineWidth=1.6;ctx.fill();ctx.stroke();ctx.fillStyle=C.ink;ctx.font=`10px ${MONO}`;ctx.textAlign="left";ctx.fillText("z₂",X(M.mstar)+8,Y(M.pstar)-6);
  },[M,Mt,traj,has2,view]);

  const drawSeries=useCallback((ctx,w,h)=>drawOrbitPanels(ctx,w,h,traj),[traj]);

  const drawBif=useCallback((ctx,w,h)=>{
    ctx.fillStyle=C.panel;ctx.fillRect(0,0,w,h);
    const xlab=bifMode==="T"?"target  M_tar  (ml)":"wait interval  T  (h)";
    const {X,Y}=frame(ctx,w,h,xlab,"asymptotic production  p_n  (ml h⁻¹)",bifData.xmax,pBifMax,bifData.xmin,0);
    ctx.save();ctx.beginPath();ctx.rect(PAD.l,PAD.t,w-PAD.l-PAD.r,h-PAD.t-PAD.b);ctx.clip();

    // Numerical asymptotic iterates (last K out of N iterations from a small perturbation of z*).
    ctx.globalAlpha=.55;
    for(const q of bifData.points){
      const xx=X(q.x);
      for(const yy of q.ys){ctx.fillStyle=yy.hitTarget?C.traj:C.mark;ctx.fillRect(xx-.85,Y(yy.y)-.85,1.7,1.7);}
    }
    ctx.globalAlpha=1;

    const branch=(stable,color,dash,width)=>{
      ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);ctx.beginPath();let pen=false;
      for(const q of bifData.fixed){
        if(q.stable!==stable){pen=false;continue;}
        const xx=X(q.x),yy=Y(q.y);
        if(!pen){ctx.moveTo(xx,yy);pen=true;}else ctx.lineTo(xx,yy);
      }
      ctx.stroke();ctx.setLineDash([]);
    };
    branch(true,C.ink,[],1.2);
    branch(false,C.unstTT,[5,3],1.4);

    // Analytical lambda=-1 crossings along this slice.
    ctx.strokeStyle=C.flip;ctx.lineWidth=1;ctx.setLineDash([2,3]);
    for(const xf of bifData.flips){const xx=X(xf);ctx.beginPath();ctx.moveTo(xx,PAD.t);ctx.lineTo(xx,h-PAD.b);ctx.stroke();}
    ctx.setLineDash([]);

    // Current selected schedule on the varied axis.
    const selected=bifMode==="T"?Mt:T;
    if(selected>=bifData.xmin&&selected<=bifData.xmax){
      const xx=X(selected);ctx.strokeStyle=C.mark;ctx.lineWidth=1.1;ctx.setLineDash([4,3]);ctx.beginPath();ctx.moveTo(xx,PAD.t);ctx.lineTo(xx,h-PAD.b);ctx.stroke();ctx.setLineDash([]);
    }
    ctx.restore();
  },[bifMode,bifData,Mt,T,pBifMax]);

  const mapC=useCanvas(300,drawMap,(px,py,w,h)=>{setT(Math.max(.15,+(((px-PAD.l)/(w-PAD.l-PAD.r))*Tmax).toFixed(2)));setMt(Math.max(.01,+(((h-PAD.b-py)/(h-PAD.t-PAD.b))*Mmax).toFixed(2)));setIc(null);});
  const phaseC=useCanvas(320,drawPhase,(px,py,w,h)=>{const m=((px-PAD.l)/(w-PAD.l-PAD.r))*view.mMax,p=((h-PAD.b-py)/(h-PAD.t-PAD.b))*view.pMax;if(m>=0&&p>=0)setIc([m,p]);});
  const seriesC=useCanvas(380,drawSeries,null);
  const bifC=useCanvas(330,drawBif,(px,py,w,h)=>{
    const q=(px-PAD.l)/(w-PAD.l-PAD.r);
    const par=bifData.xmin+Math.max(0,Math.min(1,q))*(bifData.xmax-bifData.xmin);
    if(bifMode==="T"){
      setT(Math.min(Tmax,Math.max(.15,bifT)));setMt(Math.max(.01,+par.toFixed(3)));
    }else{
      setMt(Math.min(Mmax,Math.max(.01,bifM)));setT(Math.max(.15,+par.toFixed(3)));
    }
    setIc(null);
  });

  // Figure 8(c): Hill duration stimulation, 720 T-values from 0 to 30 h,
  // 4000 iterations and the last 2000 kept for the bifurcation diagram.
  const loadDemo=()=>{
    const dp={...P,dlt:0.15,gam:0.15,amax:50,muTau:0.5,tauEmp:1/3,resp:"sig"};
    const demoT=25, demoM=0.01;
    setP(dp);setT(demoT);setMt(demoM);
    setTmax(30);setMmax(30);setZoom(1);
    setMapMode("regions");
    setBifMode("M");setBifM(demoM);setBifT(demoT);
    setBifSamples(720);setBifIter(4000);setBifKeep(2000);setPBifMax(2.5);
    setOrbitShow(180);
    const dm=modelTargetTauStim(dp,demoT,demoM);
    setIc([Math.max(0,dm.mstar+0.001),Math.max(0,dm.pstar-0.001)]);
  };
  const fate=!traj?null:(has2?"→ stable 2-cycle":M.stable?"→ fixed point":"fixed point unstable; iterate to inspect");


  return <div>
    <header className="mb-3"><h1 style={{fontFamily:SERIF,fontSize:16,color:C.ink,fontWeight:600}}>Duration-dependent stimulation, set-target removal</h1>
      <p style={{fontFamily:SANS,fontSize:12,color:C.ink2,maxWidth:860,marginTop:2}}><span style={{fontFamily:MONO}}>m′ = max(x−M_tar, 0), p′ = gp + a(τ(x)), x = dm + Gp; τ(x) = τ_emp if x ≤ M_tar, otherwise τ(x) = M_tar·τ_emp/x.</span></p></header>
    <div className="flex flex-wrap gap-3 items-start">
      <div className="w-[236px] shrink-0 rounded-sm p-3" style={{background:C.panel,border:`1px solid ${C.rule}`}}>
        <div className="flex gap-1 mb-3">{[["con","Michaelis–Menten"],["sig","Hill"]].map(([k,l])=><button key={k} onClick={()=>{set("resp")(k);setIc(null);}} className="flex-1 py-1 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,background:P.resp===k?C.ink:"transparent",color:P.resp===k?C.panel:C.ink2,border:`1px solid ${P.resp===k?C.ink:C.rule}`}}>{l}</button>)}</div>
        <ResponseFormula P={P} mode="duration" />
        <Slider label="milk decay" sym="δ" value={P.dlt} set={set("dlt")} min={.005} max={.4} step={.005}/><Slider label="production decay" sym="γ" value={P.gam} set={set("gam")} min={.005} max={.4} step={.005}/><Slider label="max stimulation" sym="a_max" value={P.amax} set={set("amax")} min={2} max={80} step={.5}/><Slider label="duration half-response" sym={<span>μ<sub>τ</sub></span>} value={P.muTau} set={set("muTau")} min={.02} max={4} step={.02}/><Slider label="emptying time" sym={<span>τ<sub>emp</sub></span>} value={P.tauEmp} set={set("tauEmp")} min={0.08} max={5} step={0.01} />
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/><Slider label="T" sym="h" value={T} set={v=>{setT(v);setIc(null);}} min={.15} max={Tmax} step={.05}/><Slider label="M_tar" sym="ml" value={Mt} set={v=>{setMt(v);setIc(null);}} min={.01} max={Mmax} step={.1}/>
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/><Slider label="T axis" sym="max" value={Tmax} set={setTmax} min={2} max={30} step={.5}/><Slider label="M axis" sym="max" value={Mmax} set={setMmax} min={10} max={1000} step={10}/>
          <Slider label="orbit iterations" sym="N" value={orbitShow} set={setOrbitShow} min={10} max={300} step={10} />
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/>
        <div className="flex gap-1 mb-2">{[["regions","Regions View"],["output","Daily Output View"]].map(([k,l])=><button key={k} onClick={()=>setMapMode(k)} className="flex-1 py-1 rounded-sm" style={{fontFamily:SANS,fontSize:11,background:mapMode===k?C.ink:"transparent",color:mapMode===k?C.panel:C.ink2,border:`1px solid ${mapMode===k?C.ink:C.rule}`}}>{l}</button>)}</div>
        {mapMode === "output" && <>
          <Slider label="heatmap max" sym="ml/day" value={outMax} set={setOutMax} min={50} max={5000} step={50}/>
          <Slider label="contour lines" sym="n" value={outContours} set={setOutContours} min={0} max={12} step={1}/>
          <div className="flex gap-1 mb-2">{[["linear","linear color"],["log","log color"]].map(([k,l])=><button key={k} onClick={()=>setHeatScale(k)} className="flex-1 py-1 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,background:heatScale===k?C.traj:"transparent",color:heatScale===k?C.panel:C.ink2,border:`1px solid ${heatScale===k?C.traj:C.rule}`}}>{l}</button>)}</div>
          <button onClick={()=>outGrid&&setOutMax(Math.max(50,Math.ceil(outGrid.observedMax/50)*50))} className="w-full py-1 rounded-sm mb-2" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,border:`1px solid ${C.rule}`}}>auto max ({outGrid?fmt(outGrid.observedMax):"—"} ml/day)</button>
          <div style={{fontFamily:SANS,fontSize:9.5,color:C.ink3,lineHeight:1.35,marginBottom:8}}>Stable and emptied states use the exact fixed-state output. In unstable regions, daily output is the numerical long-time average on the attractor reached from a small perturbation of the fixed point.</div>
        </>}
        <div className="my-2" style={{borderTop:`1px solid ${C.rule}`}}/><button onClick={loadDemo} className="w-full py-1.5 rounded-sm" style={{fontFamily:SANS,fontSize:10.5,color:C.panel,background:C.sep,border:`1px solid ${C.sep}`}}>Load non-constant attractor demo</button>
      </div>
      <div className="flex-1 min-w-[300px] flex flex-col gap-3">
        <div className="flex flex-wrap gap-3"><Panel title="Schedule space (T, M_tar)" hint={mapMode === "output" ? "long-term daily output" : "stable / unstable attractor type / emptied"}><div ref={mapC.wrap}><canvas ref={mapC.cvs} onClick={mapC.click} tabIndex={0} style={{cursor:"crosshair",display:"block"}}/></div><div className="flex flex-wrap gap-3 px-3 py-2" style={{borderTop:`1px solid ${C.rule}`}}>{mapMode === "regions" ? <>{[[C.lact,"stable target-met"],[C.unstTT,"unstable: target-met / target-met"],[C.unstET,"unstable: emptying / target-met"],[C.empt,"emptied"]].map(([c,l])=><span key={l} className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:11,height:11,background:c,border:`1px solid ${C.rule}`,display:"inline-block"}}/>{l}</span>)}<span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:14,height:0,borderTop:`2px dashed ${C.traj}`,display:"inline-block"}}/>bifurcation slice</span></> : <><span style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}>Heatmap: realized long-term daily output; {heatScale} color scale.</span><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:16,height:0,borderTop:`2.4px solid ${C.opt}`,display:"inline-block"}}/>numerical M_DAY argmax (lowest M_tar on a plateau; plotted range)</span><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:14,height:0,borderTop:`1.5px dashed ${C.ink2}`,display:"inline-block"}}/>emptying boundary</span><span className="flex items-center gap-1.5" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}><i style={{width:14,height:0,borderTop:`1.5px solid ${C.flip}`,display:"inline-block"}}/>λ = −1</span></>}</div></Panel>
        <Panel title="Phase plane (m, p)" hint="click to drop an initial condition"><div ref={phaseC.wrap}><canvas ref={phaseC.cvs} onClick={phaseC.click} tabIndex={0} style={{cursor:"crosshair",display:"block"}}/></div><div className="px-3 pt-2"><Slider label="phase zoom" sym="×" value={zoom} set={setZoom} min={0.4} max={4} step={0.1}/></div><div className="flex flex-wrap gap-4 px-3 py-2" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,borderTop:`1px solid ${C.rule}`}}><span>Dashed line: target-met / complete-removal boundary.</span><span className="flex items-center gap-1.5"><i style={{width:8,height:8,background:C.traj,display:"inline-block",borderRadius:99}}/>target-met iterate</span><span className="flex items-center gap-1.5"><i style={{width:8,height:8,background:C.mark,display:"inline-block",borderRadius:99}}/>emptying iterate</span><span>Outlined circles mark a detected attracting 2-cycle.</span></div></Panel></div>
        <div className="flex flex-wrap gap-3"><CollapsibleBifurcation title="bifurcation diagram"><Panel title={`Numerical bifurcation diagram · ${bifMode==="T"?`T = ${fmt3(bifT)} h`:`M_tar = ${fmt3(bifM)} ml`}`} hint="click to select schedule">
          <div className="flex flex-wrap gap-3 items-center px-3 pt-2 pb-1" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}>
            <div className="flex gap-1">
              {[["M","vary wait interval T"],["T","vary M_tar"]].map(([k,label])=>(
                <button key={k} onClick={()=>setBifMode(k)} className="px-2 py-1 rounded-sm"
                  style={{background:bifMode===k?C.ink:"transparent",color:bifMode===k?C.panel:C.ink2,border:`1px solid ${bifMode===k?C.ink:C.rule}`}}>{label}</button>
              ))}
            </div>
            <div className="flex-1 min-w-[160px]">
              <Slider label={bifMode==="T"?"fixed T":"fixed M_tar"} sym={bifMode==="T"?"h":"ml"}
                value={bifMode==="T"?T:Mt} set={bifMode==="T"?v=>{setBifT(v);setT(v);setIc(null);}:v=>{setBifM(v);setMt(v);setIc(null);}}
                min={bifMode==="T"?0.15:0.01} max={bifMode==="T"?Tmax:Mmax} step={bifMode==="T"?0.05:0.1}/>
            </div>
            <div className="flex-1 min-w-[150px]">
              <Slider label="p-axis maximum" sym="ml/h" value={pBifMax} set={setPBifMax} min={0.5} max={5000} step={0.5}/>
            </div>
          </div>
          <div className="flex flex-wrap gap-3 items-center px-3 pb-2" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2}}>
            <div className="flex-1 min-w-[130px]"><Slider label="parameter samples" sym="n" value={bifSamples} set={setBifSamples} min={160} max={800} step={40}/></div>
            <div className="flex-1 min-w-[130px]"><Slider label="iterations" sym="N" value={bifIter} set={setBifIter} min={100} max={5000} step={50}/></div>
            <div className="flex-1 min-w-[130px]"><Slider label="iterates kept" sym="K" value={bifKeep} set={setBifKeep} min={10} max={2000} step={10}/></div>
          </div>
          <div className="px-3 pb-2" style={{fontFamily:SANS,fontSize:9.5,color:C.ink3,lineHeight:1.35}}>
            Each sample starts near its fixed point. The last K iterates are computed; a subset is plotted when K is large.
          </div>
          <div ref={bifC.wrap} style={{width:"100%"}}><canvas ref={bifC.cvs} onClick={bifC.click} tabIndex={0} style={{cursor:"crosshair",display:"block"}}/></div><div className="flex flex-wrap gap-4 px-3 py-2" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,borderTop:`1px solid ${C.rule}`}}><span className="flex items-center gap-1.5"><i style={{width:8,height:8,background:C.traj,display:"inline-block"}}/>target-met iterate</span><span className="flex items-center gap-1.5"><i style={{width:8,height:8,background:C.mark,display:"inline-block"}}/>emptying iterate</span><span className="flex items-center gap-1.5"><i style={{width:14,height:0,borderTop:`1.5px solid ${C.ink}`,display:"inline-block"}}/>stable fixed-point branch</span><span className="flex items-center gap-1.5"><i style={{width:14,height:0,borderTop:`1.5px dashed ${C.unstTT}`,display:"inline-block"}}/>unstable fixed-point branch</span><span className="flex items-center gap-1.5"><i style={{width:14,height:0,borderTop:`1px dotted ${C.flip}`,display:"inline-block"}}/>λ = −1 crossing</span><span>showing up to 320 samples from the last {bifData.keep} of {bifData.total} iterates at each of {bifSamples} parameter values</span></div></Panel></CollapsibleBifurcation></div>
        <div className="flex flex-wrap gap-3"><Panel title="Orbit: mₙ, pₙ, Mₙ, τₙ" hint={fate||"no initial condition"}><div ref={seriesC.wrap}><canvas ref={seriesC.cvs} style={{display:"block"}}/></div><div className="flex flex-wrap gap-4 px-3 py-2" style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,borderTop:`1px solid ${C.rule}`}}><span className="flex items-center gap-1.5"><i style={{width:11,height:11,background:C.lact,border:`1px solid ${C.rule}`,display:"inline-block"}}/>target-met session</span><span className="flex items-center gap-1.5"><i style={{width:11,height:11,background:C.empt,border:`1px solid ${C.rule}`,display:"inline-block"}}/>emptying session</span><span className="flex items-center gap-1.5"><i style={{width:8,height:8,background:C.traj,display:"inline-block",borderRadius:99}}/>target-met iterate</span><span className="flex items-center gap-1.5"><i style={{width:8,height:8,background:C.mark,display:"inline-block",borderRadius:99}}/>emptying iterate</span></div></Panel>
          <FixedPoints points={[{label:"z₂",m:M.mstar,p:M.pstar,M:M.type==="emptied"?M.Memp:M.Mt}]} />
        </div>
      </div>
    </div>
  </div>;
}


/* ---------------------------------------------------------------- shell */
function LactationExplorer() {
  const [tab, setTab] = useState("target");
  const [P, setP] = useState({ gam: 0.1, dlt: 0.05, amax: 30, muM: 110, muTau: 0.1, tauEmp: 1/3, resp: "sig" });
  const [T, setT] = useState(5);
  const matrixRows = [
    {label:"Amount-dependent",short:"Amount",buttons:[["target","Set target"],["duration","Set duration"]]},
    {label:"Duration-dependent",short:"Duration",buttons:[["tautarget","Set target"],["tauduration","Set duration"]]},
  ];
  return (
    <div className="w-full min-h-screen p-4" style={{ background: C.paper }}>
      <style>{`input[type=range]:focus-visible{outline:2px solid ${C.traj};outline-offset:2px}
        canvas:focus-visible{outline:2px solid ${C.traj}}`}</style>
      <section className="mb-4 rounded-sm p-3" style={{background:C.panel,border:`1px solid ${C.rule}`}}>
        <div className="flex flex-wrap gap-5 items-start">
          <nav aria-label="Select stimulation mechanism and removal protocol" style={{minWidth:338,maxWidth:470,flex:"0 1 420px"}}>
            <div style={{fontFamily:SANS,fontSize:15,fontWeight:600,color:C.ink,marginBottom:7}}>
              Select model
            </div>
            <div style={{display:"grid",gridTemplateColumns:"126px minmax(98px,1fr) minmax(98px,1fr)",gap:5,alignItems:"stretch"}}>
              <div style={{gridColumn:"2 / span 2",fontFamily:SANS,fontSize:10.5,
                fontWeight:600,color:C.ink2,textAlign:"center"}}>Removal protocol →</div>
              <div style={{fontFamily:SANS,fontSize:10,color:C.ink3,alignSelf:"end"}}>Stimulation mechanism ↓</div>
              <div style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,textAlign:"center",alignSelf:"end"}}>Set target</div>
              <div style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,textAlign:"center",alignSelf:"end"}}>Set duration</div>
              {matrixRows.map((row,i)=><React.Fragment key={row.label}>
                <div style={{fontFamily:SANS,fontSize:10.5,color:C.ink2,alignSelf:"center"}}>{row.label}</div>
                {row.buttons.map(([key,label])=><button key={key} onClick={()=>setTab(key)}
                  aria-label={`${row.label} stimulation, ${label.toLowerCase()} removal`}
                  aria-pressed={tab===key}
                  className="px-2 py-2 rounded-sm"
                  style={{fontFamily:SANS,fontSize:10.5,fontWeight:tab===key?600:400,
                    background:tab===key?C.ink:C.paper,color:tab===key?C.panel:C.ink2,
                    border:`1px solid ${tab===key?C.ink:C.rule}`}}>
                    {`${row.short} · ${label.toLowerCase()}`}
                  </button>)}
              </React.Fragment>)}
            </div>

          </nav>
          <div className="flex-1 min-w-[320px]">
            <h1 style={{fontFamily:SERIF,fontSize:19,fontWeight:600,color:C.ink,marginBottom:5}}>
              Interactive Lactation Explorer
            </h1>
            <p style={{fontFamily:SANS,fontSize:12,color:C.ink2,lineHeight:1.55}}>
              Explore how repeated waiting and milk-removal sessions affect milk accumulation and production.
              Choose a stimulation mechanism and a removal protocol on the left, then select a schedule in the
              schedule plot and an initial condition in the phase plane. The removal law is
              <span style={{fontFamily:MONO}}> h(τ) = max(1−τ/τ_emp, 0)</span>.
            </p>
            <p style={{fontFamily:SANS,fontSize:11,color:C.ink2,lineHeight:1.5,marginTop:6}}>
              <b>Paper:</b> A. Manhart and T. Müller, <i>Milk-removal feedback shapes lactation sustainability
              and daily output in a model of human lactation</i> (manuscript, 2026).
            </p>
            <p style={{fontFamily:SANS,fontSize:11,color:C.ink2,lineHeight:1.5,marginTop:6}}>
              <b>Analytical versus numerical:</b> Fixed points, stability boundaries and steady-state output
              are analytical; orbits, non-constant attractors and numerical optimum traces are approximations.
              In bistable regions, daily output refers to the positive state and need not occur from every initial
              condition. Numerical results depend on resolution. These are not feeding or pumping recommendations.
            </p>
          </div>
        </div>
      </section>
      {tab === "target" ? <TargetView P={P} setP={setP} T={T} setT={setT} />
       : tab === "duration" ? <DurationView P={P} setP={setP} T={T} setT={setT} />
       : tab === "tautarget" ? <TauStimTargetView P={P} setP={setP} T={T} setT={setT} />
       : <TauStimDurationView P={P} setP={setP} T={T} setT={setT} />}
    </div>
  );
}

function hex(h) {
  return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
}

class ExplorerErrorBoundary extends React.Component {
  constructor(props) { super(props); this.state = { error: null }; }
  static getDerivedStateFromError(error) { return { error }; }
  componentDidCatch(error, info) { console.error("Explorer rendering error", error, info); }
  render() {
    if (this.state.error) return <div style={{padding:24,fontFamily:"system-ui",color:"#822"}}>
      <h2>Explorer could not render</h2>
      <pre style={{whiteSpace:"pre-wrap"}}>{String(this.state.error)}</pre>
    </div>;
    return this.props.children;
  }
}
createRoot(document.getElementById("root")).render(<ExplorerErrorBoundary><LactationExplorer /></ExplorerErrorBoundary>);
