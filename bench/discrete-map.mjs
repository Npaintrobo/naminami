// 円形シートを離散ピンで駆動したときの「写像できる上限」を調べる。
//
// モデル:
//   - 正方格子のピンのうち、円板内にあるものだけを使う
//   - 各ピンに進行波の鉛直変位を与える
//   - ピン間の無荷重シート形状は、曲げエネルギーを最小にする薄板スプラインで補間する
//   - 上面の面内変位は Kirchhoff-Love の u_top = -(t/2) grad(w) から評価する
//
// これは無荷重・準静的な上界評価であり、溝接触や超弾性を解くFEMではない。
// reg はピン追従性を落とす無次元の感度パラメータで、実材料との校正前は物性値ではない。

function lu(A, n) {
  const piv = Int32Array.from({ length: n }, (_, i) => i);
  for (let k = 0; k < n; k += 1) {
    let p = k, mx = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i += 1) {
      const v = Math.abs(A[i * n + k]);
      if (v > mx) { mx = v; p = i; }
    }
    if (!(mx > 1e-13)) throw new Error(`singular pivot ${k}: ${mx}`);
    if (p !== k) {
      for (let j = 0; j < n; j += 1) {
        const v = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = v;
      }
      const v = piv[k]; piv[k] = piv[p]; piv[p] = v;
    }
    const d = A[k * n + k];
    for (let i = k + 1; i < n; i += 1) {
      const f = A[i * n + k] / d;
      A[i * n + k] = f;
      for (let j = k + 1; j < n; j += 1) A[i * n + j] -= f * A[k * n + j];
    }
  }
  return piv;
}

function luSolve(A, piv, b, n) {
  const y = new Float64Array(n), x = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    let s = b[piv[i]];
    for (let j = 0; j < i; j += 1) s -= A[i * n + j] * y[j];
    y[i] = s;
  }
  for (let i = n - 1; i >= 0; i -= 1) {
    let s = y[i];
    for (let j = i + 1; j < n; j += 1) s -= A[i * n + j] * x[j];
    x[i] = s / A[i * n + i];
  }
  return x;
}

function phi(dx, dy) {
  const r2 = dx * dx + dy * dy;
  return r2 === 0 ? 0 : 0.5 * r2 * Math.log(r2);
}

function gradPhi(dx, dy) {
  const r2 = dx * dx + dy * dy;
  if (r2 === 0) return [0, 0];
  const c = Math.log(r2) + 1;
  return [c * dx, c * dy];
}

function pinSet(N) {
  const R = N / 2, pts = [];
  for (let j = 0; j < N; j += 1) for (let i = 0; i < N; i += 1) {
    const x = i - (N - 1) / 2, y = j - (N - 1) / 2;
    if (Math.hypot(x, y) <= R) pts.push([x, y]);
  }
  return { R, pts };
}

function fitPair({ N = 12, samplesPerWave = 6, psi = 0, reg = 0 } = {}) {
  const { R, pts } = pinSet(N), n = pts.length, dim = n + 3;
  const A = new Float64Array(dim * dim);
  for (let i = 0; i < n; i += 1) for (let j = 0; j < n; j += 1)
    A[i * dim + j] = phi(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]);
  for (let i = 0; i < n; i += 1) {
    A[i * dim + i] += reg;
    const [x, y] = pts[i];
    for (const [q, v] of [[n, 1], [n + 1, x], [n + 2, y]]) {
      A[i * dim + q] = v; A[q * dim + i] = v;
    }
  }
  const factored = new Float64Array(A), piv = lu(factored, dim);
  const k = 2 * Math.PI / samplesPerWave, d = [Math.cos(psi), Math.sin(psi)];
  const bc = new Float64Array(dim), bs = new Float64Array(dim);
  for (let i = 0; i < n; i += 1) {
    const s = pts[i][0] * d[0] + pts[i][1] * d[1];
    bc[i] = Math.cos(k * s); bs[i] = Math.sin(k * s);
  }
  return { R, pts, k, d, cc: luSolve(factored, piv, bc, dim), cs: luSolve(factored, piv, bs, dim) };
}

function evaluate(f, c, x, y) {
  const n = f.pts.length;
  let z = c[n] + c[n + 1] * x + c[n + 2] * y, gx = c[n + 1], gy = c[n + 2];
  for (let i = 0; i < n; i += 1) {
    const dx = x - f.pts[i][0], dy = y - f.pts[i][1], g = gradPhi(dx, dy);
    z += c[i] * phi(dx, dy); gx += c[i] * g[0]; gy += c[i] * g[1];
  }
  return { z, gx, gy };
}

function quantile(a, q) {
  const b = [...a].sort((x, y) => x - y);
  return b[Math.min(b.length - 1, Math.max(0, Math.floor(q * (b.length - 1))))];
}

export function mappingMetrics({ N = 12, samplesPerWave = 6, psi = 0, reg = 0,
                                 marginPitch = 1, evalStepPitch = 0.25 } = {}) {
  const f = fitPair({ N, samplesPerWave, psi, reg });
  const gain = [], leakage = [], vertical = [];
  const perp = [-f.d[1], f.d[0]];
  for (let y = -f.R + marginPitch; y <= f.R - marginPitch + 1e-9; y += evalStepPitch)
    for (let x = -f.R + marginPitch; x <= f.R - marginPitch + 1e-9; x += evalStepPitch) {
      if (Math.hypot(x, y) > f.R - marginPitch) continue;
      const c = evaluate(f, f.cc, x, y), s = evaluate(f, f.cs, x, y);
      const pc = c.gx * f.d[0] + c.gy * f.d[1], ps = s.gx * f.d[0] + s.gy * f.d[1];
      const qc = c.gx * perp[0] + c.gy * perp[1], qs = s.gx * perp[0] + s.gy * perp[1];
      gain.push(Math.hypot(pc, ps) / f.k);
      leakage.push(Math.hypot(qc, qs) / f.k);
      vertical.push(Math.hypot(c.z, s.z));
    }
  const g50 = quantile(gain, 0.5), g05 = quantile(gain, 0.05), g95 = quantile(gain, 0.95);
  return {
    N, samplesPerWave, psiDeg: psi * 180 / Math.PI, reg, pinCount: f.pts.length,
    gainMedian: g50, gainP05: g05, gainP95: g95, gainSpread: (g95 - g05) / g50,
    leakageP95: quantile(leakage, 0.95), verticalP05: quantile(vertical, 0.05),
  };
}

export function directionEnvelope({ directionsDeg = [0, 15, 30, 45], ...opts } = {}) {
  const rows = directionsDeg.map(deg => mappingMetrics({ ...opts, psi: deg * Math.PI / 180 }));
  return {
    rows,
    gainMedianMin: Math.min(...rows.map(r => r.gainMedian)),
    gainP05Min: Math.min(...rows.map(r => r.gainP05)),
    gainSpreadMax: Math.max(...rows.map(r => r.gainSpread)),
    leakageP95Max: Math.max(...rows.map(r => r.leakageP95)),
    verticalP05Min: Math.min(...rows.map(r => r.verticalP05)),
  };
}

export function selfCheck() {
  const f = fitPair({ N: 12, samplesPerWave: 6, psi: 0.37, reg: 0 });
  let residual = 0;
  for (let i = 0; i < f.pts.length; i += 1) {
    const [x, y] = f.pts[i], s0 = x * f.d[0] + y * f.d[1];
    residual = Math.max(residual,
      Math.abs(evaluate(f, f.cc, x, y).z - Math.cos(f.k * s0)),
      Math.abs(evaluate(f, f.cs, x, y).z - Math.sin(f.k * s0)));
  }
  const a = mappingMetrics({ psi: 0 }), b = mappingMetrics({ psi: Math.PI / 2 });
  const symmetry = Math.abs(a.gainMedian - b.gainMedian);
  const coarse = mappingMetrics({ samplesPerWave: 6, psi: Math.PI / 6, evalStepPitch: 0.5 });
  const fine = mappingMetrics({ samplesPerWave: 6, psi: Math.PI / 6, evalStepPitch: 0.125 });
  const quadrature = Math.abs(coarse.gainMedian - fine.gainMedian);
  if (residual > 1e-8) throw new Error(`pin interpolation residual ${residual}`);
  if (symmetry > 1e-8) throw new Error(`90-degree symmetry residual ${symmetry}`);
  if (quadrature > 0.01) throw new Error(`evaluation-grid residual ${quadrature}`);
  return { residual, symmetry, quadrature };
}

if (import.meta.url === new URL(`file://${process.argv[1]}`).href) {
  const check = selfCheck();
  console.log(`self-check: pin residual=${check.residual.toExponential(2)}, ` +
              `90deg symmetry=${check.symmetry.toExponential(2)}, ` +
              `grid residual=${check.quadrature.toExponential(2)}`);
  console.log('\n12×12 unit, circular sheet; worst of ψ=0/15/30/45°');
  console.log('λ/p  reg  gain50(min) gain05(min) spread(max) leak95(max) z05(min)');
  for (const samplesPerWave of [2, 3, 4, 6, 8]) for (const reg of [0, 1]) {
    const r = directionEnvelope({ N: 12, samplesPerWave, reg });
    console.log(`${String(samplesPerWave).padStart(3)}  ${String(reg).padStart(3)}      ` +
      `${r.gainMedianMin.toFixed(3)}       ${r.gainP05Min.toFixed(3)}       ` +
      `${r.gainSpreadMax.toFixed(2)}        ${r.leakageP95Max.toFixed(3)}      ${r.verticalP05Min.toFixed(3)}`);
  }
}
