// 上下変位だけを与えて、ゴム表面に面内運動が「出てくる」か。
//
// docs/membrane-wave.md は薄板理論から u̇ = −h·W·k·Ω·cos φ を導いて、それを
// そのままコードへ入力していた。つまり答えを仮定していた。ここでは
//
//   入力は下面の鉛直変位だけ。面内は一切与えない。
//
// 平面ひずみの線形弾性を有限要素で解き、上面の粒子が楕円を描くか、その面内振幅が
// 薄板の予測 a_eff = (t/2)·W·k と合うかを見る。
//
// 波は w = W·cos(kx − ωt) = W[cos(kx)·cos(ωt) + sin(kx)·sin(ωt)] と分解できるので、
// 静的な問題を2つ解けば時間発展は解析的に出る。各節点の軌道はその2解が張る楕円。
//
//   node bench/elastic.mjs

/* ── 平面ひずみ 線形弾性、双一次四辺形、2×2 Gauss ────────────── */
function elementK(w, h, E, nu) {
  /* 平面ひずみの構成則 */
  const c = E / ((1 + nu) * (1 - 2 * nu));
  const D = [[c * (1 - nu), c * nu, 0], [c * nu, c * (1 - nu), 0], [0, 0, c * (1 - 2 * nu) / 2]];
  const K = new Float64Array(64);
  const g = 1 / Math.sqrt(3);
  for (const xi of [-g, g]) for (const et of [-g, g]) {
    /* 形状関数の自然座標微分 */
    const dN = [[-(1 - et) / 4, -(1 - xi) / 4], [(1 - et) / 4, -(1 + xi) / 4],
                [(1 + et) / 4, (1 + xi) / 4], [-(1 + et) / 4, (1 - xi) / 4]];
    const J = [w / 2, h / 2];                      // 長方形要素なのでヤコビアンは対角
    const detJ = J[0] * J[1];
    const B = [new Float64Array(8), new Float64Array(8), new Float64Array(8)];
    for (let a = 0; a < 4; a += 1) {
      const dx = dN[a][0] / J[0], dz = dN[a][1] / J[1];
      B[0][a * 2] = dx;
      B[1][a * 2 + 1] = dz;
      B[2][a * 2] = dz; B[2][a * 2 + 1] = dx;
    }
    for (let i = 0; i < 8; i += 1) for (let j = 0; j < 8; j += 1) {
      let s = 0;
      for (let m = 0; m < 3; m += 1) for (let n = 0; n < 3; n += 1) s += B[m][i] * D[m][n] * B[n][j];
      K[i * 8 + j] += s * detJ;
    }
  }
  return K;
}

/* 密行列 LU（部分ピボット）。1回分解して2回解く。 */
function lu(A, n) {
  const piv = new Int32Array(n);
  for (let i = 0; i < n; i += 1) piv[i] = i;
  for (let k = 0; k < n; k += 1) {
    let p = k, mx = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i += 1) { const v = Math.abs(A[i * n + k]); if (v > mx) { mx = v; p = i; } }
    if (p !== k) {
      for (let j = 0; j < n; j += 1) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
      const t = piv[k]; piv[k] = piv[p]; piv[p] = t;
    }
    const d = A[k * n + k];
    for (let i = k + 1; i < n; i += 1) {
      const f = A[i * n + k] / d;
      A[i * n + k] = f;
      if (f === 0) continue;
      for (let j = k + 1; j < n; j += 1) A[i * n + j] -= f * A[k * n + j];
    }
  }
  return piv;
}

function luSolve(A, piv, b, n) {
  const y = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    let s = b[piv[i]];
    for (let j = 0; j < i; j += 1) s -= A[i * n + j] * y[j];
    y[i] = s;
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i -= 1) {
    let s = y[i];
    for (let j = i + 1; j < n; j += 1) s -= A[i * n + j] * x[j];
    x[i] = s / A[i * n + i];
  }
  return x;
}

/* ── 下面に鉛直変位だけを与えて解く ───────────────────── */
function solve({ lam = 100, thick = 6, W = 3, E = 5, nu = 0.45,
                 nx = 60, nz = 6, bottomLateral = 'fixed', nWave = 2 } = {}) {
  const L = lam * nWave;
  const dx = L / nx, dz = thick / nz;
  const nnx = nx + 1, nnz = nz + 1;
  const nn = nnx * nnz, ndof = nn * 2;
  const nid = (i, j) => j * nnx + i;
  const k = 2 * Math.PI / lam;

  /* 全体剛性 */
  const Kg = new Float64Array(ndof * ndof);
  const Ke = elementK(dx, dz, E, nu);
  for (let j = 0; j < nz; j += 1) for (let i = 0; i < nx; i += 1) {
    const nd = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)];
    for (let a = 0; a < 4; a += 1) for (let b = 0; b < 4; b += 1)
      for (let p = 0; p < 2; p += 1) for (let q = 0; q < 2; q += 1)
        Kg[(nd[a] * 2 + p) * ndof + nd[b] * 2 + q] += Ke[(a * 2 + p) * 8 + b * 2 + q];
  }

  /* 拘束: 下面（j=0）の鉛直は規定。水平は溝の拘束の有無で切り替える。
     左右端は x 変位を止めて端部の剛体並進を抑える。 */
  const fixed = new Uint8Array(ndof);
  const val = [new Float64Array(ndof), new Float64Array(ndof)];   // [cos成分, sin成分]
  for (let i = 0; i < nnx; i += 1) {
    const n0 = nid(i, 0), x = i * dx;
    fixed[n0 * 2 + 1] = 1;
    val[0][n0 * 2 + 1] = W * Math.cos(k * x);
    val[1][n0 * 2 + 1] = W * Math.sin(k * x);
    if (bottomLateral === 'fixed') fixed[n0 * 2] = 1;             // 溝が面内も拘束する場合
  }
  if (bottomLateral !== 'fixed') {
    for (let j = 0; j < nnz; j += 1) { fixed[nid(0, j) * 2] = 1; fixed[nid(nx, j) * 2] = 1; }
  }

  const free = [];
  for (let d = 0; d < ndof; d += 1) if (!fixed[d]) free.push(d);
  const nf = free.length;
  const A = new Float64Array(nf * nf);
  for (let a = 0; a < nf; a += 1) for (let b = 0; b < nf; b += 1) A[a * nf + b] = Kg[free[a] * ndof + free[b]];
  const piv = lu(A, nf);

  const out = [];
  for (const v of val) {
    const rhs = new Float64Array(nf);
    for (let a = 0; a < nf; a += 1) {
      let s = 0;
      for (let d = 0; d < ndof; d += 1) if (fixed[d] && v[d] !== 0) s -= Kg[free[a] * ndof + d] * v[d];
      rhs[a] = s;
    }
    const x = luSolve(A, piv, rhs, nf);
    const u = new Float64Array(ndof);
    for (let d = 0; d < ndof; d += 1) if (fixed[d]) u[d] = v[d];
    for (let a = 0; a < nf; a += 1) u[free[a]] = x[a];
    out.push(u);
  }

  /* 上面（j=nz）の中央1波長ぶんで、楕円の半径を測る */
  const [uc, us] = out;
  let ax = 0, az = 0, cnt = 0;
  const i0 = Math.floor(nnx * (nWave - 1) / (2 * nWave)), i1 = i0 + Math.floor(nnx / nWave);
  for (let i = i0; i < i1; i += 1) {
    const n = nid(i, nz);
    ax += Math.hypot(uc[n * 2], us[n * 2]);
    az += Math.hypot(uc[n * 2 + 1], us[n * 2 + 1]);
    cnt += 1;
  }
  return { ax: ax / cnt, az: az / cnt, pred: (thick / 2) * W * k, W, thick, lam, k, nf };
}

const R = r => `面内${r.ax.toFixed(3)}mm 上下${r.az.toFixed(3)}mm  薄板予測${r.pred.toFixed(3)}mm  比${(r.ax / r.pred).toFixed(2)}`;

console.log('■ 入力は下面の鉛直変位だけ。面内は一切与えていない。\n');
console.log('  上面の粒子が描く楕円の半径を測る。薄板理論の予測は a_eff = (t/2)·W·k。');
console.log('  E=5MPa（Shore A50 相当）、ν=0.45、平面ひずみ、2波長のうち中央1波長で平均。\n');

console.log('  下面の面内拘束   λ[mm]  厚[mm]  W[mm]   測定結果');
for (const bl of ['fixed', 'free']) {
  for (const [lam, thick, W] of [[100, 6, 3], [100, 6, 5], [150, 6, 5], [100, 10, 5]]) {
    const r = solve({ lam, thick, W, bottomLateral: bl });
    console.log(`  ${(bl === 'fixed' ? '溝で拘束' : '自由（滑る）').padEnd(14)} ${String(lam).padStart(4)}  ${String(thick).padStart(5)}  ${String(W).padStart(4)}   ${R(r)}`);
  }
}

console.log('\n■ 収束確認（λ=100, 厚6, W=5, 溝で拘束）');
console.log('  要素分割      面内[mm]  上下[mm]   比');
for (const [nx, nz] of [[40, 4], [60, 6], [80, 8], [100, 10]]) {
  const r = solve({ lam: 100, thick: 6, W: 5, nx, nz });
  console.log(`  ${String(nx).padStart(3)}×${String(nz).padStart(2)}（自由度${String(r.nf).padStart(4)}）  ` +
    `${r.ax.toFixed(3).padStart(7)}  ${r.az.toFixed(3).padStart(7)}  ${(r.ax / r.pred).toFixed(2)}`);
}

console.log('\n■ ポアソン比（ゴムは非圧縮に近い。双一次要素はロッキングする）');
console.log('  ν       面内[mm]  上下[mm]   比      備考');
for (const nu of [0.3, 0.40, 0.45, 0.48, 0.49]) {
  const r = solve({ lam: 100, thick: 6, W: 5, nu });
  const note = nu >= 0.48 ? 'ロッキングの疑い' : '';
  console.log(`  ${nu.toFixed(2)}    ${r.ax.toFixed(3).padStart(7)}  ${r.az.toFixed(3).padStart(7)}  ` +
    `${(r.ax / r.pred).toFixed(2)}    ${note}`);
}

console.log('\n■ 端部の影響（波長数を増やして中央で測る）');
console.log('  波長数   面内[mm]   比');
for (const nWave of [2, 3, 4]) {
  const r = solve({ lam: 100, thick: 6, W: 5, nx: 30 * nWave, nWave });
  console.log(`  ${nWave}      ${r.ax.toFixed(3).padStart(7)}   ${(r.ax / r.pred).toFixed(2)}`);
}

console.log(`
■ 分かったこと

  ・**面内運動は本当に出てくる。** 入力は下面の鉛直変位だけで、面内は一切与えていない。
    docs/membrane-wave.md は答えを仮定していたが、仮定した向きは正しかった。
  ・**下面を溝で面内拘束すると、薄板予測の 1.3〜1.8倍出る。** 下面が滑れるときは
    薄板理論とほぼ一致（比 0.93〜0.99）。溝は両方向拘束だけでなく、面内利得も稼いでいる。

■ ただし注意

  ・同心円溝が拘束するのは**半径方向だけ**。波の進む向きが半径方向と一致する場所では
    「溝で拘束」、直交する場所では「自由」に近い。**利得は場所によって 1.0〜1.7 に変わる。**
    つまり円板上で駆動が不均一になる。これは membrane.mjs にも入っていない。
  ・線形弾性・微小変形。ゴムの超弾性、W/λ が大きいときの幾何非線形は入っていない。
  ・平面ひずみの2次元断面。円板の3次元性（半径方向と円周方向で曲げ剛性が違うこと、
    docs/groove-coupling.md §5）は入っていない。
  ・接地面との接触を解いていない。上面は自由表面のまま。荷を載せると楕円は潰れる。`);
