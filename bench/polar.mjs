// 和差位相式・偏極操舵波動面を、案A（直交2ラティス）と同じ物理コアで測る。
//
//   θ_i = k·r_i − q          搬送波位相。どの接触子が荷重を受けるかを決める
//   z_i = −b·cos θ_i         上下（歯は下向き）
//   p_i =  a·sin θ_i·[cosψ, sinψ]   面内。ψ が偏極角＝推力方向
//
//   ṗ_i = −a·q̇·cos θ_i·[cosψ, sinψ] + a·ψ̇·sin θ_i·[−sinψ, cosψ]
//          ↑ 山（θ=0）で a·q̇                ↑ 操舵中だけ出る横成分。山では消える
//
//   q = (φ_A+φ_B)/2,  ψ = (φ_A−φ_B)/2
//
// このファイルが数値の出どころです。3D（sim3d-polar.html）は同じ式・同じ整定則を
// 使いますがヨーを持たないので、論文値はこちらを使ってください。
//
// 測定条件はすべての試験で共通です:
//   ・接地剛性 K は「静的釣合いさせたときの最大貫入の1周期平均が δp になる」よう二分法で解く
//   ・そこから静的釣合い高さに置いて開始する（自由落下させない）
//   ・回転数は RAMP 秒かけて立ち上げる（起動衝撃を測定に混ぜない）
//   ・暖機を捨てた区間の平均を出す
//
//   node bench/polar.mjs

const N = 12, P = 50, AMP = 6, RPM = 300;
const BREF = { v: 6 };                     // 上下振幅 b（試験で振る）
const W = RPM * Math.PI / 30;
const TIP = AMP * W;                       // 山での面内速度 = 案Aの r·ω と同じ 188 mm/s
const MU = 0.6, G = 9810, VEPS = 6, ACC = 1000;
const MASS = 20, HCG = 100;
const DP = { v: 0.3 };                     // パッドの沈み込み δp [mm]
const RAMP = 0.4;                          // 回転数の立ち上げ [s]
const DT = 5e-5, WARM = 1.2, SPAN = 2.0;
const TAU = 2 * Math.PI;

const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const wrap = v => { let d = v % TAU; if (d > Math.PI) d -= TAU; if (d < -Math.PI) d += TAU; return d; };
const rampOf = t => Math.min(1, t / RAMP);

/* ── 接触子の配置と位相場 ───────────────────────── */
function lattice(dphiDeg, { shuffle = false, seed = 7, kdir = 'diag', tol = 0 } = {}) {
  const d = dphiDeg * Math.PI / 180;
  const pads = [];
  const idx = (i, k) => (kdir === 'diag' ? i + k : i);
  for (let i = 0; i < N; i += 1) for (let k = 0; k < N; k += 1) {
    pads.push({ x: (i - (N - 1) / 2) * P, z: (k - (N - 1) / 2) * P,
                ph: idx(i, k) * d, xy: (i + k) % 2, n: 0, prev: 0 });
  }
  if (shuffle) {
    let s = seed;
    const rnd = () => (s = (s * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    const ph = pads.map(p => p.ph);
    for (let i = ph.length - 1; i > 0; i -= 1) { const j = Math.floor(rnd() * (i + 1)); [ph[i], ph[j]] = [ph[j], ph[i]]; }
    pads.forEach((p, i) => { p.ph = ph[i]; });
  }
  /* 組立公差。足ごとに固定の高さ誤差を与える（正規分布、σ = tol [mm]）。
     接地窓は δp で決まるので、公差が δp に対して無視できないと
     「一部の足だけが常に荷重を持つ」状態になる。 */
  if (tol > 0) {
    let st = seed * 7919 + 13;
    const rnd = () => (st = Math.imul(st, 1103515245) + 12345 & 0x7fffffff) / 0x7fffffff;
    for (const p of pads) {
      const u = Math.max(1e-12, rnd()), v = rnd();
      p.dz = tol * Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * v);   // Box–Muller
    }
  } else {
    for (const p of pads) p.dz = 0;
  }
  return pads;
}

/* ── 接地剛性の整定 ─────────────────────────────
   「δp」は静的に置いたときの最大貫入と定義する。位相によって振れるので1周期で平均し、
   その平均が δp になる K を二分法で解く。以前は総貫入量から一発で決めていて、
   実際に釣り合わせた貫入（0.12mm）と設定値（0.3mm）が合っていなかった。 */
/* 高さだけの釣合い。3自由度解の初期値に使う。 */
function heaveOnly(ys, K, weight) {
  let lo = -Math.max(...ys.map(Math.abs)) - 30, hi = -Math.min(...ys) + 1;
  for (let i = 0; i < 60; i += 1) {
    const mid = (lo + hi) / 2;
    let f = 0;
    for (const y of ys) f += Math.max(0, -(mid + y)) * K;
    if (f > weight) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

/* 静的平衡を 高さ Y・ピッチ tx・ロール tz で連立して解く。

   貫入 u_i = −(Y − tx·z_i + tz·x_i + y_i)、法線力 N_i = K·max(0,u_i) のとき、
   釣合い ΣN = Mg, ΣN·z = 0, ΣN·x = 0 は、弾性＋重力のポテンシャル

     E(Y,tx,tz) = Σ (K/2)·max(0,u_i)² + Mg·Y

   の停留点そのものです。E は u が引数のアフィン関数なので**凸**で、
   減衰Newton＋直線探索なら必ず収束します（接地が1〜2点しかない配置でも、
   より多くの足が当たるまで傾いていく）。

   以前は高さだけを解いて姿勢を 0 に固定していたので、「傾けば釣り合う」配置まで
   支持不能と誤判定していました。素の Newton は接地が少ない条件で発散しました。 */
function equilibrium(ys, K, weight, xs, zs) {
  if (!xs) return { Y: heaveOnly(ys, K, weight), tx: 0, tz: 0, n: 0, ok: true };
  const N = ys.length;
  let Y = heaveOnly(ys, K, weight), tx = 0, tz = 0;

  const energy = (Y_, tx_, tz_) => {
    let e = weight * Y_;
    for (let i = 0; i < N; i += 1) {
      const u = -(Y_ - tx_ * zs[i] + tz_ * xs[i] + ys[i]);
      if (u > 0) e += 0.5 * K * u * u;
    }
    return e;
  };

  for (let it = 0; it < 200; it += 1) {
    let g1 = weight, g2 = 0, g3 = 0;
    let h11 = 0, h12 = 0, h13 = 0, h22 = 0, h23 = 0, h33 = 0, cnt = 0;
    for (let i = 0; i < N; i += 1) {
      const u = -(Y - tx * zs[i] + tz * xs[i] + ys[i]);
      if (u <= 0) continue;
      cnt += 1;
      const f = K * u, x = xs[i], z = zs[i];
      g1 -= f; g2 += f * z; g3 -= f * x;
      h11 += K;        h12 += -K * z;     h13 += K * x;
      h22 += K * z * z; h23 += -K * x * z; h33 += K * x * x;
    }
    const gn = Math.hypot(g1, g2, g3);
    if (gn < weight * 1e-9) break;

    /* 接地が少ないとヘッセ行列が特異になるので正則化する */
    const lam = Math.max(1e-6, 1e-6 * (h11 + h22 + h33));
    const a11 = h11 + lam, a22 = h22 + lam, a33 = h33 + lam;
    const det = a11 * (a22 * a33 - h23 * h23) - h12 * (h12 * a33 - h23 * h13) + h13 * (h12 * h23 - a22 * h13);
    if (!isFinite(det) || Math.abs(det) < 1e-30) { Y -= 0.1; continue; }
    const b1 = -g1, b2 = -g2, b3 = -g3;
    const dY  = (b1 * (a22 * a33 - h23 * h23) - h12 * (b2 * a33 - h23 * b3) + h13 * (b2 * h23 - a22 * b3)) / det;
    const dtx = (a11 * (b2 * a33 - h23 * b3) - b1 * (h12 * a33 - h23 * h13) + h13 * (h12 * b3 - b2 * h13)) / det;
    const dtz = (a11 * (a22 * b3 - h23 * b2) - h12 * (h12 * b3 - h13 * b2) + b1 * (h12 * h23 - a22 * h13)) / det;
    if (!isFinite(dY) || !isFinite(dtx) || !isFinite(dtz)) break;

    /* 直線探索。凸なので必ず下がる歩幅がある */
    const e0 = energy(Y, tx, tz);
    let step = 1;
    let moved = false;
    for (let k = 0; k < 40; k += 1) {
      if (energy(Y + step * dY, tx + step * dtx, tz + step * dtz) < e0) { moved = true; break; }
      step *= 0.5;
    }
    if (!moved) break;
    Y += step * dY; tx += step * dtx; tz += step * dtz;
    void cnt;
    if (Math.abs(step * dY) < 1e-10 && Math.abs(step * dtx) < 1e-13 && Math.abs(step * dtz) < 1e-13) break;
  }

  let cnt = 0, R = -weight;
  for (let i = 0; i < N; i += 1) {
    const u = -(Y - tx * zs[i] + tz * xs[i] + ys[i]);
    if (u > 0) { cnt += 1; R += K * u; }
  }
  return { Y, tx, tz, n: cnt, ok: isFinite(Y) && Math.abs(R) < weight * 1e-3 };
}

/* 標本数は位相格子と互いに素な素数にする。24点だと Δφ=60°（6階級）と噛み合って
   相対位相が {0,15,30,45}° しか出ず、K が 11.5% ずれていた。 */
const TSAMP = 97;

/* 「δp」は静的平衡での最大貫入と定義し、その1周期平均が δp になる K を二分法で解く。 */
function tune(pads, b = BREF.v, dp = DP.v, mass = MASS) {
  const weight = mass * G / ACC;
  const xs = pads.map(p => p.x), zs = pads.map(p => p.z);
  const phases = [];
  for (let s = 0; s < TSAMP; s += 1) {
    const q = s / TSAMP * TAU;
    phases.push(pads.map(p => -b * Math.cos(p.ph - q)));
  }
  const at = K => {
    let pen = 0, ncs = 0, tilt = 0, bad = 0;
    for (const ys of phases) {
      const e = equilibrium(ys, K, weight, xs, zs);
      if (!e.ok) bad += 1;
      let mx = 0;
      for (let i = 0; i < ys.length; i += 1) {
        const d = -(e.Y - e.tx * zs[i] + e.tz * xs[i] + ys[i]);
        if (d > 0) { ncs += 1; if (d > mx) mx = d; }
      }
      pen += mx;
      tilt = Math.max(tilt, Math.hypot(e.tx, e.tz));
    }
    return { pen: pen / TSAMP, nc: Math.max(1, ncs / TSAMP), tilt: tilt * 180 / Math.PI, bad };
  };
  let klo = 0.05, khi = 20000;
  for (let i = 0; i < 50; i += 1) {
    const km = (klo + khi) / 2;
    if (at(km).pen > dp) klo = km; else khi = km;       // K が大きいほど貫入は小さい
  }
  const K = (klo + khi) / 2;
  const r = at(K);
  const C = 2 * 0.35 * Math.sqrt(r.nc * K * 1000 * mass) / 1000 / r.nc;
  return { K, C, nc: r.nc, pen: r.pen, tilt: r.tilt, bad: r.bad,
           window: Math.acos(clamp(1 - r.pen / b, -1, 1)) * 180 / Math.PI };
}

/* ── 支持の判定 ───────────────────────────────
   「接触が切れない」と「姿勢を保持できる」は別問題。最大位相ギャップが接地窓に
   収まっていても、ある位相で接地点が1〜2点しか無ければピッチ・ロールは支持できない。
   各位相で接地点の凸包を作り、重心の投影（原点）からの余裕を測る。 */
function hull(pts) {
  if (pts.length < 3) return null;
  const s = pts.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cr = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const build = arr => {
    const h = [];
    for (const p of arr) { while (h.length >= 2 && cr(h[h.length - 2], h[h.length - 1], p) <= 0) h.pop(); h.push(p); }
    h.pop(); return h;
  };
  const h = build(s).concat(build(s.slice().reverse()));
  return h.length >= 3 ? h : null;
}

/* 各搬送位相で 高さ・ピッチ・ロール を連立して静的平衡を解き、
   ・解が存在するか
   ・そのとき機体はどれだけ傾く必要があるか
   ・接地点の凸包が重心投影を含むか（水平姿勢を保てるか）
   を分けて返す。以前は水平姿勢に固定した接触集合で凸包だけを見ていたので、
   「傾けば釣り合う」配置を支持不能と誤判定していた（偽陰性）。 */
function supportMargin(pads, b, K) {
  const weight = MASS * G / ACC;
  const xs = pads.map(p => p.x), zs = pads.map(p => p.z);
  let worstMargin = Infinity, worstCnt = 1e9, maxTilt = 0, bad = 0;
  for (let s = 0; s < TSAMP; s += 1) {
    const q = s / TSAMP * TAU;
    const ys = pads.map(p => -b * Math.cos(p.ph - q));
    const e = equilibrium(ys, K, weight, xs, zs);
    if (!e.ok) bad += 1;
    maxTilt = Math.max(maxTilt, Math.hypot(e.tx, e.tz) * 180 / Math.PI);
    const pts = [];
    for (let i = 0; i < pads.length; i += 1) {
      if (-(e.Y - e.tx * zs[i] + e.tz * xs[i] + ys[i]) > 0) pts.push([xs[i], zs[i]]);
    }
    worstCnt = Math.min(worstCnt, pts.length);
    const h = hull(pts);
    if (!h) { worstMargin = -Infinity; continue; }
    let m = Infinity;
    for (let i = 0; i < h.length; i += 1) {
      const a = h[i], c = h[(i + 1) % h.length];
      const len = Math.hypot(c[0] - a[0], c[1] - a[1]) || 1;
      m = Math.min(m, ((c[0] - a[0]) * (0 - a[1]) - (c[1] - a[1]) * (0 - a[0])) / len);
    }
    worstMargin = Math.min(worstMargin, m);
  }
  return { margin: worstMargin, minPts: worstCnt, maxTilt, noSolution: bad };
}

function phaseGap(pads) {
  const u = [...new Set(pads.map(p => ((p.ph % TAU) + TAU) % TAU))].sort((a, b) => a - b);
  let g = 0;
  for (let i = 0; i < u.length; i += 1) {
    const d = (i === u.length - 1 ? u[0] + TAU : u[i + 1]) - u[i];
    if (d > g) g = d;
  }
  return g * 180 / Math.PI;
}

/* ── 6自由度コア ─────────────────────────────── */
function core(pads, drv, { b = BREF.v, K, C } = {}) {
  const L = N * P;
  const Ih = MASS * (L * L / 12 + HCG * HCG), Iy = MASS * (L * L / 6);
  const weight = MASS * G / ACC;

  drv.reset();
  /* 静的釣合いは運動学ごとに高さの式が違うので drv.foot から取り、
     高さ・ピッチ・ロールを連立して解いた姿勢から始める。 */
  const ys = pads.map(p => drv.foot(p)[0]);
  const xs = pads.map(p => p.x), zs = pads.map(p => p.z);
  const eq = equilibrium(ys, K, weight, xs, zs);
  let Y = eq.Y, tx = eq.tx, tz = eq.tz;
  for (let i = 0; i < pads.length; i += 1) {
    pads[i].prev = Math.max(0, -(Y - tx * zs[i] + tz * xs[i] + ys[i]));
  }

  let X = 0, Z = 0, yaw = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  /* 滑り率は「散逸 / 法線力積」で定義する。時間ステップごとに wsum で割って
     時間平均すると、荷重が 0 に近い瞬間も同じ重みで入って過小評価になる。 */
  let dx = 0, dz = 0, path = 0, slipNum = 0, slipDen = 0;
  let wt = 0, nmean = 0, nvar = 0, minC = 1e9, cntSum = 0, tilt2 = 0;
  const T = WARM + SPAN;

  for (let s = 0; s < T / DT; s += 1) {
    const t = s * DT, live = t > WARM;
    drv.step(t, DT);
    const cp = Math.cos(yaw), sp = Math.sin(yaw);
    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0, sl = 0, wsum = 0;

    for (let j = 0; j < pads.length; j += 1) {
      const p = pads[j];
      const [h, ox, oz, ux, uz] = drv.foot(p);
      const ax = p.x + ox, az = p.z + oz;
      const rx = ax * cp - az * sp, rz = ax * sp + az * cp;
      const pen = -(Y - tx * rz + tz * rx + h);
      if (pen <= 0) { p.prev = 0; continue; }
      const Nf = Math.max(0, K * pen + C * (pen - p.prev) / DT);
      p.prev = pen;
      if (Nf <= 0) continue;
      cnt += 1; SN += Nf;
      const wx_ = ux * cp - uz * sp, wz_ = ux * sp + uz * cp;
      /* 接地点は重心より HCG 下。機体の角速度 ω=(wx, vp, wz) がそこに水平速度 ω×r を作る */
      const gx = vX + vp * rz + wz * HCG + wx_;
      const gz = vZ - vp * rx - wx * HCG + wz_;
      const m = Math.max(Math.hypot(gx, gz), VEPS);
      const fx = -MU * Nf * gx / m, fz = -MU * Nf * gz / m;
      Fx += fx; Fz += fz;
      Ty += rz * fx - rx * fz;
      Tx += -HCG * fz - rz * Nf;
      Tz += rx * Nf + HCG * fx;
      sl += Nf * Math.hypot(gx, gz); wsum += Nf;
    }

    vX += ACC * Fx / MASS * DT; vZ += ACC * Fz / MASS * DT;
    vY += (ACC * SN / MASS - G) * DT;
    vp += ACC * Ty / Iy * DT; wx += ACC * Tx / Ih * DT; wz += ACC * Tz / Ih * DT;
    X += vX * DT; Z += vZ * DT; Y += vY * DT; yaw += vp * DT; tx += wx * DT; tz += wz * DT;

    if (live) {
      dx += vX * DT; dz += vZ * DT; wt += DT;
      path += Math.hypot(vX, vZ) * DT;
      slipNum += sl * DT; slipDen += wsum * DT;
      const ld = SN / weight;
      nmean += ld * DT; nvar += ld * ld * DT;
      cntSum += cnt * DT; tilt2 += (tx * tx + tz * tz) * DT;
      if (cnt < minC) minC = cnt;
    }
  }

  const mn = nmean / wt;
  return {
    /* speed は始終点間の弦長÷時間。旋回すると経路長より短く出るので pathSpeed も返す。 */
    speed: Math.hypot(dx, dz) / wt, pathSpeed: path / wt,
    dir: Math.atan2(dz, dx) * 180 / Math.PI,
    slip: slipDen > 0 ? slipNum / slipDen / TIP : 0, load: mn,
    ripple: Math.sqrt(Math.max(0, nvar / wt - mn * mn)) / Math.max(mn, 1e-9),
    contact: cntSum / wt, minC, tilt: Math.sqrt(tilt2 / wt) * 180 / Math.PI,
    yawRate: yaw / (WARM + SPAN) * 180 / Math.PI,      // 平均ヨー角速度 [°/s]
    psi: drv.psi * 180 / Math.PI
  };
}

/* ── 駆動 ───────────────────────────────────── */
function polarDrive({ shaftA = 1, shaftB = 1, hold = false, target = 0, sweep = null, swirl = null } = {}) {
  let phiA, phiB, q, qd, psi, psid;
  const d = {
    reset() { phiA = 0; phiB = 0; q = 0; qd = 0; psi = 0; psid = 0; d.q = 0; d.psi = 0; },
    step(t, dt) {
      const g = rampOf(t);
      let a = shaftA, b = shaftB;
      if (sweep) { const [ps, pd] = sweep(t); psi = ps; psid = pd; q += W * g * dt; qd = W * g; d.q = q; d.psi = psi; return; }
      if (hold) {
        const base = (shaftA + shaftB) / 2, e = wrap(target - psi);
        const dd = Math.abs(e) < 0.009 ? 0 : clamp(e * 0.6, -1.2, 1.2);
        a = base + dd; b = base - dd;
      }
      const wA = a * W * g, wB = b * W * g;
      phiA += wA * dt; phiB += wB * dt;
      q = (phiA + phiB) / 2; qd = (wA + wB) / 2;
      psi = (phiA - phiB) / 2; psid = (wA - wB) / 2;
      d.q = q; d.psi = psi;
    },
    foot(p) {
      /* swirl: 足ごとに ψ をずらす（ヨーを出せるか調べるため）。
         ψ が全足で共通なら合力は必ず ψ に共線で、ヨーモーメントは出ない。 */
      const pp = swirl ? psi + swirl(p) : psi;
      const th = p.ph - q, c = Math.cos(pp), s = Math.sin(pp);
      const sn = Math.sin(th), cs = Math.cos(th);
      return [-BREF.v * cs + (p.dz || 0), AMP * sn * c, AMP * sn * s,
              -AMP * qd * cs * c - AMP * psid * sn * s,
              -AMP * qd * cs * s + AMP * psid * sn * c];
    },
    q: 0, psi: 0
  };
  return d;
}

function shareDrive(cmd, retract = true) {
  /* 位相は積分して作る。rate×t で作ると、ランプ中は実際の微分が 2倍になり
     接触位置と接触速度が食い違う。 */
  let px = 0, pz = 0, wx = 0, wz = 0;
  const d = {
    reset() { px = 0; pz = 0; wx = 0; wz = 0; d.psi = cmd; },
    step(t, dt) {
      const g = rampOf(t);
      wx = W * g * Math.cos(cmd); wz = W * g * Math.sin(cmd);
      px += wx * dt; pz += wz * dt;
    },
    foot(p) {
      /* 指令が 0 の軸だけを退避させる。以前は速度どうしを比較していて、
         起動時（両方 0）に両ラティスが持ち上がり発散していた。 */
      const liftX = retract && Math.abs(Math.cos(cmd)) < 1e-9 ? 2 * BREF.v : 0;
      const liftZ = retract && Math.abs(Math.sin(cmd)) < 1e-9 ? 2 * BREF.v : 0;
      if (p.xy === 0) { const th = p.ph - px; return [-BREF.v * Math.cos(th) + liftX, 0, 0, -AMP * wx * Math.cos(th), 0]; }
      const th = p.ph - pz; return [-BREF.v * Math.cos(th) + liftZ, 0, 0, 0, -AMP * wz * Math.cos(th)];
    },
    psi: cmd
  };
  return d;
}

export { lattice, tune, core, polarDrive, shareDrive, supportMargin, phaseGap, BREF, DP, TIP, N, P, AMP, W };

/* ================================================================
   卓上（有限の荷）。接触子は地面に固定され、荷がその上を渡る。

   位相は φ(i,k) = (i+k)·Δφ を添字から直接作る。以前は 12×12 の表を mod 12 で
   引いていたので、境界ごとに位相が跳んで（Δφ=15° では 180°）不連続面ができ、
   荷は開始時からそれを跨いでいた。荷も静的釣合いから始めて回転数をランプする。
   ================================================================ */
function tableRun(dphiDeg, { kdir = 'diag', shuffle = false, psi = Math.PI / 4, q0 = 0,
                             b = 1.5, S = 150, M = 0.5, HCG2 = 15, dp = 0.4, Kfix = null } = {}) {
  const d = dphiDeg * Math.PI / 180;
  /* シャッフルは「同じ位相集合を別の場所に配る」でなければ比較にならない。
     階級数 M = 360/Δφ の集合 {m·Δφ} から一様に引く（以前は連続一様乱数で、
     Δφ を変えても同じ位相場になっていた）。 */
  const cls = Math.max(1, Math.round(360 / dphiDeg));
  const hash = (i, k) => {
    let h = (i * 73856093) ^ (k * 19349663);
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    return ((h >>> 0) / 4294967296);
  };
  const at = (i, k) => (shuffle ? Math.floor(hash(i, k) * cls) * d
                                : (kdir === 'diag' ? i + k : i) * d) + q0;

  const Ih = M * (S * S / 12 + HCG2 * HCG2), Iy = M * (S * S / 6);
  const weight = M * G / ACC;
  const TIPPED = 8 * Math.PI / 180;
  const c = Math.cos(psi), sn = Math.sin(psi);

  /* 荷の下に入る接触子の高さ top = b(1+cos θ) から静的釣合い高さを解く */
  const under = (X, Z) => {
    const out = [];
    for (let i = Math.floor((X - S / 2) / P); i <= Math.ceil((X + S / 2) / P); i += 1)
      for (let k = Math.floor((Z - S / 2) / P); k <= Math.ceil((Z + S / 2) / P); k += 1) {
        if (Math.abs(i * P - X) > S / 2 || Math.abs(k * P - Z) > S / 2) continue;
        out.push([i, k]);
      }
    return out;
  };
  const list0 = under(0, 0);
  const K = (() => {
    const meanMax = KK => {
      let acc = 0; const SS = 24;
      for (let s = 0; s < SS; s += 1) {
        const q = s / SS * TAU;
        const tops = list0.map(([i, k]) => b + b * Math.cos(at(i, k) - q));
        let lo = Math.min(...tops) - 30, hi = Math.max(...tops) + 1;
        for (let n = 0; n < 60; n += 1) {
          const mid = (lo + hi) / 2;
          let f = 0; for (const tp of tops) f += Math.max(0, tp - mid) * KK;
          if (f > weight) lo = mid; else hi = mid;
        }
        const Y = (lo + hi) / 2;
        let mx = 0; for (const tp of tops) mx = Math.max(mx, tp - Y);
        acc += mx;
      }
      return acc / SS;
    };
    if (Kfix) return Kfix;
    let klo = 0.001, khi = 500;
    for (let n = 0; n < 50; n += 1) { const km = (klo + khi) / 2; if (meanMax(km) > dp) klo = km; else khi = km; }
    return (klo + khi) / 2;
  })();
  const tops0 = list0.map(([i, k]) => b + b * Math.cos(at(i, k)));
  let lo = Math.min(...tops0) - 30, hi = Math.max(...tops0) + 1;
  for (let n = 0; n < 60; n += 1) {
    const mid = (lo + hi) / 2;
    let f = 0; for (const tp of tops0) f += Math.max(0, tp - mid) * K;
    if (f > weight) lo = mid; else hi = mid;
  }
  const C = 2 * 0.35 * Math.sqrt(Math.max(1, list0.length) * K * 1000 * M) / 1000 / Math.max(1, list0.length);

  let X = 0, Z = 0, Y = (lo + hi) / 2, yaw = 0, tx = 0, tz = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  const prev = new Map();
  for (const [i, k] of list0) prev.set(i * 1e5 + k, Math.max(0, b + b * Math.cos(at(i, k)) - Y));
  let dx = 0, dz = 0, wt = 0, drop = 0, tilt2 = 0;
  const T = WARM + SPAN;

  let qq = 0;                       // 搬送位相。rate×t ではなく積分して作る
  for (let s = 0; s < T / DT; s += 1) {
    const t = s * DT, live = t > WARM, g = rampOf(t);
    qq += W * g * DT;
    if (Math.hypot(tx, tz) > TIPPED) return { fell: t, K };
    const cp = Math.cos(yaw), spz = Math.sin(yaw);
    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0;

    for (const [i, k] of under(X, Z)) {
      const rx = i * P - X, rz = k * P - Z;
      const bx = rx * cp + rz * spz, bz = -rx * spz + rz * cp;
      if (Math.abs(bx) > S / 2 || Math.abs(bz) > S / 2) continue;
      const th = at(i, k) - qq;
      const pen = b + b * Math.cos(th) - (Y - tx * rz + tz * rx);
      const key = i * 1e5 + k;
      if (pen <= 0) { prev.set(key, 0); continue; }
      const Nf = Math.max(0, K * pen + C * (pen - (prev.get(key) || 0)) / DT);
      prev.set(key, pen);
      if (Nf <= 0) continue;
      cnt += 1; SN += Nf;
      const ux = -AMP * W * g * Math.cos(th) * c, uz = -AMP * W * g * Math.cos(th) * sn;
      const gx = vX + vp * rz + wz * HCG2 - ux;
      const gz = vZ - vp * rx - wx * HCG2 - uz;
      const m = Math.max(Math.hypot(gx, gz), VEPS);
      const fx = -MU * Nf * gx / m, fz = -MU * Nf * gz / m;
      Fx += fx; Fz += fz; Ty += rz * fx - rx * fz;
      Tx += -HCG2 * fz - rz * Nf; Tz += rx * Nf + HCG2 * fx;
    }

    vX += ACC * Fx / M * DT; vZ += ACC * Fz / M * DT;
    vY += (ACC * SN / M - G) * DT;
    vp += ACC * Ty / Iy * DT; wx += ACC * Tx / Ih * DT; wz += ACC * Tz / Ih * DT;
    X += vX * DT; Z += vZ * DT; Y += vY * DT; yaw += vp * DT; tx += wx * DT; tz += wz * DT;
    if (Y < -2) { Y = -2; if (vY < 0) vY = 0; }
    if (live) { dx += vX * DT; dz += vZ * DT; wt += DT; tilt2 += (tx * tx + tz * tz) * DT; if (cnt === 0) drop += DT; }
  }
  return { speed: Math.hypot(dx, dz) / wt, dir: Math.atan2(dz, dx) * 180 / Math.PI,
           drop: drop / wt * 100, tilt: Math.sqrt(tilt2 / wt) * 180 / Math.PI, K };
}

export { tableRun };

/* ================================================================
   試験
   ================================================================ */
if (import.meta.url === `file://${process.argv[1]}`) {
  const pct = v => (v * 100).toFixed(1).padStart(5) + '%';
  /* 弦（始終点間）と経路長の両方を出す。旋回や蛇行があると前者は小さく出る。
     滑りは「散逸 / 法線力積」。 */
  const fmt = r => `弦${r.speed.toFixed(0).padStart(4)}/経路${r.pathSpeed.toFixed(0).padStart(4)} mm/s ` +
    `(${(r.pathSpeed / TIP * 100).toFixed(0).padStart(3)}%)  ${r.dir.toFixed(1).padStart(7)}°  ` +
    `滑り${pct(r.slip)}  荷重${r.load.toFixed(2)}±${pct(r.ripple)}  ` +
    `傾き${r.tilt.toFixed(2).padStart(5)}°  接地${r.contact.toFixed(1).padStart(4)}(最小${r.minC})`;

  const CFG = [
    ['Δφ=60° k=斜め', 60, 'diag'], ['Δφ=60° k=軸平行', 60, 'x'],
    ['Δφ=15° k=斜め', 15, 'diag'], ['Δφ=15° k=軸平行', 15, 'x'],
  ];

  console.log(`■ 条件（12×12接触子, ピッチ${P}mm, a=b=${AMP}mm, ${RPM}rpm, μ=${MU}, δp=${DP.v}mm）`);
  console.log(`  山での面内速度 a·q̇ = ${TIP.toFixed(0)} mm/s。速度の括弧内はこれに対する比。`);
  console.log(`  K は「静的釣合いの最大貫入の1周期平均 = δp」を二分法で解いた値。\n`);
  console.log('  条件               K[N/mm] 実測貫入  接地窓  位相ギャップ  静的平衡の傾き  最小接地  水平時の余裕');
  const T = {};
  for (const [lab, dphi, kdir] of CFG) {
    const pads = lattice(dphi, { kdir });
    const t = tune(pads);
    const sm = supportMargin(pads, BREF.v, t.K);
    T[lab] = { pads, t, sm, gap: phaseGap(pads) };
    console.log(`  ${lab.padEnd(16)} ${t.K.toFixed(1).padStart(6)} ${t.pen.toFixed(3)}mm  ±${t.window.toFixed(1)}°` +
      `   ${T[lab].gap.toFixed(0).padStart(3)}°        ${sm.maxTilt.toFixed(3).padStart(6)}°      ` +
      `${String(sm.minPts).padStart(3)} 点   ` +
      `${(sm.margin === -Infinity ? '  含まず' : sm.margin.toFixed(0).padStart(5) + 'mm')}` +
      `${sm.noSolution ? '  ⚠解なし' + sm.noSolution : ''}`);
  }
  console.log(`
  ※ 静的平衡は高さ・ピッチ・ロールを連立して解いている（${TSAMP} 位相で走査）。
    「最大位相ギャップ < 接地窓」は接触が切れない条件でしかなく、支持条件ではない。
    支持できるかは「解が存在するか」で、できるとしても機体は上表のぶん傾く。
    最後の列は水平姿勢に固定した場合の凸包余裕で、参考値（傾けば釣り合う配置でも「含まず」になる）。
`);

  console.log('■ 試験1: 推力方向が ψ に追従するか（駆動は A=B、方位保持なし）');
  const S1 = {};
  for (const [lab, dphi, kdir] of CFG) {
    const dirs = lab === 'Δφ=60° k=斜め' ? [0, 45, 90, 135, 180, 225, 270, 315] : [0, 45, 90, 135];
    let worst = 0;
    console.log(`  ${lab}`);
    for (const c of dirs) {
      const pads = lattice(dphi, { kdir });
      const t = T[lab].t;
      const r = core(pads, polarDrive({ shaftA: 1, shaftB: 1, hold: true, target: c * Math.PI / 180 }), { K: t.K, C: t.C });
      let e = r.dir - c; while (e > 180) e -= 360; while (e < -180) e += 360;
      worst = Math.max(worst, Math.abs(e));
      console.log(`    ψ*=${String(c).padStart(3)}°  ${fmt(r)}  誤差${e.toFixed(2).padStart(7)}°`);
    }
    S1[lab] = worst;
    console.log(`    → 最大操舵誤差 ${worst.toFixed(2)}°\n`);
  }

  console.log('  直交2ラティス（案A, Δφ=60° k=斜め, 停止軸リトラクトあり）');
  {
    const pads = lattice(60, { kdir: 'diag' });
    const t = tune(pads);
    let worst = 0;
    for (const c of [0, 45, 90, 135]) {
      const r = core(pads, shareDrive(c * Math.PI / 180), { K: t.K, C: t.C });
      let e = r.dir - c; while (e > 180) e -= 360; while (e < -180) e += 360;
      worst = Math.max(worst, Math.abs(e));
      console.log(`    指令=${String(c).padStart(3)}°  ${fmt(r)}  誤差${e.toFixed(2).padStart(7)}°`);
    }
    console.log(`    → 最大操舵誤差 ${worst.toFixed(2)}°\n`);
    S1['案A'] = worst;
  }

  console.log('■ 試験2: 空間位相をシャッフルすると推力が落ちるか');
  console.log('  ケースA 車両（大きな剛体の下に144個が全部入る、Δφ=60° k=斜め, ψ=45°）');
  for (const sh of [false, true]) {
    const pads = lattice(60, { kdir: 'diag', shuffle: sh });
    const t = tune(pads);
    const r = core(pads, polarDrive({ shaftA: 1, shaftB: 1, hold: true, target: Math.PI / 4 }), { K: t.K, C: t.C });
    console.log(`    ${(sh ? '同じ位相集合をシャッフル' : '空間波（勾配あり）').padEnd(24)} ${fmt(r)}`);
  }

  console.log('\n  ケースB 卓上（接触子は地面に固定、150mm角0.5kgの荷、b=1.5mm, k=斜め, ψ=45°）');
  console.log('    位相は φ(i,k)=(i+k)Δφ を添字から直接作る（周期境界なし）。荷も静的釣合いから始める。');
  console.log('    位相原点 q0 を 8 通り走査する。1点だけ見ると起動位相の当たり外れを');
  console.log('    「空間配置の効果」と取り違える。K は各行で空間波側の値に固定し、');
  console.log('    シャッフルは同じ位相集合 {m·Δφ} から一様に引く（位相配置だけの比較にする）。\n');
  console.log('    Δφ    波長(対角)   空間波: 完走/8  速度            シャッフル: 完走/8  速度');
  const Q0 = [0, 45, 90, 135, 180, 225, 270, 315];
  for (const dphi of [180, 120, 90, 60, 30, 15]) {
    const lam = (360 / dphi) * P * Math.SQRT1_2;
    const K0 = tableRun(dphi, { kdir: 'diag', q0: 0 }).K;
    const roll = shuffle => {
      let ok = 0; const sp = [];
      for (const q0 of Q0) {
        const r = tableRun(dphi, { kdir: 'diag', shuffle, q0: q0 * Math.PI / 180, Kfix: K0 });
        if (!r.fell) { ok += 1; sp.push(r.speed); }
      }
      return { ok, txt: sp.length ? `${Math.min(...sp).toFixed(0)}〜${Math.max(...sp).toFixed(0)} mm/s` : '—' };
    };
    const a = roll(false), b = roll(true);
    console.log(`    ${String(dphi).padStart(3)}°  ${lam.toFixed(0).padStart(5)} mm   ` +
      `${String(a.ok).padStart(8)}/8  ${a.txt.padEnd(16)}${String(b.ok).padStart(10)}/8  ${b.txt}`);
  }
  console.log('    → 荷の差し渡しは 150mm。完走 = 3.2 秒間 傾き 8° を超えなかった。');

  console.log('\n■ 試験3: 走行中に ψ を 0°→90° へ振る（Δφ=60° k=斜め）');
  console.log('    操舵レート    速度                方向     滑り   荷重');
  for (const rate of [0, 45, 180, 720]) {
    const R = rate * Math.PI / 180;
    const pads = lattice(60, { kdir: 'diag' });
    const t = T['Δφ=60° k=斜め'].t;
    const sweep = tt => { const raw = Math.max(0, tt - WARM) * R; return raw >= Math.PI / 2 ? [Math.PI / 2, 0] : [raw, R]; };
    const r = core(pads, polarDrive({ sweep }), { K: t.K, C: t.C });
    console.log(`    ${String(rate).padStart(4)}°/s   ${r.speed.toFixed(0).padStart(4)} mm/s (${(r.speed / TIP * 100).toFixed(0)}%)  ` +
      `${r.dir.toFixed(1).padStart(7)}°  ${pct(r.slip)}  ${r.load.toFixed(2)}±${pct(r.ripple)}`);
  }

  console.log('\n■ 試験4: 任意の固定方位を、実際のモータ入力で保持できるか（同じ6自由度コア）');
  console.log('    差動で ψ を目標へ寄せたあと、両軸を同速に戻して保持する制御。');
  console.log('    目標 ψ*   保持した ψ   推力方向    ψ*との差   速度');
  for (const c of [0, 45, 90, 135, 180, 225, 270, 315]) {
    const pads = lattice(60, { kdir: 'diag' });
    const t = T['Δφ=60° k=斜め'].t;
    const r = core(pads, polarDrive({ shaftA: 1, shaftB: 1, hold: true, target: c * Math.PI / 180 }), { K: t.K, C: t.C });
    let e = r.dir - c; while (e > 180) e -= 360; while (e < -180) e += 360;
    let dp2 = r.psi - c; while (dp2 > 180) dp2 -= 360; while (dp2 < -180) dp2 += 360;
    console.log(`    ${String(c).padStart(4)}°   ${dp2.toFixed(2).padStart(7)}°ずれ  ${r.dir.toFixed(1).padStart(7)}°   ` +
      `${e.toFixed(2).padStart(6)}°   ${r.speed.toFixed(0)} mm/s`);
  }

  console.log('\n■ 上下振幅 b を下げると跳ねは減るか（Δφ=60° k=斜め, a=6mm 固定, ψ=0）');
  console.log('    b と一緒に K・C を動かすと効果が分離できないので、両方出す。');
  console.log('    b       K,C を b ごとに整定                       K,C を b=6 の値に固定');
  const base = T['Δφ=60° k=斜め'].t;
  for (const b of [6, 3, 1.5, 0.8]) {
    BREF.v = b;
    const pads = lattice(60, { kdir: 'diag' });
    const ta = tune(pads, b);
    const ra = core(pads, polarDrive({ shaftA: 1, shaftB: 1 }), { b, K: ta.K, C: ta.C });
    const rb = core(lattice(60, { kdir: 'diag' }), polarDrive({ shaftA: 1, shaftB: 1 }), { b, K: base.K, C: base.C });
    const g = r => `${r.speed.toFixed(0).padStart(4)}mm/s 滑り${pct(r.slip)} 荷重${r.load.toFixed(2)}±${pct(r.ripple)}`;
    console.log(`    ${(b + 'mm').padEnd(6)}  K=${ta.K.toFixed(1).padStart(5)} ${g(ra)}   K=${base.K.toFixed(1)} ${g(rb)}`);
  }
  BREF.v = 6;

  console.log('\n■ 組立公差に対する頑健性（Δφ=60° k=斜め, ψ=0, δp=0.3mm）');
  console.log('    足ごとに固定の高さ誤差を与える（正規分布 σ）。接地窓は δp で決まるので、');
  console.log('    公差が δp に対して無視できなければ「一部の足だけが常に荷重を持つ」はず、');
  console.log('    という予想を確かめる。\n');
  console.log('    高さ誤差 σ   速度              滑り    荷重リップル  平均接地  最小接地');
  for (const bb of [6, 1.5]) {
    BREF.v = bb;
    console.log(`    --- 上下振幅 b = ${bb} mm（ストローク ${2 * bb} mm）`);
    for (const tol of [0, 0.05, 0.1, 0.2, 0.5]) {
      const pads = lattice(60, { kdir: 'diag', tol });
      const tt = tune(pads, bb);
      const r = core(pads, polarDrive({ shaftA: 1, shaftB: 1 }), { b: bb, K: tt.K, C: tt.C });
      console.log(`    σ=${String(tol).padStart(4)}mm    ${r.pathSpeed.toFixed(0).padStart(4)} mm/s ` +
        `(${(r.pathSpeed / TIP * 100).toFixed(0).padStart(3)}%)  ${(r.slip * 100).toFixed(1).padStart(5)}%  ` +
        `${(r.ripple * 100).toFixed(0).padStart(6)}%   ${r.contact.toFixed(1).padStart(5)}    ${String(r.minC).padStart(3)}`);
    }
  }
  BREF.v = 6;
  console.log(`
    → **予想は外れた。** δp=0.3mm に対して σ=0.5mm、つまり公差が接地窓より大きくても
      速度は 94〜97% 出る。b=6mm では跳ねがむしろ収まる（リップル 120% → 70%）。
      高さのばらつきが同期した離床を崩すため。組立公差はこの機構の律速ではない。
`);

  console.log(`
■ この試験が言えること / 言えないこと

  言えること
   ・操舵は構造的に厳密。接地している接触子の摩擦力が全部 ψ に共線なので、
     ψ に直交する速度成分は必ず減衰する。案Aの過拘束（X が (r·ωx,0)、Y が (0,r·ωy) を
     要求して競合する）が存在せず、停止軸もないのでリトラクト機構も要らない。
   ・差動で追い込んでから両軸を揃える制御で、任意の固定方位を保持できる（試験4）。
   ・k を軸に平行にすると同位相の接触子が1列に並んで線接触になり、機体が首を振る。
   ・車両（荷 ≫ λ）では位相をシャッフルしても推力が落ちない。ここでは波ではなく当番表。

  言えないこと
   ・**部品数が重い。** 足あたり最小3部品（偏心・上下カム・回る案内溝）で、
     12×12 なら 432 部品。逆旋回2偏心で作るなら 1296（bench/mechanism.mjs）。
   ・**厳密には非ホロノミックで、2モータではヨーが出ない**（bench/mechanism.mjs）。
     ψ を空間分布させればヨーは作れるが、滑りが 0.6% → 15.1% になる。
   ・**跳ねている。** b=6mm では最小接地数が 0、荷重リップルが 100% を超える。
     上の速度はどれも「周期的に完全離床する機体」の値。b を下げれば収まるが、
     K・C も一緒に動くので効果の分離は上の表を見ること。
   ・**卓上（有限の荷）の設計則は立っていない。** 周期境界を外して静的釣合いから
     始めると Δφ 依存が単調でなくなる。以前ここに書いていた「λ ≲ 荷 が要る」は撤回。
     言えるのは「位相をシャッフルすると必ず転倒する」ことまで。
   ・μ・接触剛性・接地窓はすべて未測定の仮定値。
   ・差動歯車のバックラッシュ、位相保持トルク、連続スキンの面内せん断、接触子間干渉は
     どのモデルにも入っていない。
   ・2入力では並進2自由度まで。ヨーを含む平面3自由度には一般化力が3系統要る。`);
}
