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
function lattice(dphiDeg, { shuffle = false, seed = 7, kdir = 'diag' } = {}) {
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
  return pads;
}

/* ── 接地剛性の整定 ─────────────────────────────
   「δp」は静的に置いたときの最大貫入と定義する。位相によって振れるので1周期で平均し、
   その平均が δp になる K を二分法で解く。以前は総貫入量から一発で決めていて、
   実際に釣り合わせた貫入（0.12mm）と設定値（0.3mm）が合っていなかった。 */
function equilibrium(ys, K, weight) {
  let lo = -Math.max(...ys.map(Math.abs)) - 30, hi = -Math.min(...ys) + 1;
  for (let i = 0; i < 60; i += 1) {
    const mid = (lo + hi) / 2;
    let f = 0;
    for (const y of ys) f += Math.max(0, -(mid + y)) * K;
    if (f > weight) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

function tune(pads, b = BREF.v, dp = DP.v, mass = MASS) {
  const weight = mass * G / ACC;
  const S = 24;
  const phases = [];
  for (let s = 0; s < S; s += 1) {
    const q = s / S * TAU;
    phases.push(pads.map(p => -b * Math.cos(p.ph - q)));
  }
  const meanMax = K => {
    let acc = 0;
    for (const ys of phases) {
      const Y = equilibrium(ys, K, weight);
      let mx = 0;
      for (const y of ys) mx = Math.max(mx, -(Y + y));
      acc += mx;
    }
    return acc / S;
  };
  let klo = 0.05, khi = 20000;
  for (let i = 0; i < 50; i += 1) {
    const km = (klo + khi) / 2;
    if (meanMax(km) > dp) klo = km; else khi = km;      // K が大きいほど貫入は小さい
  }
  const K = (klo + khi) / 2;
  const actual = meanMax(K);

  let ncSum = 0;
  for (const ys of phases) {
    const Y = equilibrium(ys, K, weight);
    for (const y of ys) if (-(Y + y) > 0) ncSum += 1;
  }
  const nc = Math.max(1, ncSum / S);
  const C = 2 * 0.35 * Math.sqrt(nc * K * 1000 * mass) / 1000 / nc;
  return { K, C, nc, pen: actual, window: Math.acos(clamp(1 - actual / b, -1, 1)) * 180 / Math.PI };
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

function supportMargin(pads, b, K, dp) {
  const weight = MASS * G / ACC;
  const S = 72;
  let worstMargin = Infinity, worstCnt = 1e9;
  for (let s = 0; s < S; s += 1) {
    const q = s / S * TAU;
    const ys = pads.map(p => -b * Math.cos(p.ph - q));
    const Y = equilibrium(ys, K, weight);
    const pts = [];
    for (let i = 0; i < pads.length; i += 1) if (-(Y + ys[i]) > 0) pts.push([pads[i].x, pads[i].z]);
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
  void dp;
  return { margin: worstMargin, minPts: worstCnt };
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
  /* 静的釣合いは運動学ごとに高さの式が違うので、drv.foot から取る（偏極式を決め打ちしない） */
  const ys = pads.map(p => drv.foot(p)[0]);
  let Y = equilibrium(ys, K, weight);
  for (let i = 0; i < pads.length; i += 1) pads[i].prev = Math.max(0, -(Y + ys[i]));

  let X = 0, Z = 0, yaw = 0, tx = 0, tz = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  let dx = 0, dz = 0, slip = 0, wt = 0, nmean = 0, nvar = 0, minC = 1e9, cntSum = 0, tilt2 = 0;
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
      slip += (wsum > 0 ? sl / wsum : 0) * DT;
      const ld = SN / weight;
      nmean += ld * DT; nvar += ld * ld * DT;
      cntSum += cnt * DT; tilt2 += (tx * tx + tz * tz) * DT;
      if (cnt < minC) minC = cnt;
    }
  }

  const mn = nmean / wt;
  return {
    speed: Math.hypot(dx, dz) / wt, dir: Math.atan2(dz, dx) * 180 / Math.PI,
    slip: slip / wt / TIP, load: mn,
    ripple: Math.sqrt(Math.max(0, nvar / wt - mn * mn)) / Math.max(mn, 1e-9),
    contact: cntSum / wt, minC, tilt: Math.sqrt(tilt2 / wt) * 180 / Math.PI,
    psi: drv.psi * 180 / Math.PI
  };
}

/* ── 駆動 ───────────────────────────────────── */
function polarDrive({ shaftA = 1, shaftB = 1, hold = false, target = 0, sweep = null } = {}) {
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
      const th = p.ph - q, c = Math.cos(psi), s = Math.sin(psi);
      const sn = Math.sin(th), cs = Math.cos(th);
      return [-BREF.v * cs, AMP * sn * c, AMP * sn * s,
              -AMP * qd * cs * c - AMP * psid * sn * s,
              -AMP * qd * cs * s + AMP * psid * sn * c];
    },
    q: 0, psi: 0
  };
  return d;
}

function shareDrive(cmd, retract = true) {
  let t0 = 0;
  const d = {
    reset() { t0 = 0; d.psi = cmd; },
    step(t) { t0 = t; },
    foot(p) {
      const g = rampOf(t0);
      const wx = W * g * Math.cos(cmd), wz = W * g * Math.sin(cmd);
      /* 指令が 0 の軸だけを退避させる。以前は速度どうしを比較していて、
         起動時（両方 0）に両ラティスが持ち上がり発散していた。 */
      const liftX = retract && Math.abs(Math.cos(cmd)) < 1e-9 ? 2 * BREF.v : 0;
      const liftZ = retract && Math.abs(Math.sin(cmd)) < 1e-9 ? 2 * BREF.v : 0;
      if (p.xy === 0) { const th = p.ph - wx * t0; return [-BREF.v * Math.cos(th) + liftX, 0, 0, -AMP * wx * Math.cos(th), 0]; }
      const th = p.ph - wz * t0; return [-BREF.v * Math.cos(th) + liftZ, 0, 0, 0, -AMP * wz * Math.cos(th)];
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
function tableRun(dphiDeg, { kdir = 'diag', shuffle = false, psi = Math.PI / 4,
                             b = 1.5, S = 150, M = 0.5, HCG2 = 15, dp = 0.4 } = {}) {
  const d = dphiDeg * Math.PI / 180;
  const hash = (i, k) => {
    let h = (i * 73856093) ^ (k * 19349663);
    h = (h ^ (h >>> 13)) * 1274126177;
    return ((h >>> 0) / 4294967296) * TAU;
  };
  const at = (i, k) => (shuffle ? hash(i, k) : (kdir === 'diag' ? i + k : i) * d);

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

  for (let s = 0; s < T / DT; s += 1) {
    const t = s * DT, live = t > WARM, g = rampOf(t);
    if (Math.hypot(tx, tz) > TIPPED) return { fell: t, K };
    const cp = Math.cos(yaw), spz = Math.sin(yaw);
    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0;

    for (const [i, k] of under(X, Z)) {
      const rx = i * P - X, rz = k * P - Z;
      const bx = rx * cp + rz * spz, bz = -rx * spz + rz * cp;
      if (Math.abs(bx) > S / 2 || Math.abs(bz) > S / 2) continue;
      const th = at(i, k) - W * g * t;
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
  const pct = v => (v * 100).toFixed(0).padStart(3) + '%';
  const fmt = r => `${r.speed.toFixed(0).padStart(4)} mm/s (${(r.speed / TIP * 100).toFixed(0).padStart(3)}%)  ` +
    `${r.dir.toFixed(1).padStart(7)}°  滑り${pct(r.slip)}  荷重${r.load.toFixed(2)}±${pct(r.ripple)}  ` +
    `傾き${r.tilt.toFixed(2).padStart(5)}°  接地${r.contact.toFixed(1).padStart(4)}(最小${r.minC})`;

  const CFG = [
    ['Δφ=60° k=斜め', 60, 'diag'], ['Δφ=60° k=軸平行', 60, 'x'],
    ['Δφ=15° k=斜め', 15, 'diag'], ['Δφ=15° k=軸平行', 15, 'x'],
  ];

  console.log(`■ 条件（12×12接触子, ピッチ${P}mm, a=b=${AMP}mm, ${RPM}rpm, μ=${MU}, δp=${DP.v}mm）`);
  console.log(`  山での面内速度 a·q̇ = ${TIP.toFixed(0)} mm/s。速度の括弧内はこれに対する比。`);
  console.log(`  K は「静的釣合いの最大貫入の1周期平均 = δp」を二分法で解いた値。\n`);
  console.log('  条件                K[N/mm]  実測貫入  接地窓   最大位相ギャップ  静的支持余裕  最小接地点数');
  const T = {};
  for (const [lab, dphi, kdir] of CFG) {
    const pads = lattice(dphi, { kdir });
    const t = tune(pads);
    const sm = supportMargin(pads, BREF.v, t.K, DP.v);
    T[lab] = { pads, t, sm, gap: phaseGap(pads) };
    console.log(`  ${lab.padEnd(16)} ${t.K.toFixed(1).padStart(6)}  ${t.pen.toFixed(3)}mm  ±${t.window.toFixed(1)}°` +
      `      ${T[lab].gap.toFixed(0).padStart(3)}°        ` +
      `${(sm.margin === -Infinity ? '  支持不能' : sm.margin.toFixed(0).padStart(6) + 'mm')}      ${String(sm.minPts).padStart(3)} 点`);
  }
  console.log('\n  ※ 「最大位相ギャップ < 接地窓の全幅」は接触が切れない条件でしかない。');
  console.log('    姿勢を保持できるかは、接地点の凸包が重心投影を含むか（静的支持余裕 > 0）で決まる。\n');

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
  console.log('    Δφ    波長(対角)   空間波                          シャッフル');
  for (const dphi of [180, 120, 90, 60, 30, 15]) {
    const lam = (360 / dphi) * P * Math.SQRT1_2;
    const a = tableRun(dphi, { kdir: 'diag' }), b = tableRun(dphi, { shuffle: true });
    const f = r => (r.fell ? `転倒 t=${r.fell.toFixed(2)}s`
      : `${r.speed.toFixed(0)}mm/s ${r.dir.toFixed(0)}° 無接地${r.drop.toFixed(0)}% 傾き${r.tilt.toFixed(1)}°`).padEnd(30);
    console.log(`    ${String(dphi).padStart(3)}°  ${lam.toFixed(0).padStart(5)} mm   ${f(a)}  ${f(b)}`);
  }
  console.log('    → 荷の差し渡しは 150mm。');

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
