// 全方式を1つの物理コアで比較する。
//
// 未測定のパラメータ（摩擦係数 μ、パッドの沈み込み δp）が結論を左右するので、
// 1点で比べずに振って、どの方式が頑健に勝つかを見る。
//
// 物理コア: 各足の対地滑りに逆らうクーロン摩擦を法線荷重で重み付けして足し、
//          車体（または荷）の速度が釣り合うところまで積分する。

export const BASE = {
  rails: 24, dphiDeg: 15, crank: 6, rpm: 300,
  discR: 30, alphaDeg: 15,
  mu: 0.6, comp: 0.30,          // comp = パッドの沈み込み δp [mm]
  veps: 6, dt: 2e-4, T: 12, settle: 4
};

const G = 9810;

/* 方式ごとに「足の高さ」と「足の対車体速度」を与える。あとは共通。 */
const ARCH = {
  // 直交2ラティス・荷重分担（案A そのもの）
  shareXY: {
    label: '直交2ラティス 荷重分担', actuators: 2, rails: 2,
    lateral: true, yaw: false, decoupled: false,
    feet(p, t, cmd) {
      const out = [];
      const wx = p.w*Math.cos(cmd), wy = p.w*Math.sin(cmd);
      for (let j = 0; j < p.rails; j += 1) {
        const a = j*p.dphi - wx*t, b = j*p.dphi - wy*t;
        out.push({ h: -p.crank*Math.cos(a), t: [-p.crank*wx*Math.cos(a), 0] });
        out.push({ h: -p.crank*Math.cos(b), t: [0, -p.crank*wy*Math.cos(b)] });
      }
      return out;
    }
  },

  // 直交2ラティス・時間交替（片方ずつ接地させる）
  altXY: {
    label: '直交2ラティス 時間交替', actuators: 3, rails: 2,
    lateral: true, yaw: false, decoupled: false,
    feet(p, t, cmd) {
      const out = [];
      const wx = p.w*Math.cos(cmd), wy = p.w*Math.sin(cmd);
      const s = Math.max(-1, Math.min(1, Math.cos(2*Math.PI*p.altHz*t)*2.5))*0.5 + 0.5;
      const liftX = 2*p.crank*(1-s), liftY = 2*p.crank*s;
      for (let j = 0; j < p.rails; j += 1) {
        const a = j*p.dphi - wx*t, b = j*p.dphi - wy*t;
        out.push({ h: -p.crank*Math.cos(a) + liftX, t: [-p.crank*wx*Math.cos(a), 0] });
        out.push({ h: -p.crank*Math.cos(b) + liftY, t: [0, -p.crank*wy*Math.cos(b)] });
      }
      return out;
    }
  },

  // 1ラティスをまるごと指令方向へ旋回（オフセット無し）
  steer1: {
    label: '1ラティス旋回', actuators: 2, rails: 1,
    lateral: 'slew', yaw: false, decoupled: false,
    feet(p, t, cmd) {
      const out = [];
      const c = Math.cos(cmd), s = Math.sin(cmd);
      for (let j = 0; j < p.rails; j += 1) {
        const a = j*p.dphi - p.w*t;
        const m = -p.crank*p.w*Math.cos(a);
        out.push({ h: -p.crank*Math.cos(a), t: [m*c, m*s] });
      }
      return out;
    }
  },

  // 先端を傾斜円盤にして、ピン自転をクランクに 1:1 同期、位相差 δ で操舵
  disc: {
    label: '傾斜円盤＋位相操舵', actuators: 2, rails: 1,
    lateral: true, yaw: false, decoupled: true,
    feet(p, t, cmd) {
      const out = [];
      const e = p.discR*Math.cos(p.alpha);
      for (let j = 0; j < p.rails; j += 1) {
        const a = j*p.dphi - p.w*t, chi = a + cmd;
        out.push({ h: -p.crank*Math.cos(a),
                   t: [p.w*e*Math.sin(chi), -p.w*e*Math.cos(chi)] });
      }
      return out;
    }
  }
};

function prep(o) {
  const p = { ...BASE, ...o };
  p.dphi = p.dphiDeg*Math.PI/180;
  p.alpha = p.alphaDeg*Math.PI/180;
  p.w = p.rpm*Math.PI/30;
  p.altHz = p.altHz ?? 3;
  p.tip = p.arch === 'disc' ? p.discR*Math.cos(p.alpha)*p.w : p.crank*p.w;
  return p;
}

export function simulate(arch, opts = {}, cmd = Math.PI/4) {
  const p = prep({ ...opts, arch });
  const A = ARCH[arch];
  const muG = p.mu*G;
  let vx = 0, vz = 0, dx = 0, dz = 0, slip = 0, wt = 0;

  for (let s = 0; s < p.T/p.dt; s += 1) {
    const t = s*p.dt;
    const F = A.feet(p, t, cmd);

    let lo = Infinity;
    for (const f of F) if (f.h < lo) lo = f.h;
    const plane = lo + p.comp;

    let tot = 0;
    for (const f of F) { f.n = Math.max(0, plane - f.h); tot += f.n; }
    if (tot <= 0) continue;

    let fx = 0, fz = 0, sl = 0;
    for (const f of F) {
      if (f.n <= 0) continue;
      const n = f.n/tot;
      const gx = vx + f.t[0], gz = vz + f.t[1];
      const m = Math.max(Math.hypot(gx, gz), p.veps);
      fx -= n*gx/m; fz -= n*gz/m;
      sl += n*Math.hypot(gx, gz);
    }
    vx += muG*fx*p.dt; vz += muG*fz*p.dt;

    if (t > p.settle) { dx += vx*p.dt; dz += vz*p.dt; slip += sl*p.dt; wt += p.dt; }
  }

  const span = p.T - p.settle;
  return { speed: Math.hypot(dx, dz)/span, dir: Math.atan2(dz, dx)*180/Math.PI,
           slip: slip/wt/p.tip, tip: p.tip, ratio: Math.hypot(dx, dz)/span/p.tip };
}

/* 操舵の直線性: 指令を振って、最良の線形当てはめからのずれの最大値 */
export function steering(arch, opts = {}) {
  const cmds = [];
  for (let d = 5; d <= 85; d += 10) cmds.push(d);
  const rows = cmds.map(d => ({ cmd: d, got: simulate(arch, opts, d*Math.PI/180).dir }));
  let best = Infinity, off = 0;
  for (let o = -180; o <= 180; o += 0.5) {
    let m = 0;
    for (const r of rows) { let e = r.got - r.cmd - o; while (e>180) e-=360; while (e<-180) e+=360; m = Math.max(m, Math.abs(e)); }
    if (m < best) { best = m; off = o; }
  }
  return { maxErr: best, offset: off, rows };
}

export { ARCH, prep };
