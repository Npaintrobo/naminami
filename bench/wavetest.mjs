// 「これはまだ波動機構か」を判定する。
//
// 進行波駆動なら、位相が空間的に勾配していること（＝波が伝わること）が本質のはず。
// 位相をランダムに入れ替えても同じ性能が出るなら、位相差は「常にどれかが接地している」
// ための当番表でしかなく、波ではない。

const N = 24, RPM = 300, CRANK = 6, DISC_E = 25;
const MU_G = 0.6 * 9810, VEPS = 6, COMP = 0.3;
const W = RPM * Math.PI / 30;

function run(phases, delta = Math.PI/4, T = 12) {
  let vx = 0, vz = 0, ph = 0, dx = 0, dz = 0, slip = 0, wt = 0, minContact = 1e9;
  const dt = 2e-4;

  for (let s = 0; s < T/dt; s += 1) {
    ph += W * dt;
    const h = [], t = [];
    for (let j = 0; j < N; j += 1) {
      const th = phases[j] - ph, chi = th + delta;
      h.push(-CRANK * Math.cos(th));
      t.push([W*DISC_E*Math.sin(chi), -W*DISC_E*Math.cos(chi)]);
    }
    const plane = Math.min(...h) + COMP;
    let tot = 0, cnt = 0;
    const load = h.map(v => { const b = Math.max(0, plane - v); if (b > 0) cnt += 1; tot += b; return b; });
    if (s*dt > 2 && cnt < minContact) minContact = cnt;

    let fx = 0, fz = 0, sl = 0;
    for (let j = 0; j < N; j += 1) {
      if (load[j] <= 0) continue;
      const n = load[j]/tot;
      const gx = vx + t[j][0], gz = vz + t[j][1];
      const m = Math.max(Math.hypot(gx, gz), VEPS);
      fx -= n*gx/m; fz -= n*gz/m; sl += n*Math.hypot(gx, gz);
    }
    vx += MU_G*fx*dt; vz += MU_G*fz*dt;
    if (s*dt > 2) { dx += vx*dt; dz += vz*dt; slip += sl*dt; wt += dt; }
  }
  const tip = DISC_E*W;
  return { speed: Math.hypot(dx,dz)/(T-2), dir: Math.atan2(dz,dx)*180/Math.PI,
           slip: slip/wt/tip, minContact };
}

const TAU = 2*Math.PI;
const linear  = Array.from({length: N}, (_, i) => i * TAU/N);              // 進行波（1波長）
const linear2 = Array.from({length: N}, (_, i) => i * 2*TAU/N);            // 進行波（2波長）
const groups3 = Array.from({length: N}, (_, i) => (i % 3) * TAU/3);        // 3相グループ
const groups2 = Array.from({length: N}, (_, i) => (i % 2) * TAU/2);        // 2相グループ
const same    = Array.from({length: N}, () => 0);                          // 全部同位相

// 乱数（再現性のため固定シード）
let seed = 12345;
const rnd = () => (seed = (seed*1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
const random = Array.from({length: N}, () => rnd() * TAU);
// 空間的な勾配を完全に壊す: 線形位相をシャッフル
const shuffled = linear.slice();
for (let i = N-1; i > 0; i -= 1) { const j = Math.floor(rnd()*(i+1)); [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]]; }

console.log('■ 位相の配り方を変える（傾斜円盤 24足, δ=45°）');
console.log('  配り方                      速度        方向      滑り   最小接地数');
for (const [label, ph] of [
  ['進行波・1波長（勾配あり）', linear],
  ['進行波・2波長（勾配あり）', linear2],
  ['線形位相をシャッフル',      shuffled],
  ['完全ランダム',              random],
  ['3相グループ',               groups3],
  ['2相グループ',               groups2],
  ['全部同位相',                same],
]) {
  const r = run(ph);
  console.log(`  ${label.padEnd(24)} ${r.speed.toFixed(0).padStart(4)} mm/s  ${r.dir.toFixed(1).padStart(6)}°  ${(r.slip*100).toFixed(0).padStart(3)}%   ${String(r.minContact).padStart(3)}`);
}
