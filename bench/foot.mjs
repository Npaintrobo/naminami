// 足1本の設計。傾斜円盤を単体機構として見たときの寸法決め。
//
// 円盤を軸に対して α 傾けて剛結し、その軸を鉛直から β 傾けて自転させる。
// 円盤の法線はピン軸まわりに歳差するので、鉛直との角 γ は |α-β| 〜 α+β を往復する。
// 最下点の深さは R·sinγ なので、1回転で「上下（接地・離地）」と「掃引（推力）」が
// 両方出る。上下用のクランクも、自転用のギアも、位相シフタも要らない。
//
//   ストローク = 2·R·cosα·sinβ
//   張り出し   e = R·cos(α+β)      接地している瞬間の、ピン軸から接地点までの水平距離
//   掃引速度   v ≈ 0.9·R·ω
//   制約       β < α               満たさないと接地点が一周せず、全方位へ向けられない
//              e < ピッチ/2         隣の足とぶつからない
//
// 操舵はハウジングを鉛直まわりに φ 回すだけ。全ユニットで同じ φ なら接地点の
// オフセットも同一方向に動くので、パッチ全体を回すことにならず、スクラブが出ない。

const deg = d => d * Math.PI / 180;
const dot = (a, b) => a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const scal = (a, k) => [a[0]*k, a[1]*k, a[2]*k];
const sub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const norm = a => { const m = Math.hypot(a[0], a[1], a[2]); return [a[0]/m, a[1]/m, a[2]/m]; };

/* Rodrigues: v を単位軸 k まわりに th 回す */
function rot(v, k, th) {
  const c = Math.cos(th), s = Math.sin(th);
  const kv = cross(k, v), kd = dot(k, v);
  return [v[0]*c + kv[0]*s + k[0]*kd*(1-c),
          v[1]*c + kv[1]*s + k[1]*kd*(1-c),
          v[2]*c + kv[2]*s + k[2]*kd*(1-c)];
}

const UP = [0, 1, 0];

/* 1回転ぶんの接地点の軌跡を出し、接地窓の中の平均と、ばらつき（＝滑り）を返す。
   squash はパッドの沈み込み [mm]。深さが最大から squash 以内なら接地とみなす。 */
export function foot({ R = 30, alphaDeg = 30, betaDeg = 25, phiDeg = 0, rpm = 300, squash = 0.3, steps = 3000 } = {}) {
  const a = deg(alphaDeg), b = deg(betaDeg), f = deg(phiDeg), w = rpm * Math.PI / 30;

  const p  = [Math.sin(b)*Math.cos(f), Math.cos(b), Math.sin(b)*Math.sin(f)];       // ピン軸
  const n0 = [Math.sin(b+a)*Math.cos(f), Math.cos(b+a), Math.sin(b+a)*Math.sin(f)]; // 法線の初期姿勢

  const rows = [];
  for (let i = 0; i < steps; i += 1) {
    const n = rot(n0, p, i / steps * 2 * Math.PI);
    const g = Math.acos(Math.max(-1, Math.min(1, dot(n, UP))));      // 鉛直との角
    const down = scal(norm(sub(UP, scal(n, dot(UP, n)))), -1);        // 円盤面内で最も下る向き
    const L = scal(down, R);                                         // 最下点（円盤中心が原点）
    rows.push({ depth: R * Math.sin(g), L, v: scal(cross(p, L), w) });
  }

  const dmax = Math.max(...rows.map(r => r.depth));
  const dmin = Math.min(...rows.map(r => r.depth));
  const win = rows.filter(r => r.depth > dmax - squash);

  let sx = 0, sz = 0;
  for (const r of win) { sx += r.v[0]; sz += r.v[2]; }
  sx /= win.length; sz /= win.length;

  let sl = 0;
  for (const r of win) sl += Math.hypot(r.v[0] - sx, r.v[2] - sz);
  sl /= win.length;

  const mean = Math.hypot(sx, sz);
  return {
    stroke: dmax - dmin,
    duty: win.length / steps,
    sweep: mean,
    slip: sl / Math.max(mean, 1e-9),
    reach: Math.hypot(win[0].L[0], win[0].L[2]),
    dir: Math.atan2(sz, sx) * 180 / Math.PI
  };
}

if (import.meta.url === `file://${process.argv[1]}`) {
  console.log('■ α と β の選び方（R=30mm, 300rpm, 沈み込み 0.3mm）');
  console.log('  α    β    ストローク  接地率  掃引速度   滑り   張り出し e');
  for (const [al, be] of [[15,5],[15,10],[20,10],[20,15],[30,20],[30,25]]) {
    const r = foot({ alphaDeg: al, betaDeg: be });
    console.log(`  ${String(al).padStart(2)}°  ${String(be).padStart(2)}°   ${r.stroke.toFixed(1).padStart(5)} mm  ${(r.duty*100).toFixed(0).padStart(4)}%  ${r.sweep.toFixed(0).padStart(5)} mm/s  ${(r.slip*100).toFixed(0).padStart(3)}%   ${r.reach.toFixed(0)} mm`);
  }
  console.log('  → β を α に近づけるほどストロークが増えて滑りが減る（深さプロファイルが尖る）。');
  console.log('    かわりに張り出し e が縮むので掃引速度は落ちる。推奨は α=30°, β=25°。\n');

  console.log('■ 操舵: ハウジングを鉛直まわりに φ 回す（α=30°, β=25°）');
  console.log('  φ       推力の向き   φ からの差   掃引速度');
  const ref = foot({ phiDeg: 0 }).dir;
  for (const f of [0, 45, 90, 180, 270]) {
    const r = foot({ phiDeg: f });
    let d = r.dir - ref - f;
    while (d > 180) d -= 360;
    while (d < -180) d += 360;
    console.log(`  ${String(f).padStart(3)}°   ${r.dir.toFixed(1).padStart(7)}°   ${d.toFixed(2).padStart(6)}°   ${r.sweep.toFixed(0)} mm/s`);
  }
  console.log('  → 1:1 で追従。全ユニットで φ を揃えれば接地点は同一方向に動くのでスクラブが出ない。\n');

  console.log('■ 接地率と滑りのトレードオフ（α=30°, β=25°）');
  console.log('  沈み込み  接地率   滑り   連続接地に要る足の数');
  for (const sq of [0.1, 0.3, 0.6, 1.2, 2.4]) {
    const r = foot({ squash: sq });
    console.log(`  ${sq.toFixed(1)} mm    ${(r.duty*100).toFixed(0).padStart(3)}%   ${(r.slip*100).toFixed(0).padStart(3)}%        ${Math.ceil(1/r.duty)} 個以上`);
  }
  console.log('  → 足1本では接地率が低いので単体では成立しない。位相をずらした配列が要る。\n');

  console.log('■ 推奨仕様（R=30mm, α=30°, β=25°, 300rpm, 沈み込み 0.3mm）');
  const r = foot();
  console.log(`  ストローク ${r.stroke.toFixed(0)} mm / 掃引速度 ${r.sweep.toFixed(0)} mm/s / 滑り ${(r.slip*100).toFixed(0)}%`);
  console.log(`  張り出し ${r.reach.toFixed(0)} mm → 足のピッチは ${(2*r.reach).toFixed(0)} mm 以上、実用上は 40 mm`);
  console.log(`  接地率 ${(r.duty*100).toFixed(0)}% → 安定支持には 25〜36 個`);
}
