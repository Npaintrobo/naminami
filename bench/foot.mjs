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
//
// ピッチの下限は「接地している瞬間の張り出し」では決まりません。離地中も円盤は回り続けていて、
// 干渉判定に要るのは1回転ぶんの掃引体積の外形（envelope）＝縁の全点・全回転角での最大水平半径です。
//
// これは解析的に必ず R になります。どんな傾きでも、円盤面と水平面の交線（節線）上にある
// 縁の2点は水平面内にあり、中心からの水平距離はちょうど R だからです。α, β には依りません。
// したがって安全側のピッチ下限は 2R。これより詰めるなら、隣接ユニットの位相差に依存する
// 対ごとの掃引体積干渉と、軸・ハウジング・支持部材まで含めた検査が要ります（本ファイル未対応）。
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
    rows.push({ depth: R * Math.sin(g), L, v: scal(cross(p, L), w), reach: Math.hypot(L[0], L[2]) });
  }

  /* 掃引体積の外形: 縁の全点を全回転角で見たときの、ピン軸まわりの最大水平半径。 */
  let envelope = 0;
  for (let i = 0; i < steps; i += 1) {
    const n = rot(n0, p, i / steps * 2 * Math.PI);
    const e1 = norm(cross(n, Math.abs(n[1]) > 0.9 ? [1, 0, 0] : UP));
    const e2 = cross(n, e1);
    for (let u = 0; u < 180; u += 1) {
      const a2 = u / 180 * 2 * Math.PI, c = Math.cos(a2) * R, s2 = Math.sin(a2) * R;
      const x = e1[0] * c + e2[0] * s2, z = e1[2] * c + e2[2] * s2;
      const h = Math.hypot(x, z);
      if (h > envelope) envelope = h;
    }
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
    reach: Math.max(...win.map(r => r.reach)),        // 接地窓の中での最大張り出し
    reachAll: Math.max(...rows.map(r => r.reach)),    // 1回転での接地点の最大張り出し
    envelope,                                          // 円盤そのものの掃引外形半径
    dir: Math.atan2(sz, sx) * 180 / Math.PI
  };
}

if (import.meta.url === `file://${process.argv[1]}`) {
  console.log('■ α と β の選び方（R=30mm, 300rpm, 沈み込み 0.3mm）');
  console.log('  α    β    ストローク  接地率  掃引速度   滑り   接地窓e   掃引外形');
  for (const [al, be] of [[15,5],[15,10],[20,10],[20,15],[30,20],[30,25]]) {
    const r = foot({ alphaDeg: al, betaDeg: be });
    console.log(`  ${String(al).padStart(2)}°  ${String(be).padStart(2)}°   ${r.stroke.toFixed(1).padStart(5)} mm  ${(r.duty*100).toFixed(0).padStart(4)}%  ${r.sweep.toFixed(0).padStart(5)} mm/s  ${(r.slip*100).toFixed(0).padStart(3)}%  ${r.reach.toFixed(1).padStart(5)} mm  ${r.envelope.toFixed(1).padStart(5)} mm`);
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
  console.log(`  接地窓での張り出し ${r.reach.toFixed(1)} mm / 1回転での接地点の張り出し ${r.reachAll.toFixed(1)} mm`);
  console.log(`  掃引体積の外形 ${r.envelope.toFixed(1)} mm（= R。傾きに依らない）`);
  console.log(`  → 足のピッチは 2R = ${(2*r.envelope).toFixed(0)} mm 以上。`);
  console.log(`    以前ここに書いていた「実用上 40 mm」は誤りでした。接地窓の張り出し ${r.reach.toFixed(1)} mm`);
  console.log(`    だけを見ていて、離地中に円盤が水平へ寝るぶんを数えていませんでした。`);
  console.log(`    40 mm に詰めたいなら R=20mm へ落とすことになり、掃引速度は 2/3 になります。`);
  console.log(`  接地率 ${(r.duty*100).toFixed(0)}% → 安定支持には 25〜36 個`);
}
