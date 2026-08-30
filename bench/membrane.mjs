// 上下だけの進行波を、厚みのあるゴム膜に通したら推力は出るのか。
//
// このプロジェクトは最初に「上下だけでは運べない（円軌道が要る）」と結論した。
// しかしあれは**剛体の歯**の話。厚みのある弾性膜なら、曲げ波の表面粒子は楕円を描く。
// 中立面から h 離れた表面の面内変位は
//
//     u = −h·∂w/∂x           w = W·cos(φ),  φ = k·s − Ω·t
//       = h·W·k·sin(φ)
//     u̇ = −h·W·k·Ω·cos(φ)    山（φ=0）で −h·W·k·Ω
//
// つまり山での面内速度は  a_eff = h·W·k  に相当する。超音波進行波モータと同じ原理。
// これが本当なら、偏極操舵で必要だった「足ごとの面内駆動」が丸ごと要らなくなる。
//
//   node bench/membrane.mjs

import { lattice, tune, core, TIP, AMP, P, N } from './polar.mjs';

/* 膜は連続体なので、接触の離散化は接触子ピッチと無関係に細かく取る必要がある。
   ピッチ50mmの格子で λ=100mm の波を標本化すると1波長2点しかなく（ナイキスト限界）、
   山の高さも面内速度も正しく出ない。盤面 600mm 角を pitch mm 刻みで張る。 */
/* 膜の接地剛性。位相場が空間で連続なので tune() の二分法は要らない。
   振幅 W の余弦波が δp 沈むとき、接地するのは |φ| < acos(1−δp/W) の範囲。
   平均貫入はその半分程度として K を置く。 */
function membraneK(pads, Wamp, dp = 0.3, mass = 20) {
  const frac = Math.acos(Math.max(-1, 1 - dp / Wamp)) / Math.PI;
  const nc = Math.max(1, pads.length * frac);
  const K = (mass * 9810 / 1000) / (nc * dp * 0.5);
  const C = 2 * 0.35 * Math.sqrt(nc * K * 1000 * mass) / 1000 / nc;
  return { K, C, nc };
}

function fineGrid(pitch) {
  const n = Math.round(N * P / pitch);
  const pads = [];
  for (let i = 0; i < n; i += 1) for (let k = 0; k < n; k += 1) {
    pads.push({ x: (i - (n - 1) / 2) * pitch, z: (k - (n - 1) / 2) * pitch,
                ph: 0, xy: (i + k) % 2, n: 0, prev: 0, dz: 0 });
  }
  return pads;
}

const RPM = 300, OMEGA = RPM * Math.PI / 30;

/* 膜の表面点を、接触子と同じインタフェースで返す駆動。
   height は下向き正の慣習に合わせて符号を反転。 */
function membraneDrive({ h, Wamp, lam, psi = 0, rpm = RPM }) {
  const k = 2 * Math.PI / lam;
  const om = rpm * Math.PI / 30;
  let t = 0;
  const c = Math.cos(psi), s = Math.sin(psi);
  const d = {
    reset() { t = 0; d.psi = psi; },
    step(tt) { t = tt; },
    foot(p) {
      const sPos = p.x * c + p.z * s;              // 波の進む向きの座標
      const phi = k * sPos - om * t;
      const cs = Math.cos(phi), sn = Math.sin(phi);
      const off = h * Wamp * k * sn;               // 面内変位（ψ方向）
      const vel = -h * Wamp * k * om * cs;         // 面内速度
      return [-Wamp * cs, off * c, off * s, vel * c, vel * s];
    },
    psi
  };
  return d;
}

const aEff = (h, Wamp, lam) => h * Wamp * 2 * Math.PI / lam;
const strain = (h, Wamp, lam) => h * Wamp * (2 * Math.PI / lam) ** 2;

console.log('■ 1. 上下だけの波は、膜を通せば推力を出すか\n');
console.log('  厚み 2h のゴム膜。山での面内速度は理論上 a_eff = h·W·k。');
console.log('  剛体の歯（面内駆動なし）と比べる。\n');
console.log('  条件                                 a_eff    予測速度   実測(経路)   滑り    膜のひずみ');
{
  const pads = lattice(60, { kdir: 'diag' });
  const t = tune(pads, 6);
  void t;

  /* 比較用: 上下だけ・面内ゼロ（＝剛体の歯） */
  const rigid = {
    reset() { rigid.psi = 0; }, step() {},
    foot(p) { return [-6 * Math.cos(p.ph), 0, 0, 0, 0]; }, psi: 0
  };
  const r0 = core(pads, rigid, { b: 6, K: t.K, C: t.C });
  console.log(`  剛体の歯・上下のみ（面内駆動なし）        0.00 mm   0 mm/s   ${r0.pathSpeed.toFixed(1).padStart(6)} mm/s   ` +
    `${(r0.slip * 100).toFixed(0).padStart(3)}%       —`);

  for (const [h, Wamp, lam] of [[3, 3, 300], [3, 5, 300], [3, 5, 150], [3, 5, 100], [5, 6, 100], [8, 8, 100]]) {
    const a = aEff(h, Wamp, lam);
    const pred = a * OMEGA;
    /* 膜の波は空間座標で決まる。1波長あたり最低12点になるよう刻む。 */
    const flat = fineGrid(Math.max(12.5, lam / 10));
    const tt = membraneK(flat, Wamp);
    const r = core(flat, membraneDrive({ h, Wamp, lam }), { b: Wamp, K: tt.K, C: tt.C });
    console.log(`  h=${String(h).padStart(2)}mm W=${String(Wamp).padStart(2)}mm λ=${String(lam).padStart(3)}mm      ` +
      `${a.toFixed(2).padStart(5)} mm  ${pred.toFixed(0).padStart(4)} mm/s   ${r.pathSpeed.toFixed(1).padStart(6)} mm/s   ` +
      `${(r.slip * 100).toFixed(0).padStart(3)}%     ${(strain(h, Wamp, lam) * 100).toFixed(1).padStart(5)}%`);
  }
}

console.log(`
  → 剛体の歯は上下だけでは動かない（このプロジェクトの最初の結論どおり）。
    **膜を挟むと動く。** 面内駆動の機構が一切ないのに推力が出る。
    ただし a_eff = h·W·k は小さい。同じ 300rpm では偏極操舵の a = ${AMP}mm
    （→ ${TIP.toFixed(0)} mm/s）に対して 1/3 以下。ただし速度は rpm に比例するので、
    回転数で埋められるかは次節で見る。`);

console.log('\n■ 2. 何が速度を決めるか\n');
{
  console.log('  a_eff = h·W·(2π/λ) なので、厚く・大振幅・短波長ほど速い。');
  console.log('  ただし膜のひずみ ε = h·W·k² も同じ方向に増える。\n');
  console.log('  λ[mm]   h=3,W=5      h=5,W=6      h=8,W=8       （上段 a_eff、下段 ひずみ）');
  for (const lam of [300, 200, 150, 100, 60]) {
    const row = [[3, 5], [5, 6], [8, 8]].map(([h, w]) =>
      `${aEff(h, w, lam).toFixed(2)}mm/${(strain(h, w, lam) * 100).toFixed(0)}%`.padEnd(13));
    console.log(`  ${String(lam).padStart(3)}     ${row.join('')}`);
  }
  console.log(`
  → ひずみ 10% を上限とすると λ=100mm で a_eff は 1mm 前後。300rpm で 30 mm/s、
    900rpm で 90 mm/s。**卓上コンベアの速度域**で、車両の速度域ではない。`);
}

console.log('\n■ 2b. 回転数で速度を埋められるか\n');
{
  /* a_eff は幾何だけで決まり回転数に依らないので、速度は rpm に比例する。
     偏心式と違いストロークを増やさずに速度を上げられる（振動は W·ω² で増える）。 */
  console.log('  h=5mm W=6mm λ=100mm（ひずみ 11.8% 固定）');
  console.log('    rpm    速度(経路)   滑り    上下加速度 W·ω²');
  const pads = fineGrid(12.5);
  const tt = membraneK(pads, 6);
  for (const rpm of [300, 600, 900, 1200]) {
    const om = rpm * Math.PI / 30;
    const r = core(pads, membraneDrive({ h: 5, Wamp: 6, lam: 100, rpm }), { b: 6, K: tt.K, C: tt.C });
    console.log(`    ${String(rpm).padStart(4)}   ${r.pathSpeed.toFixed(0).padStart(5)} mm/s  ` +
      `${(r.slip * 100).toFixed(0).padStart(3)}%      ${(6 * om * om / 9810).toFixed(1).padStart(5)} G`);
  }
  console.log(`
  → 速度は rpm に比例する。**膜のひずみは rpm に依らない**（幾何だけで決まる）ので、
    偏心式のように「速度を上げるとストロークが伸びて振動と皮膜ひずみが悪化する」
    という結合がない。ここが偏極操舵との決定的な差。`);
}

console.log('\n■ 3. 操舵は効くか（波の向きを変えるだけ）\n');
{
  console.log('  接触子は動かさず、波の向き ψ だけを変える。');
  console.log('    ψ*      推力方向    誤差     速度(経路)   滑り');
  const pads = fineGrid(12.5);
  for (const psiDeg of [0, 45, 90, 135, 180, 270]) {
    const tt = membraneK(pads, 5);
    const r = core(pads, membraneDrive({ h: 5, Wamp: 5, lam: 100, psi: psiDeg * Math.PI / 180 }),
                   { b: 5, K: tt.K, C: tt.C });
    let e = r.dir - psiDeg; while (e > 180) e -= 360; while (e < -180) e += 360;
    console.log(`    ${String(psiDeg).padStart(3)}°   ${r.dir.toFixed(1).padStart(7)}°  ${e.toFixed(2).padStart(6)}°   ` +
      `${r.pathSpeed.toFixed(1).padStart(5)} mm/s  ${(r.slip * 100).toFixed(0).padStart(3)}%`);
  }
  console.log(`
  → 波の向きを変えるだけで推力方向が追従する。接触子は一切動かない。
    偏極操舵で必要だった「足ごとの面内駆動」も「全足共通の ψ 軸」も要らない。`);
}

console.log(`
■ 4. 機構として何が要るか

  波は  w = W·cos(k·s − Ω·t)  。位相 k·s は**空間座標だけ**で決まる。
  つまり足ごとの位相クロッキングが要らない。螺旋がその位相配りを丸ごと担う。

  SAW（Single Actuator Wave、Zarrouk ら）の構成そのもの:
    ・ピッチ p の螺旋を1本、波の進む向きに寝かせて回す
    ・螺旋に乗ったフォロワが、位置 s に応じた位相で上下する
    ・波長 λ = p、波速 = p·Ω/2π

  操舵は、螺旋アセンブリごと鉛直軸まわりに回す。**膜も接地面も回らない**
  （波は変形パターンであって 材料の回転ではないため）。

  部品の見積り（直径600mm、λ=100mm、レール間隔25mm）:
    ・螺旋 1本 ＋ レール 24本 ＋ フォロワ 24個 ＋ 操舵の回転台 1  =  **約50部品**
    ・偏極操舵の 432〜1296 部品に対して 1桁少ない

■ 未検証（重要）

  ・膜の力学を「中立面から h の点の面内変位 = −h·∂w/∂x」という薄板の線形式で
    置いただけ。ゴムの大変形・非圧縮性・接触での平坦化は入っていない。
    ひずみ 10% はこの線形式が怪しくなる領域。
  ・膜と接地面の間の実際の接触（面圧分布、滑り域と固着域の分離）を解いていない。
    超音波モータではここが効率を決める。
  ・螺旋→レール→膜 の伝達（フォロワの摩擦、レールの案内、螺旋の軸力）は未検討。
  ・**SAW は既存**。Zarrouk らの単一アクチュエータ波動ロボットが先行技術。
    残り得るのは「波の向きを steer する」構成だが、未調査。`);
