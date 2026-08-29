// 直交2ラティスの操舵誤差を決めているのは何か。
//
// 「卓上（荷が数列だけ）なら成立、全面接地（車両）だと壊れる」と書いていましたが誤りでした。
// 支配しているのは接地の局所性ではなく**位相ステップ Δφ** です。
//
// 過拘束（X歯は vx=r·ωx かつ vz=0、Y歯は vz=r·ωy かつ vx=0 を要求する）が実際に害になるのは、
// 両ラティスの歯が「同時に」荷重を持っているときだけです。同時に接地する歯数は
//
//     同時接地歯数 ≒ 接地窓の角度幅 / Δφ,   接地窓 = 2·acos(1 − δp/r)
//
// で決まります。これが 1 を超えると両方の拘束が常時同時に効き、クーロン摩擦は滑り速度に
// 依らないので遅いほうがフルの摩擦力でブレーキになり、進行方向が 45° 方向へ引き寄せられます。
//
//   node bench/phasestep.mjs

import { simulate, steering } from './compare.mjs';

const r = 6, rpm = 300, tip = r * rpm * Math.PI / 30;
const teeth = (dphi, comp) => 2 * Math.acos(Math.max(-1, 1 - comp / r)) * 180 / Math.PI / dphi;

console.log('■ Δφ を振る（r=6mm, δp=0.3mm, μ=0.6, 300rpm, 全レール接地）');
console.log('  Δφ    λ=(360/Δφ)·p  同時接地歯数  最大操舵誤差  45°方向の荷速(r·ω比)');
for (const dphi of [120, 90, 60, 45, 30, 20, 15]) {
  const rails = Math.max(3, Math.round(360 / dphi) * 2);
  const o = { rails, dphiDeg: dphi, crank: r, rpm, comp: 0.3 };
  const st = steering('shareXY', o);
  const sp = simulate('shareXY', o, Math.PI / 4);
  console.log(`  ${String(dphi).padStart(3)}°   ${String(Math.round(360 / dphi)).padStart(2)}·p        ` +
    `${teeth(dphi, 0.3).toFixed(2).padStart(6)}      ${st.maxErr.toFixed(1).padStart(5)}°        ` +
    `${(sp.speed / tip * 100).toFixed(0).padStart(3)}%`);
}
console.log('  → 誤差は Δφ で 3.3°〜46.9° と一桁動くのに、荷速は 41〜50% でほぼ動かない。');
console.log('    Δφ を粗くするのは速度を犠牲にしない。ただ波長が短くなる（下記のジレンマ）。\n');

console.log('■ Δφ と δp を両方振って、同時接地歯数で並べ直す');
console.log('  同時接地歯数   Δφ    δp      最大操舵誤差');
const rows = [];
for (const dphi of [120, 90, 60, 45, 30, 20, 15]) {
  for (const comp of [0.15, 0.3, 0.6, 1.2]) {
    const rails = Math.max(3, Math.round(360 / dphi) * 2);
    rows.push({ t: teeth(dphi, comp), dphi, comp, e: steering('shareXY', { rails, dphiDeg: dphi, crank: r, rpm, comp }).maxErr });
  }
}
rows.sort((a, b) => a.t - b.t);
for (const x of rows) {
  console.log(`  ${x.t.toFixed(2).padStart(9)}   ${String(x.dphi).padStart(3)}°  ${String(x.comp).padStart(4)}mm   ${x.e.toFixed(1).padStart(5)}°`);
}
console.log(`
  → 誤差は Δφ と δp のどちらか一方ではなく、その比＝同時接地歯数でおおむね整理できます。
    ばらつきは残る（0.86 のところで 14.9° と 28.7°）ので、単一パラメータではありません。

    同時接地歯数 < 0.5  → 誤差 1〜7°     使える
    0.5 〜 1.0          → 誤差 3〜29°    δp の実測が要る
    1.0 〜 2.5          → 誤差 20〜47°   壊れている
    > 2.5               → 誤差 22〜36°   飽和（常時両方が効くので、もう悪くならない）

■ 設計上のジレンマ

  操舵は Δφ を粗くしたい。λ = (360/Δφ)·p なので、Δφ ≥ 60° は λ ≤ 6·p を意味します。
  一方 bench/spatial.mjs のケースB が示すとおり、有限寸法の荷を支え続けるには
  波長が荷の差し渡しより十分長いことが要ります（λ が荷の 2 倍では 0.83 秒で転倒した）。

  つまり     荷の差し渡し  ≪  λ  =  (360/Δφ)·p  ≤  6·p

  ピッチ 50mm なら λ ≤ 300mm、支えられる荷はせいぜい 100mm 級です。
  600mm の車体は λ の 2 倍あり、この不等式を満たしようがありません。
  **直交2ラティスが卓上で成立して車両で壊れる本当の理由はこれです。**
  「全面接地だから」ではありません（全面接地でも Δφ=60° なら誤差 3.3°）。`);
