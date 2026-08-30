// sim3d-polar.html を固定刻みでヘッドレス実行し、周期平均を記録する。
//
// 画面の読み出しは瞬時値なので、そこから数字を拾うと過渡や位相のどこを見たかで振れる。
// ここは静的釣合いから始めて固定刻みで回し、暖機を捨てた区間の平均を出す。
//
//   node bench/polar3d.mjs
//
// 3D 側はヨーを持たない（軸平行の箱で描いているため）。ヨーを含む完全な6自由度は
// bench/polar.mjs。両者は接地点速度 (ω×r) と ψ̇ 項の扱いを揃えてある。

import { chromium } from 'playwright';
import { pathToFileURL } from 'url';
import { resolve } from 'path';

const DT = 2e-4, WARM = 2.0, SPAN = 3.0;

const browser = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome' });
const page = await browser.newPage();
await page.goto(pathToFileURL(resolve('sim3d-polar.html')).href);
await page.waitForFunction(() => !!window.__probe);

const run = (cfg) => page.evaluate(({ cfg, DT, WARM, SPAN }) => {
  window.__probe.config(cfg);
  return window.__probe.run(SPAN, WARM, DT);
}, { cfg, DT, WARM, SPAN });

const row = r => `${r.speed.toFixed(0).padStart(4)} mm/s  ${r.dir.toFixed(1).padStart(7)}°  ` +
  `${r.slip === null ? '  — ' : (r.slip * 100).toFixed(0).padStart(3) + '%'}  ` +
  `荷重 ${r.load.toFixed(2)}±${(r.ripple * 100).toFixed(0)}%  傾き ${r.tilt.toFixed(2).padStart(5)}°  ` +
  `接地 ${r.contact.toFixed(1).padStart(4)}(最小${r.minC})`;

console.log(`■ sim3d-polar.html 固定刻み計測（dt=${DT}s、静的釣合いから暖機${WARM}s、平均${SPAN}s）`);
console.log(`  荷重は「法線力の合計 / 重量」。1.00 なら定常。\n`);

console.log('■ k の向きと Δφ（A=B=2、方位保持なし）');
console.log('  条件                  速度        方向     滑り   荷重          傾き     接地');
for (const [lab, dphiDeg, kdiag] of [
  ['Δφ=60° k=斜め', 60, true], ['Δφ=60° k=軸平行', 60, false],
  ['Δφ=15° k=斜め', 15, true], ['Δφ=15° k=軸平行', 15, false],
]) {
  const r = await run({ dphiDeg, kdiag, shaftA: 2, shaftB: 2 });
  console.log(`  ${lab.padEnd(18)} ${row(r)}`);
  if (r.gap > 2 * r.window) console.log(`    ↑ 最大位相ギャップ ${r.gap.toFixed(0)}° > 接地窓 ${(2*r.window).toFixed(0)}°（支持が切れる条件）`);
}

console.log('\n■ 上下振幅 b（Δφ=60° k=斜め、a=6mm 固定）');
console.log('  b        速度        方向     滑り   荷重          傾き     接地');
for (const b of [6, 3, 1.5, 0.8]) {
  const r = await run({ dphiDeg: 60, kdiag: true, b, shaftA: 2, shaftB: 2 });
  console.log(`  ${(b + 'mm').padEnd(18)} ${row(r)}`);
}

console.log('\n■ 方位保持（Δφ=60° k=斜め、差動で追い込んでから両軸を揃えて保持）');
console.log('  目標 ψ*   実際の ψ   推力方向    ψ* との差   速度      滑り');
for (const t of [0, 45, 90, 135, 180, 225, 270, 315]) {
  const r = await run({ dphiDeg: 60, kdiag: true, shaftA: 2, shaftB: 2, hold: true, targetDeg: t });
  let e = r.dir - t; while (e > 180) e -= 360; while (e < -180) e += 360;
  let dp = r.psi - t; while (dp > 180) dp -= 360; while (dp < -180) dp += 360;
  console.log(`  ${String(t).padStart(4)}°   ${(r.psi).toFixed(1).padStart(7)}°  ${r.dir.toFixed(1).padStart(7)}°   ` +
    `${e.toFixed(2).padStart(6)}°   ${r.speed.toFixed(0).padStart(4)} mm/s  ${(r.slip * 100).toFixed(0).padStart(3)}%  (ψ誤差 ${dp.toFixed(2)}°)`);
}
console.log('  → 差動を打って ψ を目標へ寄せたあと両軸を同速に戻すので、固定方位のまま走り続ける。');
console.log('    「操舵中」プリセット（A=2.5, B=1.5）は差動が掛かりっぱなしで ψ が 50rpm で回り続け、');
console.log('    固定方位走行にはならない。任意方位を出すにはこの保持制御が要る。');

await browser.close();
