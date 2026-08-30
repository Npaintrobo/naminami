// index.html の実機近似モードをヘッドレスで走らせて、荷の搬送速度と操舵の線形性を測る。
//
// docs/hardware-approximation.md §7 の表はここで出しています。
// 実機の測定値ではなく、あくまで index.html に入っている接触モデルの中の値です。
//
//   node bench/table.mjs

import { launch } from './browser.mjs';
import { pathToFileURL } from 'url';
import { resolve } from 'path';

const URL = pathToFileURL(resolve('index.html')).href;

const browser = await launch();
const page = await browser.newPage();
await page.goto(URL);
await page.waitForFunction(() => typeof window.__probe !== 'undefined' || document.readyState === 'complete');

/* ページ内で物理だけを回す。描画も rAF も通さず advanceHardware を直接叩く。 */
async function run({ wx, wy, mu = 0.6, squash = 0.3, retract = true, T = 5, settle = 1 }) {
  return page.evaluate(({ wx, wy, mu, squash, retract, T, settle }) => {
    const h = window.__probe;
    h.reset(mu, retract);
    h.setSquash(squash);
    h.setMotor(wx, wy);
    const dt = 1 / 600;
    let x0 = 0, y0 = 0, t = 0, sx = 0, sy = 0, n = 0;
    for (let i = 0; i < T / dt; i += 1) {
      h.step(dt);
      t += dt;
      if (t >= settle) {
        if (n === 0) { x0 = h.puck().x; y0 = h.puck().y; }
        sx += h.vel().x; sy += h.vel().y; n += 1;
      }
    }
    const p = h.puck();
    return { dx: p.x - x0, dy: p.y - y0, span: T - settle, vx: sx / n, vy: sy / n, tip: h.tip() };
  }, { wx, wy, mu, squash, retract, T, settle });
}

console.log('■ index.html 実機近似モード（n=12, p=50, Δφ=60°, r=6, 300rpm, μ=0.6）');
console.log('  クーロン摩擦モデル。η は使っていません。\n');

console.log('■ 搬送速度');
console.log('  条件                          荷の速度   r·ω 比');
for (const [label, wx, wy, retract] of [
  ['片軸のみ・停止軸リトラクト', 3, 0, true],
  ['片軸のみ・停止軸は接地したまま', 3, 0, false],
  ['両軸運転（45°方向）', 3*Math.SQRT1_2, 3*Math.SQRT1_2, true],
]) {
  const r = await run({ wx, wy, retract });
  const v = Math.hypot(r.dx, r.dy) / r.span;
  console.log(`  ${label.padEnd(28)} ${v.toFixed(1).padStart(6)} mm/s  ${(v / r.tip * 100).toFixed(0).padStart(3)}%`);
}

console.log('\n■ 操舵の線形性（両軸運転・過渡1秒を除いて4秒積分）');
const cmds = [5, 15, 25, 35, 45, 55, 65, 75, 85];
const got = [];
for (const c of cmds) {
  const a = c * Math.PI / 180;
  const r = await run({ wx: 3 * Math.cos(a), wy: 3 * Math.sin(a) });
  got.push(Math.atan2(r.dy, r.dx) * 180 / Math.PI);
}
console.log('  指令   ' + cmds.map(c => String(c).padStart(6)).join(''));
console.log('  結果   ' + got.map(g => g.toFixed(1).padStart(6)).join(''));
console.log('  誤差   ' + got.map((g, i) => (g - cmds[i]).toFixed(1).padStart(6)).join(''));
const err = got.map((g, i) => g - cmds[i]);
console.log(`  最大 |誤差| ${Math.max(...err.map(Math.abs)).toFixed(1)}°  平均 ${(err.reduce((a, b) => a + b, 0) / err.length).toFixed(2)}°（系統的な偏りがあればここに出る）`);

console.log('\n■ パッドの沈み込み δp を振る（これが未実測の最大の不確かさ）');
console.log('  δp[mm]  最大|操舵誤差|   片軸リトラクト時の荷速（r·ω比）');
for (const sq of [0.15, 0.3, 0.6, 1.0, 1.5, 2.0]) {
  const errs = [];
  for (const c of cmds) {
    const a = c * Math.PI / 180;
    const r = await run({ wx: 3 * Math.cos(a), wy: 3 * Math.sin(a), squash: sq });
    errs.push(Math.atan2(r.dy, r.dx) * 180 / Math.PI - c);
  }
  const solo = await run({ wx: 3, wy: 0, squash: sq });
  const v = Math.hypot(solo.dx, solo.dy) / solo.span;
  console.log(`  ${String(sq).padStart(5)}   ${Math.max(...errs.map(Math.abs)).toFixed(1).padStart(6)}°        ${(v / solo.tip * 100).toFixed(0).padStart(3)}%`);
}
console.log('  → このモデルでは δp を 13 倍振っても誤差 9.5〜10.5°、荷速 21% で動かない。');
console.log('    操舵誤差を決めているのは δp ではなく位相ステップ Δφ です（bench/phasestep.mjs）。');

await browser.close();
