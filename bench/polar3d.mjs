// sim3d-polar.html を固定刻みでヘッドレス実行して、可視化が物理と食い違っていないかを見る。
//
// **これは論文値の出どころではありません。** 数値は bench/polar.mjs（ヨーを含む6自由度）が
// 正です。3D 側はヨーを持たない（軸平行の箱で描いているため）ので、同じ条件でも
// わずかに違う答えを出します。ここはその差が小さいことを確認するための一致検査です。
//
// 画面の読み出しは瞬時値なので、そこから数字を拾うと位相のどこを見たかで振れます。
// ここは静的釣合いから始めて固定刻みで回し、暖機を捨てた区間の平均を出します。
//
//   node bench/polar3d.mjs
//   CHROMIUM_PATH=/path/to/chrome node bench/polar3d.mjs   （別のバイナリを使うとき）
//
// 各ケースは毎回すべてのパラメータを既定値へ戻してから設定します（試験順に依存しません）。

import { launch } from './browser.mjs';
import { pathToFileURL } from 'url';
import { resolve } from 'path';

const DT = 2e-4, WARM = 2.0, SPAN = 3.0;

const browser = await launch();
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

console.log(`■ sim3d-polar.html 一致検査（dt=${DT}s、静的釣合いから暖機${WARM}s、平均${SPAN}s）`);
console.log(`  荷重は「法線力の合計 / 重量」。1.00 なら定常。`);
console.log(`  数値の出どころは bench/polar.mjs です。ここは可視化が物理と食い違っていないかの確認。`);
console.log(`  各ケースは全パラメータを明示して設定します（a=6mm, δp=0.3mm, A=B=2 が既定）。\n`);

console.log('■ k の向きと Δφ（A=B=2、方位保持なし）');
console.log('  条件                  速度        方向     滑り   荷重          傾き     接地');
for (const [lab, dphiDeg, kdiag] of [
  ['Δφ=60° k=斜め', 60, true], ['Δφ=60° k=軸平行', 60, false],
  ['Δφ=15° k=斜め', 15, true], ['Δφ=15° k=軸平行', 15, false],
]) {
  const r = await run({ dphiDeg, kdiag, a: 6, b: 6, dp: 0.3, shaftA: 2, shaftB: 2 });
  console.log(`  ${lab.padEnd(18)} ${row(r)}`);
  if (r.gap > 2 * r.window) console.log(`    ↑ 最大位相ギャップ ${r.gap.toFixed(0)}° > 接地窓 ${(2*r.window).toFixed(0)}°（支持が切れる条件）`);
}

console.log('\n■ 上下振幅 b（Δφ=60° k=斜め、a=6mm 固定）');
console.log('  b        速度        方向     滑り   荷重          傾き     接地');
for (const b of [6, 3, 1.5, 0.8]) {
  const r = await run({ dphiDeg: 60, kdiag: true, a: 6, dp: 0.3, b, shaftA: 2, shaftB: 2 });
  console.log(`  ${(b + 'mm').padEnd(18)} ${row(r)}`);
}

console.log('\n■ 方位保持（Δφ=60° k=斜め、差動で追い込んでから両軸を揃えて保持）');
console.log('  目標 ψ*   実際の ψ   推力方向    ψ* との差   速度      滑り');
for (const t of [0, 45, 90, 135, 180, 225, 270, 315]) {
  const r = await run({ dphiDeg: 60, kdiag: true, a: 6, b: 6, dp: 0.3, shaftA: 2, shaftB: 2, hold: true, targetDeg: t });
  let e = r.dir - t; while (e > 180) e -= 360; while (e < -180) e += 360;
  let dp = r.psi - t; while (dp > 180) dp -= 360; while (dp < -180) dp += 360;
  console.log(`  ${String(t).padStart(4)}°   ${(r.psi).toFixed(1).padStart(7)}°  ${r.dir.toFixed(1).padStart(7)}°   ` +
    `${e.toFixed(2).padStart(6)}°   ${r.speed.toFixed(0).padStart(4)} mm/s  ${(r.slip * 100).toFixed(0).padStart(3)}%  (ψ誤差 ${dp.toFixed(2)}°)`);
}
console.log('  → 差動を打って ψ を目標へ寄せたあと両軸を同速に戻すので、固定方位のまま走り続ける。');
console.log('    差動を掛けっぱなしにすると ψ が回り続けるだけで、固定方位走行にはならない。');
console.log('\n  ※ 速度が bench/polar.mjs と違うのは、こちらが A=B=2（200rpm）で回しているため。');
console.log('    比で見ること。ヨーが無いぶんの差も残る。数値は bench/polar.mjs を使うこと。');


/* ── 回帰検査 ──────────────────────────────────
   数値だけ見ていると「描画が真っ暗」に気づけない。実際、静的釣合いから始めるように
   したときに rAF の初回 dt=0 で (pen−prev)/dt が 0/0 になり、初期状態が丸ごと NaN に
   なっていたのを、HUD の文字だけ読んでいて見逃した（毎回リセットを押していたため）。
   読み込んだだけの状態で、キャンバスに実際に絵が出ているかを確かめる。 */
console.log('\n■ 回帰検査（読み込んだ直後、操作なし）');
{
  const p2 = await browser.newPage({ viewport: { width: 1280, height: 800 } });
  const bad = [];
  p2.on('pageerror', e => bad.push('例外: ' + e.message));
  p2.on('console', m => { if (m.type() === 'error') bad.push('console.error: ' + m.text()); });
  await p2.goto(pathToFileURL(resolve('sim3d-polar.html')).href);
  await p2.waitForTimeout(1500);

  const r = await p2.evaluate(() => {
    const c = document.getElementById('view');
    const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
    let ink = 0, opaque = 0;
    for (let i = 0; i < d.length; i += 4) {
      /* 透明画素 (0,0,0,0) は背景色と数値が違うので、α を見ないと「未描画」を
         「描画済み」と数えてしまい、真っ白な canvas でも 100% で合格してしまう。 */
      if (d[i + 3] < 250) continue;
      opaque += 1;
      if (Math.abs(d[i] - 10) > 6 || Math.abs(d[i + 1] - 14) > 6 || Math.abs(d[i + 2] - 13) > 6) ink += 1;
    }
    const px = c.width * c.height;
    const txt = id => document.getElementById(id).textContent;
    const num = id => parseFloat(txt(id));
    return { ink: ink / px * 100, opaque: opaque / px * 100,
             nan: ['rV', 'rPsi', 'rTilt'].some(id => /NaN|Infinity/.test(txt(id))),
             v: num('rV'), q: num('rQ') };
  });

  /* 背景そのものが塗られていなければ描画ループが回っていない */
  if (r.opaque < 99) bad.push(`背景が塗られていない（不透明画素 ${r.opaque.toFixed(1)}%）`);
  if (r.ink < 5) bad.push(`キャンバスがほぼ背景のみ（描画 ${r.ink.toFixed(1)}%）`);
  if (r.nan) bad.push('読み出しに NaN / Infinity');
  /* 物理の粗い回帰も見る。既定 A=B=2（200rpm）で a·q̇ = 126 mm/s、その 60〜105% */
  if (!(r.q > 190 && r.q < 210)) bad.push(`波の回転数が想定外（${r.q}）`);
  if (!(r.v > 75 && r.v < 133)) bad.push(`速度が想定範囲外（${r.v} mm/s、期待 75〜133）`);
  console.log(`  描画 ${r.ink.toFixed(1)}% / 不透明 ${r.opaque.toFixed(1)}% / 速度 ${r.v} mm/s / 波 ${r.q} rpm`);
  console.log(bad.length ? '  ✗ ' + bad.join(' / ') : '  ✓ 読み込み直後から描画・読み出しとも正常');
  await p2.close();
  if (bad.length) process.exitCode = 1;
}

await browser.close();
