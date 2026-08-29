import { simulate, steering, ARCH } from './compare.mjs';

const ARCHS = ['shareXY', 'altXY', 'steer1', 'disc'];
const NAME = Object.fromEntries(ARCHS.map(a => [a, ARCH[a].label]));

console.log('='.repeat(78));
console.log('1. 基準点（μ=0.6, δp=0.30mm, 24本, Δφ=15°, r=6mm, 300rpm, 円盤R=30 α=15°）');
console.log('='.repeat(78));
console.log('  方式                    速度      v/歯先   滑り    操舵誤差');
const base = {};
for (const a of ARCHS) {
  const r = simulate(a);
  const s = steering(a);
  base[a] = { ...r, err: s.maxErr };
  console.log(`  ${NAME[a].padEnd(20)} ${r.speed.toFixed(0).padStart(5)} mm/s  ${r.ratio.toFixed(3)}  ${(r.slip*100).toFixed(0).padStart(3)}%   ${s.maxErr.toFixed(1).padStart(5)}°`);
}

console.log('\n' + '='.repeat(78));
console.log('2. 未測定パラメータを振る — 順位はひっくり返るか');
console.log('='.repeat(78));
console.log('\n[滑り率]  行=パッド沈み込み δp,  列=方式');
console.log('  δp      窓幅    ' + ARCHS.map(a => NAME[a].slice(0,8).padStart(10)).join(''));
for (const comp of [0.06, 0.15, 0.30, 0.60, 1.50]) {
  const win = Math.acos(Math.max(-1, 1-comp/6))*180/Math.PI;
  const cells = ARCHS.map(a => ((simulate(a, { comp }).slip*100).toFixed(0)+'%').padStart(10));
  console.log(`  ${comp.toFixed(2)}mm  ±${win.toFixed(0).padStart(2)}°  ` + cells.join(''));
}

console.log('\n[滑り率]  摩擦係数 μ を振る（δp=0.30mm）');
console.log('  μ      ' + ARCHS.map(a => NAME[a].slice(0,8).padStart(10)).join(''));
for (const mu of [0.3, 0.6, 1.0]) {
  const cells = ARCHS.map(a => ((simulate(a, { mu }).slip*100).toFixed(0)+'%').padStart(10));
  console.log(`  ${mu.toFixed(1)}   ` + cells.join(''));
}

console.log('\n[操舵誤差]  パッド沈み込みを振る');
console.log('  δp      ' + ARCHS.map(a => NAME[a].slice(0,8).padStart(10)).join(''));
for (const comp of [0.15, 0.30, 0.60, 1.50]) {
  const cells = ARCHS.map(a => (steering(a, { comp }).maxErr.toFixed(1)+'°').padStart(10));
  console.log(`  ${comp.toFixed(2)}mm  ` + cells.join(''));
}

console.log('\n' + '='.repeat(78));
console.log('3. 目標速度を出すために何を払うか（ここが効く）');
console.log('='.repeat(78));
console.log('速度はクランク式では v = ratio·r·ω なので、速くするには r（ストローク）を上げるしかない。');
console.log('振動 a = (2π)⁴n²r/16 も皮膜ひずみ ε = 2πr/L も r で決まるので、速度と一緒に悪化する。');
console.log('円盤式は v = ratio·e·ω で e は円盤半径。r は接地クリアランスだけで決めてよい。\n');

const W = 300*Math.PI/30, L = 1500, N_WAVE = 1;
for (const V of [300, 600, 1000]) {
  console.log(`■ 目標 ${V} mm/s（300 rpm 固定）`);
  console.log('  方式                    必要な r     ストローク   振動      皮膜ひずみ');
  for (const a of ARCHS) {
    const ratio = base[a].ratio;
    let r, note = '';
    if (a === 'disc') { r = 6; note = `（円盤 e=${(V/(ratio*W)).toFixed(0)}mm）`; }
    else r = V/(ratio*W);
    const acc = Math.pow(2*Math.PI,4)*Math.pow(300/60,2)*(r/1000)/16/9.81;
    const eps = 2*Math.PI*N_WAVE*r/L*100;
    console.log(`  ${NAME[a].padEnd(20)} ${r.toFixed(1).padStart(5)} mm  ${(2*r).toFixed(0).padStart(6)} mm  ${acc.toFixed(1).padStart(6)} G  ${eps.toFixed(1).padStart(6)} %  ${note}`);
  }
  console.log('');
}

console.log('='.repeat(78));
console.log('4. 設計上の性質（シミュレーションでは出ない項目）');
console.log('='.repeat(78));
console.log('  方式                  アクチュエータ  ラティス  定常横行  ヨー  速度とストローク');
for (const a of ARCHS) {
  const A = ARCH[a];
  console.log(`  ${NAME[a].padEnd(20)} ${String(A.actuators).padStart(8)}      ${A.rails}        ` +
    `${(A.lateral === true ? '可' : A.lateral === 'slew' ? '要旋回' : '不可').padEnd(6)}  ${A.yaw ? '可' : '不可'}  ${A.decoupled ? '分離' : '連動'}`);
}
