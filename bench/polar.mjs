// 和差位相式・偏極操舵波動面を、案A（直交2ラティス）と同じ物理コアで比べる。
//
//   θ_i = k·r_i − q          搬送波位相。どの接触子が荷重を受けるかを決める
//   z_i = b·cos θ_i          上下
//   p_i = a·sin θ_i ·[cosψ, sinψ]   面内。ψ が偏極角＝推力方向
//
// 面内速度 = d/dt p_i = −a·q̇·cos θ_i·[cosψ, sinψ] + a·ψ̇·sin θ_i·[−sinψ, cosψ]
//   第1項: 山（θ=0）でちょうど a·q̇。案Aの r·ω に対応する
//   第2項: 操舵中だけ出る横成分。山では sin θ = 0 なので消える
//
// 逆旋回2偏心の合成は u₊(θ+ψ)+u₋(θ−ψ) = 2cos θ·[cosψ,sinψ] で cos になる。
// z と同位相になって推力ゼロなので、両偏心を −π/2 クロッキングして sin にしてある。
//
// 案A との構造的な違い: 接地している接触子の要求速度が全部 ψ 方向で揃う。
// 案Aは X が (r·ωx, 0)、Y が (0, r·ωy) を要求する過拘束だった。

const N = 12, P = 50, AMP = 6, RPM = 300;
const BREF = { v: 6 };        // 上下振幅 b（試験で振る）
const W = RPM * Math.PI / 30;
const TIP = AMP * W;                       // 山での面内速度 = 案Aの r·ω と同じ 188 mm/s
const MU = 0.6, G = 9810, VEPS = 6, ACC = 1000;
const DT = 2e-5, WARM = 1.0;
const TREF = { v: 3.0 };

const TAU = 2 * Math.PI;

/* ── 接触子の配置と位相場 ──────────────────────────── */
function lattice(dphiDeg, { shuffle = false, seed = 7, kdir = 'x' } = {}) {
  const d = dphiDeg * Math.PI / 180;
  const pads = [];
  /* 波数ベクトル k の向き。'x' は軸に平行 → 同位相の接触子が1列に並んで線接触になり、
     機体がその線まわりに首を振る。'diag' は斜め → 同位相の接触子が反対角線に散る。 */
  const step = (i, k) => (kdir === 'diag' ? i + k : i);
  for (let i = 0; i < N; i += 1) for (let k = 0; k < N; k += 1) {
    pads.push({ x: (i - (N - 1) / 2) * P, z: (k - (N - 1) / 2) * P, ph: step(i, k) * d, xy: (i + k) % 2 });
  }
  if (shuffle) {
    let s = seed;
    const rnd = () => (s = (s * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    const ph = pads.map(p => p.ph);
    for (let i = ph.length - 1; i > 0; i -= 1) { const j = Math.floor(rnd() * (i + 1)); [ph[i], ph[j]] = [ph[j], ph[i]]; }
    pads.forEach((p, i) => { p.ph = ph[i]; });
  }
  return pads;
}

/* ── 6自由度コア。kin(pad, t) が [高さ, 面内x, 面内z, 速度x, 速度z] を返す ── */
function core(pads, kin, { M = 20, HCG = 100, K = 47, C = 0.15, load = null } = {}) {
  const L = N * P;
  const Ih = M * (L * L / 12 + HCG * HCG), Iy = M * (L * L / 6);
  let X = 0, Z = 0, Y = BREF.v + 0.4, psi = 0, tx = 0, tz = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  const prev = new Float64Array(pads.length);
  let dx = 0, dz = 0, slip = 0, wt = 0, nvar = 0, nmean = 0, minC = 1e9;

  for (let s = 0; s < TREF.v / DT; s += 1) {
    const t = s * DT, live = t > WARM;
    const cp = Math.cos(psi), sp = Math.sin(psi);
    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0, sl = 0, wsum = 0;

    for (let j = 0; j < pads.length; j += 1) {
      const p = pads[j];
      const [h, ox, oz, ux, uz] = kin(p, t);
      const ax = p.x + ox, az = p.z + oz;
      const rx = ax * cp - az * sp, rz = ax * sp + az * cp;
      const pen = -(Y - tx * rz + tz * rx + h);
      if (pen <= 0) { prev[j] = 0; continue; }
      const Nf = Math.max(0, K * pen + C * (pen - prev[j]) / DT);
      prev[j] = pen;
      if (Nf <= 0) continue;
      cnt += 1; SN += Nf;
      const wx_ = ux * cp - uz * sp, wz_ = ux * sp + uz * cp;
      const gx = vX + vp * rz + wx_, gz = vZ - vp * rx + wz_;
      const m = Math.max(Math.hypot(gx, gz), VEPS);
      const fx = -MU * Nf * gx / m, fz = -MU * Nf * gz / m;
      Fx += fx; Fz += fz;
      Ty += rz * fx - rx * fz;
      Tx += -HCG * fz - rz * Nf;
      Tz += rx * Nf + HCG * fx;
      sl += Nf * Math.hypot(gx, gz); wsum += Nf;
    }

    vX += ACC * Fx / M * DT; vZ += ACC * Fz / M * DT;
    vY += (ACC * SN / M - G) * DT;
    vp += ACC * Ty / Iy * DT; wx += ACC * Tx / Ih * DT; wz += ACC * Tz / Ih * DT;
    X += vX * DT; Z += vZ * DT; Y += vY * DT; psi += vp * DT; tx += wx * DT; tz += wz * DT;

    if (live) {
      dx += vX * DT; dz += vZ * DT; wt += DT;
      slip += (wsum > 0 ? sl / wsum : 0) * DT;
      nmean += SN * DT; nvar += SN * SN * DT;
      if (cnt < minC) minC = cnt;
    }
  }
  const mn = nmean / wt;
  void load;
  return {
    speed: Math.hypot(dx, dz) / wt, dir: Math.atan2(dz, dx) * 180 / Math.PI,
    slip: slip / wt / TIP, minC,
    nripple: Math.sqrt(Math.max(0, nvar / wt - mn * mn)) / Math.max(mn, 1e-9)
  };
}

/* ── 2つの運動学 ──────────────────────────────── */
const polar = (pads, psiOf) => core(pads, (p, t) => {
  const th = p.ph - W * t, ps = psiOf(t), c = Math.cos(ps), s = Math.sin(ps);
  const sn = Math.sin(th), cs = Math.cos(th);
  return [-BREF.v * cs, AMP * sn * c, AMP * sn * s, -AMP * W * cs * c, -AMP * W * cs * s];
});

const shareXYr = (pads, cmd, retract) => {
  const wx = W * Math.cos(cmd), wz = W * Math.sin(cmd);
  return core(pads, (p, t) => {
    const lift = a => (retract && Math.abs(a) < 1e-9 * W ? 2 * BREF.v : 0);
    if (p.xy === 0) { const th = p.ph - wx * t; return [-BREF.v * Math.cos(th) + lift(wx), 0, 0, -AMP * wx * Math.cos(th), 0]; }
    const th = p.ph - wz * t; return [-BREF.v * Math.cos(th) + lift(wz), 0, 0, 0, -AMP * wz * Math.cos(th)];
  });
};
const shareXY = (pads, cmd) => shareXYr(pads, cmd, true);

/* ── 試験1: ψ に推力方向が追従するか ── */
function sweep(label, dphi, run, kdir = 'x') {
  const pads = lattice(dphi, { kdir });
  const cmds = [0, 45, 90, 135, 180, 225, 270, 315];
  const rows = cmds.map(c => ({ c, r: run(pads, c * Math.PI / 180) }));
  let worst = 0;
  console.log(`  ${label}`);
  console.log('    指令    推力方向    誤差    速度       滑り  最小接地  荷重リップル');
  for (const { c, r } of rows) {
    let e = r.dir - c; while (e > 180) e -= 360; while (e < -180) e += 360;
    worst = Math.max(worst, Math.abs(e));
    console.log(`    ${String(c).padStart(3)}°  ${r.dir.toFixed(1).padStart(7)}°  ${e.toFixed(2).padStart(6)}°  ` +
      `${r.speed.toFixed(0).padStart(4)} mm/s  ${(r.slip * 100).toFixed(0).padStart(3)}%  ` +
      `${String(r.minC).padStart(4)}    ${(r.nripple * 100).toFixed(0).padStart(3)}%`);
  }
  console.log(`    → 最大操舵誤差 ${worst.toFixed(2)}°\n`);
  return worst;
}

console.log(`■ 偏極操舵 vs 直交2ラティス（12×12接触子, ピッチ50mm, b=a=6mm, 300rpm, μ=0.6）`);
console.log(`  山での面内速度 a·q̇ = ${TIP.toFixed(0)} mm/s（案Aの r·ω と同値）`);
console.log(`  速度は「r·ω 比」で読めます。接触子は両方式とも 144 個。`);
console.log(`  案Aはそれを市松で X/Y に折半、偏極式は全部が同じ1群です。\n`);

console.log('■ 試験1: 推力方向が指令に追従するか（案Aは停止軸リトラクトあり）');
const R = {};
for (const [tag, dphi, kdir] of [['60x', 60, 'x'], ['60d', 60, 'diag'], ['15d', 15, 'diag'], ['15x', 15, 'x']]) {
  R['P' + tag] = sweep(`偏極操舵 Δφ=${dphi}° k=${kdir}`, dphi, (p, c) => polar(p, () => c), kdir);
}
for (const [tag, dphi, kdir] of [['60d', 60, 'diag'], ['15d', 15, 'diag']]) {
  R['S' + tag] = sweep(`直交2ラティス Δφ=${dphi}° k=${kdir}`, dphi, shareXY, kdir);
}

/* ── 試験2: 空間位相をシャッフルすると推力が落ちるか（反証可能な波動性の判定） ── */
console.log('■ 試験2: 空間位相をシャッフルすると推力が落ちるか');
console.log('  ケースA 車両（12×12接触子が大きな剛体の下に全部入る）Δφ=60°, k=斜め, ψ=45°');
console.log('    位相の配り方        速度      方向     滑り  最小接地  荷重リップル');
for (const [lab, sh] of [['空間波（勾配あり）', false], ['同じ位相集合をシャッフル', true]]) {
  const r = polar(lattice(60, { kdir: 'diag', shuffle: sh }), () => Math.PI / 4);
  console.log(`    ${lab.padEnd(20)} ${r.speed.toFixed(0).padStart(4)} mm/s  ${r.dir.toFixed(1).padStart(6)}°  ` +
    `${(r.slip * 100).toFixed(0).padStart(3)}%  ${String(r.minC).padStart(4)}    ${(r.nripple * 100).toFixed(0).padStart(3)}%`);
}

/* ケースB 卓上: 接触子は地面に固定、有限寸法の荷が渡る（局所支持） */
function tableRun(pads, psi, { S = 150, M = 0.5, HCG = 15, K = 2.2, C = 0.006 } = {}) {
  const Ih = M * (S * S / 12 + HCG * HCG), Iy = M * (S * S / 6);
  const TIPPED = 8 * Math.PI / 180;
  /* ラティスは周期 N·P で無限に続くとみなし、格子添字を mod N で引く。
     pads は i*N+k の順に生成してあるのでそのまま添字で取れる。 */
  const at = (i, k) => pads[(((i % N) + N) % N) * N + (((k % N) + N) % N)].ph;
  let X = 0, Z = 0, Y = BREF.v - 0.15, ps = 0, tx = 0, tz = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  const prev = new Map();
  let dx = 0, dz = 0, wt = 0, drop = 0, tilt2 = 0;
  const c = Math.cos(psi), sn = Math.sin(psi);

  for (let s = 0; s < TREF.v / DT; s += 1) {
    const t = s * DT, live = t > WARM;
    if (Math.hypot(tx, tz) > TIPPED) return { fell: t };
    const cp = Math.cos(ps), spz = Math.sin(ps);
    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0;
    const i0 = Math.floor((X - S / 2) / P), i1 = Math.ceil((X + S / 2) / P);
    const k0 = Math.floor((Z - S / 2) / P), k1 = Math.ceil((Z + S / 2) / P);
    for (let i = i0; i <= i1; i += 1) for (let k = k0; k <= k1; k += 1) {
      const rx = i * P - X, rz = k * P - Z;
      const bx = rx * cp + rz * spz, bz = -rx * spz + rz * cp;
      if (Math.abs(bx) > S / 2 || Math.abs(bz) > S / 2) continue;
      const th = at(i, k) - W * t;
      const top = BREF.v + BREF.v * Math.cos(th);                 // 接触子先端の世界高さ
      const pen = top - (Y - tx * rz + tz * rx);
      const key = i * 1e5 + k;
      if (pen <= 0) { prev.set(key, 0); continue; }
      const Nf = Math.max(0, K * pen + C * (pen - (prev.get(key) || 0)) / DT);
      prev.set(key, pen);
      if (Nf <= 0) continue;
      cnt += 1; SN += Nf;
      const ux = -AMP * W * Math.cos(th) * c, uz = -AMP * W * Math.cos(th) * sn;
      const gx = vX + vp * rz - ux, gz = vZ - vp * rx - uz;   // 荷 − 接触子
      const m = Math.max(Math.hypot(gx, gz), VEPS);
      const fx = -MU * Nf * gx / m, fz = -MU * Nf * gz / m;
      Fx += fx; Fz += fz; Ty += rz * fx - rx * fz;
      Tx += -HCG * fz - rz * Nf; Tz += rx * Nf + HCG * fx;
    }
    vX += ACC * Fx / M * DT; vZ += ACC * Fz / M * DT;
    vY += (ACC * SN / M - G) * DT;
    vp += ACC * Ty / Iy * DT; wx += ACC * Tx / Ih * DT; wz += ACC * Tz / Ih * DT;
    X += vX * DT; Z += vZ * DT; Y += vY * DT; ps += vp * DT; tx += wx * DT; tz += wz * DT;
    if (Y < -2) { Y = -2; if (vY < 0) vY = 0; }
    if (live) { dx += vX * DT; dz += vZ * DT; wt += DT; tilt2 += (tx * tx + tz * tz) * DT; if (cnt === 0) drop += DT; }
  }
  return { speed: Math.hypot(dx, dz) / wt, dir: Math.atan2(dz, dx) * 180 / Math.PI,
           drop: drop / wt * 100, tilt: Math.sqrt(tilt2 / wt) * 180 / Math.PI };
}

console.log('\n  ケースB 卓上（接触子は地面に固定、150mm角の荷が渡る。下に入るのは3×4個）Δφ=60°, k=斜め, ψ=45°');
console.log('    位相の配り方        結果');
for (const [lab, sh] of [['空間波（勾配あり）', false], ['同じ位相集合をシャッフル', true]]) {
  const r = tableRun(lattice(60, { kdir: 'diag', shuffle: sh }), Math.PI / 4);
  console.log(`    ${lab.padEnd(20)} ` + (r.fell
    ? `荷が転倒（t=${r.fell.toFixed(2)}s）`
    : `${r.speed.toFixed(0)} mm/s  ${r.dir.toFixed(1)}°  無接地 ${r.drop.toFixed(1)}%  傾き ${r.tilt.toFixed(2)}°`));
}

/* ── 試験3: 操舵の過渡。ψ を走行中に振ったとき荷重と滑りが暴れないか ── */
console.log('\n■ 試験3: 走行中に ψ を 0°→90° へ振る（Δφ=60°, k=斜め）');
console.log('    操舵レート    速度      方向     滑り   荷重リップル');
for (const rate of [0, 45, 180, 720]) {
  const psiOf = t => Math.min(Math.PI / 2, Math.max(0, (t - WARM)) * rate * Math.PI / 180);
  const r = polar(lattice(60, { kdir: 'diag' }), psiOf);
  console.log(`    ${String(rate).padStart(4)}°/s   ${r.speed.toFixed(0).padStart(4)} mm/s  ${r.dir.toFixed(1).padStart(6)}°  ` +
    `${(r.slip * 100).toFixed(0).padStart(3)}%   ${(r.nripple * 100).toFixed(0).padStart(3)}%`);
}
console.log('    ※ 0°/s は ψ=0 固定の基準。720°/s は 0.125 秒で 90° 振る想定。\n');

console.log('■ まとめ: 最大操舵誤差');
console.log('  方式               Δφ=60° k斜め   Δφ=15° k斜め   Δφ=60° k軸平行  Δφ=15° k軸平行');
console.log(`  偏極操舵          ${R.P60d.toFixed(2).padStart(7)}°      ${R.P15d.toFixed(2).padStart(7)}°       ${R.P60x.toFixed(2).padStart(7)}°       ${R.P15x.toFixed(2).padStart(7)}°`);
console.log(`  直交2ラティス      ${R.S60d.toFixed(2).padStart(7)}°      ${R.S15d.toFixed(2).padStart(7)}°`);

/* ── 追加: 卓上（有限の荷）で Δφ を振る。案Aは Δφ を細かくすると操舵が壊れたが、
      偏極式は過拘束がないので細かくできるはず。それが本当なら案Aのジレンマが解ける。 ── */
TREF.v = 2.0;
console.log('\n■ 卓上（150mm角の荷）で Δφ を振る（k=斜め, ψ=45°）。λ と荷の大小に注目');
console.log('    Δφ    波長(対角)   空間波                    シャッフル');
for (const dphi of [180, 120, 90, 60, 30, 15]) {
  const lam = (360 / dphi) * P * Math.SQRT1_2;
  BREF.v = 1.5;
  const a = tableRun(lattice(dphi, { kdir: 'diag' }), Math.PI / 4);
  const b = tableRun(lattice(dphi, { kdir: 'diag', shuffle: true }), Math.PI / 4);
  const fmt = r => r.fell ? `転倒 t=${r.fell.toFixed(2)}s`.padEnd(24)
    : `${r.speed.toFixed(0)}mm/s ${r.dir.toFixed(0)}° 無接地${r.drop.toFixed(0)}% 傾き${r.tilt.toFixed(1)}°`.padEnd(24);
  console.log(`    ${String(dphi).padStart(3)}°  ${lam.toFixed(0).padStart(5)} mm   ${fmt(a)}  ${fmt(b)}`);
}
TREF.v = 3.0; BREF.v = 6;
console.log('    → 荷の差し渡しは 150mm。上下振幅 b=1.5mm（跳ねを抑えた値）。');

/* ── 追加: 上下振幅 b を下げると跳ねが収まるか（b と a は独立に選べる） ── */
console.log('\n■ 上下振幅 b を下げる（a=6mm 固定, Δφ=15°, k=斜め, ψ=45°, 車両）');
console.log('    b      速度      滑り   最小接地  荷重リップル');
for (const b of [6, 3, 1.5, 0.8]) {
  BREF.v = b;
  const r = polar(lattice(15, { kdir: 'diag' }), () => Math.PI / 4);
  console.log(`    ${String(b).padStart(3)}mm  ${r.speed.toFixed(0).padStart(4)} mm/s  ${(r.slip * 100).toFixed(0).padStart(3)}%  ` +
    `${String(r.minC).padStart(4)}    ${(r.nripple * 100).toFixed(0).padStart(3)}%`);
}
BREF.v = 6;
console.log('    → b は接地クリアランスだけ確保すればよく、推力は a·q̇ で決まる（傾斜円盤と同じ分離）。');

console.log(`
■ この試験が言えること

  1. 操舵は構造的に厳密。Δφ=60°・k斜めで8方位すべて誤差 0.00°、速度 184 mm/s（a·q̇ の 98%）、滑り 0%。
     接地している接触子の摩擦力が全部 ψ に共線なので、ψ に直交する速度成分は必ず減衰する。
     案Aの過拘束（X が (r·ωx,0)、Y が (0,r·ωy) を要求して競合する）が構造的に存在しない。
     停止軸もないので、案Aで必須だったリトラクトが要らない。

  2. 波数ベクトル k は軸に平行にしてはいけない。同位相の接触子が1列に並んで線接触になり、
     機体がその線まわりに首を振る。k を斜めにするだけで速度 50→184 mm/s、誤差 6.8→0.0°。

  3. 上下振幅 b と面内振幅 a は独立に選べる。推力は a·q̇ で決まるので、b は接地クリアランス分だけでよい。
     b を 6→0.8mm に下げると荷重リップルが 127%→6% に落ちて速度は変わらない。跳ねはこれで消える。

  4. 空間波が要るかどうかは荷の寸法で決まる。ただし傾斜円盤とは逆向き。
     ・車両（荷 ≫ λ）: シャッフルしても 185 対 184 mm/s。波は不要、位相差は当番表。
     ・卓上（荷が有限）: **λ ≲ 荷** が要る。荷の下で高さに差がないと、
       荷が下降中の接触子（面内速度が逆向き）にも乗り続けて推力が相殺される。
       λ=849mm で転倒、λ=424mm で 4mm/s、λ=106mm で 152mm/s。
     傾斜円盤は λ ≫ 荷 を要求した。こちらは λ ≲ 荷。設計指針が正反対になる。

■ 言えないこと

  ・実機がこう動くとは言えない。μ・接触剛性・接地窓はすべて未測定の仮定値。
  ・差動機構のバックラッシュ、連続スキンの面内せん断、接触子間干渉はモデルに入っていない。
    偏極式は面内振幅 a=6mm ぶんスキンをせん断するので、案Aより皮膜の負担は大きい。
  ・2入力では並進2自由度まで。ヨーを含む平面3自由度には一般化力が3系統要る。`);
