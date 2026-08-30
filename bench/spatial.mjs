// 「これは波動機構か」を、足の空間座標を持つモデルで判定し直す。
//
// bench/wavetest.mjs は無効だった。あのモデルは足の位置を持たず、機体を質点として
// 扱っていたので、各時刻の力は「位相の集合」だけの対称関数になる。位相をどう並べ替えても
// 結果が一致するのは数式上の恒等式であって、測定ではない。
//
// ここでは足に (x, z) を与え、剛体を6自由度（並進3・ヨー・ピッチ・ロール）で積分する。
// 位置が入ると、法線力が姿勢に、姿勢が接地分布に、接地分布がモーメントに戻る閉ループができる。
// 位相の空間配置はこのループを通してのみ効く。効くか効かないかを、ここで測る。
//
//   ケースA 車両: 足は機体に付き、地面を蹴る
//   ケースB 卓上: 足は地面に固定され、有限寸法の荷がその上を渡っていく（局所支持）

const N = 12;            // ラティス N×N
const P = 50;            // 足ピッチ [mm]
const CRANK = 6;         // クランク半径 [mm]（上下ストローク = 2·CRANK）
const DISC_E = 25;       // 傾斜円盤の張り出し [mm]（掃引速度 = DISC_E·ω）
const RPM = 300;
const W = RPM * Math.PI / 30;
const DELTA = Math.PI / 4;   // 掃引方向を決める位相差 δ
const MU = 0.6, G = 9810, VEPS = 6;
const ACC = 1000;        // N, kg → mm/s²

const DT0 = 2e-5, T = 3.0, WARM = 1.0;
const DTREF = { v: DT0 };          // 収束確認で差し替える
const TIPPED = 8 * Math.PI / 180;   // これを超えたら小角近似が破れる → 転倒扱い
const TAU = 2 * Math.PI;

/* 足 j の、機体基準での先端高さと先端接線速度。
   h = -CRANK·cos(θ) は最下点で -CRANK。χ は最下点の方位。 */
function tip(theta) {
  const chi = theta + DELTA + Math.PI / 2;
  return [-CRANK * Math.cos(theta), W * DISC_E * Math.sin(chi), -W * DISC_E * Math.cos(chi)];
}

/* ── ケースA: 車両 ──────────────────────────────────────────
   足は機体に固定。機体はヨーで回るので足の世界座標も回る。 */
function vehicle(phase) {
  const DT = DTREF.v;
  const M = 20, L = N * P, HCG = 100;
  const Ih = M * (L * L / 12 + HCG * HCG);   // ピッチ・ロール慣性 [kg·mm²]
  const Iy = M * (L * L / 6);                // ヨー慣性
  const K = 47, C = 0.15;                    // 接地ばね [N/mm], 減衰 [N·s/mm]

  const ax = [], az = [];
  for (let i = 0; i < N; i += 1) for (let k = 0; k < N; k += 1) {
    ax.push((i - (N - 1) / 2) * P); az.push((k - (N - 1) / 2) * P);
  }
  const nf = ax.length;

  let X = 0, Z = 0, Y = CRANK + 0.4, psi = 0, tx = 0, tz = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  const prevPen = new Float64Array(nf);

  let acc = { dx: 0, dz: 0, slip: 0, wt: 0, ay2: 0, tilt2: 0, psi0: 0, minC: 1e9, drop: 0 };

  for (let s = 0; s < T / DT; s += 1) {
    const t = s * DT, live = t > WARM;
    if (live && acc.wt === 0) { acc.psi0 = psi; }
    const cp = Math.cos(psi), sp = Math.sin(psi);

    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0, sl = 0, wsum = 0;
    for (let j = 0; j < nf; j += 1) {
      const [h, ux, uz] = tip(phase[j] - W * t);
      const rx = ax[j] * cp - az[j] * sp, rz = ax[j] * sp + az[j] * cp;
      const y = Y - tx * rz + tz * rx + h;
      const pen = -y;
      if (pen <= 0) { prevPen[j] = 0; continue; }
      const rate = (pen - prevPen[j]) / DT;
      prevPen[j] = pen;
      const Nf = Math.max(0, K * pen + C * rate);
      if (Nf <= 0) continue;
      cnt += 1; SN += Nf;
      // 先端接線速度を世界系へ（ヨーで回る）
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
      acc.dx += vX * DT; acc.dz += vZ * DT; acc.wt += DT;
      acc.slip += (wsum > 0 ? sl / wsum : 0) * DT;
      const ay = ACC * SN / M - G;
      acc.ay2 += ay * ay * DT;
      acc.tilt2 += (tx * tx + tz * tz) * DT;
      if (cnt < acc.minC) acc.minC = cnt;
      if (cnt === 0) acc.drop += DT;
    }
  }
  const wt = acc.wt;
  return {
    speed: Math.hypot(acc.dx, acc.dz) / wt,
    dir: Math.atan2(acc.dz, acc.dx) * 180 / Math.PI,
    slip: acc.slip / wt / (DISC_E * W),
    yaw: (psi - acc.psi0) / wt * 180 / Math.PI,
    vib: Math.sqrt(acc.ay2 / wt) / 9810,
    tilt: Math.sqrt(acc.tilt2 / wt) * 180 / Math.PI,
    minC: acc.minC, drop: acc.drop / wt * 100
  };
}

/* ⚠ ケースB の結論（有限の荷には空間波が要る）は、後の検証で成立しないことが
   分かっています。同じ欠陥が bench/polar.mjs の卓上試験にもありました:
     ・格子を mod N で引くので位相場に周期境界の不連続がある
     ・荷が接触面より低い位置から始まり、静的釣合いでも回転数ランプでもない
     ・位相原点を1点しか見ていない
   位相原点を8通り走査すると、シャッフルでも 2/8〜7/8 は完走し、逆に空間波が
   0/8 になる Δφ もあります（docs/polarization-steering.md 試験2ケースB）。
   ここの「シャッフルは転倒する」は**採用しないでください**。ケースA（車両）の
   結論は影響を受けません。

   ── ケースB: 卓上 ──────────────────────────────────────────
   足は地面に固定。有限寸法の荷がその上を渡る。荷の下に来る足だけが効く（局所支持）。
   ラティスは周期 N·P で無限に続くとみなし、添字を mod N で引く。 */
function table(phase) {
  const DT = DTREF.v;
  const M = 0.5, S = 150, HCG = 15;          // 荷: 0.5kg, 150mm角, 重心高さ15mm
  const Ih = M * (S * S / 12 + HCG * HCG);
  const Iy = M * (S * S / 6);
  const K = 2.2, C = 0.006;
  const SPAN = N * P;
  const at = (i, k) => phase[((i % N) + N) % N * N + ((k % N) + N) % N];

  let X = 0, Z = 0, Y = CRANK - 0.15, psi = 0, tx = 0, tz = 0;
  let vX = 0, vZ = 0, vY = 0, vp = 0, wx = 0, wz = 0;
  const prevPen = new Map();

  let acc = { dx: 0, dz: 0, slip: 0, wt: 0, ay2: 0, tilt2: 0, psi0: 0, minC: 1e9, drop: 0 };
  let fell = 0;

  for (let s = 0; s < T / DT; s += 1) {
    const t = s * DT, live = t > WARM;
    // 傾きが TIPPED を超えたら小角近似も「荷が載っている」前提も破れる。転倒として打ち切る。
    if (Math.hypot(tx, tz) > TIPPED) { fell = t; break; }
    if (live && acc.wt === 0) acc.psi0 = psi;
    const cp = Math.cos(psi), sp = Math.sin(psi);
    const i0 = Math.floor((X - S / 2) / P), i1 = Math.ceil((X + S / 2) / P);
    const k0 = Math.floor((Z - S / 2) / P), k1 = Math.ceil((Z + S / 2) / P);

    let SN = 0, Fx = 0, Fz = 0, Ty = 0, Tx = 0, Tz = 0, cnt = 0, sl = 0, wsum = 0;
    for (let i = i0; i <= i1; i += 1) for (let k = k0; k <= k1; k += 1) {
      const gx0 = i * P, gz0 = k * P;
      const rx = gx0 - X, rz = gz0 - Z;                     // 荷の中心から見た接地点
      const bx = rx * cp + rz * sp, bz = -rx * sp + rz * cp; // 荷の機体座標
      if (Math.abs(bx) > S / 2 || Math.abs(bz) > S / 2) continue;
      const key = i * 100000 + k;
      const [h, ux, uz] = tip(at(i, k) - W * t);
      const top = CRANK + h;                                // 足先の世界高さ（最下点で0）
      const under = Y - tx * rz + tz * rx;                  // 荷の下面
      const pen = top - under;
      if (pen <= 0) { prevPen.set(key, 0); continue; }
      const rate = (pen - (prevPen.get(key) || 0)) / DT;
      prevPen.set(key, pen);
      const Nf = Math.max(0, K * pen + C * rate);
      if (Nf <= 0) continue;
      cnt += 1; SN += Nf;
      const gx = vX + vp * rz - ux, gz = vZ - vp * rx - uz;  // 荷 − 足
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
    if (Y < -2) { Y = -2; if (vY < 0) vY = 0; }             // 落ち切りは床で止める

    if (live) {
      acc.dx += vX * DT; acc.dz += vZ * DT; acc.wt += DT;
      acc.slip += (wsum > 0 ? sl / wsum : 0) * DT;
      const ay = ACC * SN / M - G;
      acc.ay2 += ay * ay * DT;
      acc.tilt2 += (tx * tx + tz * tz) * DT;
      if (cnt < acc.minC) acc.minC = cnt;
      if (cnt === 0) acc.drop += DT;
    }
    void SPAN;
  }
  if (fell) return { fell };
  const wt = acc.wt;
  return {
    speed: Math.hypot(acc.dx, acc.dz) / wt,
    dir: Math.atan2(acc.dz, acc.dx) * 180 / Math.PI,
    slip: acc.slip / wt / (DISC_E * W),
    yaw: (psi - acc.psi0) / wt * 180 / Math.PI,
    vib: Math.sqrt(acc.ay2 / wt) / 9810,
    tilt: Math.sqrt(acc.tilt2 / wt) * 180 / Math.PI,
    minC: acc.minC, drop: acc.drop / wt * 100
  };
}

/* ── 位相の配り方（同じ位相の集合を、違う場所に配る） ── */
const idx = (i, k) => i * N + k;
const build = f => { const a = new Float64Array(N * N); for (let i = 0; i < N; i += 1) for (let k = 0; k < N; k += 1) a[idx(i, k)] = f(i, k); return a; };

let seed = 20240131;
const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;

const waveX    = build(i => i * TAU / N);                       // x方向へ進む波（1波長）
const waveDiag = build((i, k) => (i + k) * TAU / N);            // 斜め45°へ進む波
const wave2    = build(i => i * 2 * TAU / N);                   // x方向 2波長
const group3   = build((i, k) => ((i + k) % 3) * TAU / 3);      // 3相の市松
const same     = build(() => 0);                                // 全部同位相

// waveX と同じ位相の集合を、空間だけシャッフルする（波としての勾配だけを壊す）
const shuffled = waveX.slice();
for (let i = shuffled.length - 1; i > 0; i -= 1) {
  const j = Math.floor(rnd() * (i + 1));
  const t = shuffled[i]; shuffled[i] = shuffled[j]; shuffled[j] = t;
}
const random = build(() => rnd() * TAU);

const CASES = [
  ['波・x方向 1波長',   waveX],
  ['波・斜め45° 1波長', waveDiag],
  ['波・x方向 2波長',   wave2],
  ['同じ位相集合を空間シャッフル', shuffled],
  ['完全ランダム',      random],
  ['3相の市松',         group3],
  ['全部同位相',        same],
];

function table1(title, fn, note) {
  console.log(title);
  console.log('  配り方                          速度      方向     滑り   ヨー      上下振動  傾き    最小接地  無接地');
  for (const [label, ph] of CASES) {
    const r = fn(ph);
    if (r.fell) { console.log(`  ${label.padEnd(30)} ── 荷が転倒（t=${r.fell.toFixed(2)}s、傾き8°超で打ち切り）──`); continue; }
    console.log(
      `  ${label.padEnd(30)} ${r.speed.toFixed(1).padStart(6)}  ${r.dir.toFixed(1).padStart(7)}°  ` +
      `${(r.slip * 100).toFixed(0).padStart(3)}%  ${r.yaw.toFixed(1).padStart(7)}°/s  ` +
      `${r.vib.toFixed(2).padStart(5)} G  ${r.tilt.toFixed(3).padStart(6)}°  ` +
      `${String(r.minC).padStart(4)}  ${r.drop.toFixed(1).padStart(5)}%`);
  }
  console.log(note + '\n');
}

/* 刻み幅を変えて、どの量が収束していてどの量が収束していないかを出す。
   結論に使ってよいのは収束している量だけ。 */
function convergence() {
  console.log('■ 収束確認（dt を 4 倍まで細かくして値が動くか見る）');
  const probes = [['A 波・x方向1波長', vehicle, waveX], ['A 空間シャッフル', vehicle, shuffled],
                  ['B 波・x方向1波長', table, waveX],   ['B 空間シャッフル', table, shuffled]];
  for (const [label, fn, ph] of probes) {
    const line = [];
    for (const dt of [2e-5, 1e-5, 5e-6]) {
      DTREF.v = dt;
      const r = fn(ph);
      line.push(r.fell ? `転倒 t=${r.fell.toFixed(2)}s`
                       : `${r.speed.toFixed(0)}mm/s 傾き${r.tilt.toFixed(2)}° 振動${r.vib.toFixed(2)}G 無接地${r.drop.toFixed(1)}%`);
    }
    DTREF.v = DT0;
    console.log(`  ${label.padEnd(18)} dt=2e-5: ${line[0]}`);
    console.log(`  ${''.padEnd(18)} dt=1e-5: ${line[1]}`);
    console.log(`  ${''.padEnd(18)} dt=5e-6: ${line[2]}`);
  }
  console.log('');
}

console.log(`■ 空間座標つき6自由度モデル（${N}×${N}足, ピッチ${P}mm, クランク${CRANK}mm, 円盤張出し${DISC_E}mm, ${RPM}rpm, δ=45°, μ=${MU})`);
console.log(`  積分 dt=${DT0}s, 計測 ${WARM}〜${T}s。速度は mm/s。\n`);
table1('■ ケースA 車両（足が機体に付き、地面を蹴る。全足接地しうる）', vehicle,
  '  ※ ヨーは自由。方向は世界系なのでヨーが大きい行では意味が薄れる。');
convergence();
table1('■ ケースB 卓上（足は地面に固定、150mm角・0.5kgの荷が渡る。局所支持）', table,
  '  ※ 荷の下に入る足は 3×3〜4×4 程度。無接地率は「支持が完全に切れた時間の割合」。\n' +
  '  ※ ケースBの上下振動は衝突が支配していて刻み幅に収束しない。結論には使わないこと。');

console.log(`■ この試験が言えること / 言えないこと

  言えること
   ・ケースA（大きな剛体の下に144足が全部入る）では、支持が切れない配り方どうしなら
     速度 776〜782 mm/s・方向 ±3° で一致する。空間シャッフルは斜め波と同等以下の
     振動（0.19G 対 0.18G）で走る。ここでは位相差は当番表であって波ではない。
     しかも「x方向1波長」は同位相の足が1列に並ぶので線接触になり、機体がその線まわりに
     首を振る（傾き1.44°、無接地10.8%、1.34G）。波らしい配り方がむしろ一番悪い。
   ・ケースB（荷が有限で、下に3×4足しか入らない）では逆で、空間勾配がないと荷は
     0.1〜0.22秒で転倒する。2波長でも 0.83秒で転倒する。荷の差し渡しに対して
     波長が十分長い勾配だけが荷を支え続ける。ここでは位相の空間配置が本質。
   ・したがって「波かどうか」は足の配り方ではなく荷の寸法で決まる。
     荷 ≫ 波長なら当番表に、荷 ≪ 波長なら波になる。車両の車体は前者の極限。

  言えないこと
   ・実機がこう動くとは言えない。μ、パッド剛性 K、接地窓はどれも未測定の仮定値。
     とくにケースBの結果は K（＝沈み込み量）で接地窓の幅が変わるので敏感。
   ・ケースBの絶対値（速度・滑り・ヨー）は刻み幅で ±10% 動く。
     収束しているのは転倒するかどうかの判定と転倒時刻（±0.01秒）まで。
   ・円盤の縁が線接触することによる応力・摩耗はこのモデルに入っていない。`);
