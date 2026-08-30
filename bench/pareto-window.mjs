// 円環写像案の多目的設計領域探索。
// 写像上界に、板の曲げ力・積載圧・慣性・準静的性を加え、用途別の代表点を抽出する。
import { evaluateDesign } from './design-window.mjs';
import { selfCheck } from './discrete-map.mjs';

const G = 9.80665;

function pinCount(N) {
  let n = 0, R = N / 2;
  for (let j = 0; j < N; j += 1) for (let i = 0; i < N; i += 1) {
    const x = i - (N - 1) / 2, y = j - (N - 1) / 2;
    if (Math.hypot(x, y) <= R) n += 1;
  }
  return n;
}

export function addPhysicalEstimates(r, { E = 5, nu = 0.49, density = 1100,
                                           payloadMass = 20, grooveMu = 0.2,
                                           steeringRateDeg = 90 } = {}) {
  const k = 2 * Math.PI / r.wavelength, omega = r.rpm * Math.PI / 30;
  const bendingRigidity = E * r.thickness ** 3 / (12 * (1 - nu * nu)); // N mm
  const bendPressure = bendingRigidity * k ** 4 * r.amplitude;          // N/mm²
  const payloadPressure = payloadMass * G / (Math.PI * (r.diameter / 2) ** 2);
  const inertiaPressure = density * (r.thickness / 1000) * (r.amplitude / 1000) *
                          omega ** 2 / 1e6;                             // N/mm²
  const forcePerPin = (bendPressure + payloadPressure + inertiaPressure) * r.pitch ** 2;
  const sheetMass = density * Math.PI * (r.diameter / 2000) ** 2 * (r.thickness / 1000);
  const totalNormal = (payloadMass + sheetMass) * G;
  const steeringTorque = grooveMu * totalNormal * (r.diameter / 1000) / 3; // μ N 2R/3
  const steeringPower = steeringTorque * steeringRateDeg * Math.PI / 180;
  const shearModulusPa = E * 1e6 / (2 * (1 + nu));
  const shearWaveSpeed = Math.sqrt(shearModulusPa / density);
  const quasiStaticRatio = omega * (r.thickness / 1000) / shearWaveSpeed;
  return { ...r, E, nu, density, payloadMass, pinCount: pinCount(r.N), bendingRigidity,
    bendPressure, payloadPressure, inertiaPressure, forcePerPin, sheetMass,
    steeringTorque, steeringPower, shearWaveSpeed, quasiStaticRatio };
}

export function selfCheckPhysical() {
  const opts = { N: 12, pitch: 20, samplesPerWave: 6, thickness: 10, amplitude: 2, rpm: 400,
    regStress: 1, minWavesAcross: 2, slopeLimit: 0.2, strainLimit: 0.03,
    verticalAccelerationLimitG: 0.5 };
  const base = evaluateDesign(opts), e5 = addPhysicalEstimates(base, { E: 5, payloadMass: 5 });
  const e10 = addPhysicalEstimates(base, { E: 10, payloadMass: 5 });
  if (!base.pass) throw new Error('reference prototype unexpectedly rejected');
  if (e5.pinCount !== 112) throw new Error(`pin count regression: ${e5.pinCount}`);
  if (Math.abs(e10.bendPressure / e5.bendPressure - 2) > 1e-12)
    throw new Error('bending pressure must scale linearly with E');
  if (!(e5.forcePerPin > 3 && e5.forcePerPin < 5))
    throw new Error(`reference pin-force regression: ${e5.forcePerPin}`);
  const invalidThinPlate = evaluateDesign({ ...opts, pitch: 15, samplesPerWave: 8, thickness: 20 });
  if (invalidThinPlate.mechanicsPass) throw new Error('t/lambda validity gate failed');
  return { referenceSpeed: e5.surfaceSpeed, referenceForce: e5.forcePerPin };
}

function enumerate(scenario) {
  const out = [];
  for (const N of [8, 12, 16, 20]) for (const pitch of [12.5, 15, 20, 25, 30, 40, 50])
    for (const samplesPerWave of [4, 6, 8, 10]) for (const thickness of [4, 6, 8, 10, 12, 15, 20])
      for (const amplitude of [1, 2, 3, 4, 5, 6, 8, 10]) for (const rpm of [60, 100, 150, 200, 300, 400, 500]) {
        const base = evaluateDesign({ N, pitch, samplesPerWave, thickness, amplitude, rpm,
          regStress: scenario.regStress, minWavesAcross: scenario.minWavesAcross,
          slopeLimit: scenario.slopeLimit, strainLimit: scenario.strainLimit,
          verticalAccelerationLimitG: scenario.verticalAccelerationLimitG,
          thicknessToWavelengthLimit: scenario.thicknessToWavelengthLimit });
        if (!base.pass) continue;
        const r = addPhysicalEstimates(base, scenario.physical);
        if (r.forcePerPin > scenario.maxForcePerPin || r.quasiStaticRatio > 0.1) continue;
        out.push(r);
      }
  return out;
}

function bestWithin(rows, maxDiameter, minSpeed = 0) {
  const q = rows.filter(r => r.diameter <= maxDiameter && r.surfaceSpeed >= minSpeed);
  q.sort((a, b) => b.surfaceSpeed - a.surfaceSpeed || a.surfaceStrain - b.surfaceStrain ||
                    a.forcePerPin - b.forcePerPin);
  return q[0];
}

function bestForN(rows, N) {
  const q = rows.filter(r => r.N === N);
  q.sort((a, b) => b.surfaceSpeed - a.surfaceSpeed || a.surfaceStrain - b.surfaceStrain ||
                    a.forcePerPin - b.forcePerPin);
  return q[0];
}

function leastPins(rows, minSpeed) {
  const q = rows.filter(r => r.surfaceSpeed >= minSpeed);
  q.sort((a, b) => a.pinCount - b.pinCount || a.diameter - b.diameter ||
                    a.surfaceStrain - b.surfaceStrain || b.surfaceSpeed - a.surfaceSpeed);
  return q[0];
}

function lowestStress(rows, maxDiameter, minSpeed) {
  const q = rows.filter(r => r.diameter <= maxDiameter && r.surfaceSpeed >= minSpeed);
  q.sort((a, b) => a.surfaceStrain - b.surfaceStrain || a.verticalAccelerationG - b.verticalAccelerationG ||
                    a.forcePerPin - b.forcePerPin || a.diameter - b.diameter);
  return q[0];
}

function dominates(a, b) {
  const noWorse = a.diameter <= b.diameter && a.surfaceSpeed >= b.surfaceSpeed &&
                  a.surfaceStrain <= b.surfaceStrain && a.forcePerPin <= b.forcePerPin;
  const better = a.diameter < b.diameter || a.surfaceSpeed > b.surfaceSpeed ||
                 a.surfaceStrain < b.surfaceStrain || a.forcePerPin < b.forcePerPin;
  return noWorse && better;
}

function pareto(rows) {
  return rows.filter((r, i) => !rows.some((q, j) => i !== j && dominates(q, r)));
}

function show(label, r) {
  if (!r) return console.log(`${label}: 該当なし`);
  console.log(`${label}: N=${r.N}, p=${r.pitch}mm, λ/p=${r.samplesPerWave}, ` +
    `D=${r.diameter}mm, λ=${r.wavelength}mm, t=${r.thickness}mm, W=${r.amplitude}mm, ` +
    `${r.rpm}rpm | v上界=${r.surfaceSpeed.toFixed(1)}mm/s, strain=${(100*r.surfaceStrain).toFixed(1)}%, ` +
    `az=${r.verticalAccelerationG.toFixed(2)}g, Fpin≈${r.forcePerPin.toFixed(1)}N, ` +
    `Tsteer≈${r.steeringTorque.toFixed(1)}Nm`);
}

if (import.meta.url === new URL(`file://${process.argv[1]}`).href) {
  selfCheck();
  const physicalCheck = selfCheckPhysical();
  console.log(`self-check: reference v=${physicalCheck.referenceSpeed.toFixed(2)}mm/s, ` +
              `Fpin=${physicalCheck.referenceForce.toFixed(2)}N, scaling/validity gates=pass`);
  const scenarios = {
    prototype: { regStress: 1, minWavesAcross: 2, slopeLimit: 0.20, strainLimit: 0.03,
      verticalAccelerationLimitG: 0.5, thicknessToWavelengthLimit: 0.1, maxForcePerPin: 10,
      physical: { payloadMass: 5, grooveMu: 0.2, steeringRateDeg: 30 } },
    proof: { regStress: 3, minWavesAcross: 2.5, slopeLimit: 0.20, strainLimit: 0.03,
      verticalAccelerationLimitG: 0.5, thicknessToWavelengthLimit: 0.1, maxForcePerPin: 10,
      physical: { payloadMass: 5, grooveMu: 0.2, steeringRateDeg: 30 } },
    balanced: { regStress: 1, minWavesAcross: 2, slopeLimit: 0.25, strainLimit: 0.05,
      verticalAccelerationLimitG: 0.8, thicknessToWavelengthLimit: 0.1, maxForcePerPin: 15,
      physical: { payloadMass: 20, grooveMu: 0.2, steeringRateDeg: 90 } },
    preloaded: { regStress: 1, minWavesAcross: 2, slopeLimit: 0.35, strainLimit: 0.08,
      verticalAccelerationLimitG: 1.2, thicknessToWavelengthLimit: 0.1, maxForcePerPin: 25,
      physical: { payloadMass: 20, grooveMu: 0.2, steeringRateDeg: 90 } },
  };

  for (const [name, scenario] of Object.entries(scenarios)) {
    const rows = enumerate(scenario), front = pareto(rows);
    console.log(`\n■ ${name}: ${rows.length}候補、Pareto ${front.length}点`);
    show('D≤250 最速', bestWithin(rows, 250));
    show('D≤300 最速', bestWithin(rows, 300));
    show('D≤400 最速', bestWithin(rows, 400));
    show('D≤300・20mm/s以上で最低ひずみ', lowestStress(rows, 300, 20));
    show('12×12 最速', bestForN(rows, 12));
    show('16×16 最速', bestForN(rows, 16));
    show('20×20 最速', bestForN(rows, 20));
    show('20mm/s以上で最少ピン', leastPins(rows, 20));
  }
}
