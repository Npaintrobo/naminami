// 離散ピン写像の上界モデルと、モデルに依存しにくい必要条件を合わせて寸法領域を走査する。
import { directionEnvelope, selfCheck } from './discrete-map.mjs';

const G = 9.80665;
const cache = new Map();

function mapEnvelope(N, m, reg) {
  const key = `${N}/${m}/${reg}`;
  if (!cache.has(key)) cache.set(key, directionEnvelope({ N, samplesPerWave: m, reg }));
  return cache.get(key);
}

export function evaluateDesign({ N, pitch, samplesPerWave, thickness, amplitude, rpm,
                                 regStress = 1, minWavesAcross = 2,
                                 slopeLimit = 0.25, strainLimit = 0.05,
                                 verticalAccelerationLimitG = 0.8,
                                 thicknessToWavelengthLimit = 0.1 } = {}) {
  const wavelength = pitch * samplesPerWave, diameter = N * pitch;
  const k = 2 * Math.PI / wavelength, omega = rpm * Math.PI / 30;
  const ideal = mapEnvelope(N, samplesPerWave, 0), compliant = mapEnvelope(N, samplesPerWave, regStress);
  const slope = k * amplitude;
  const surfaceStrain = thickness * amplitude * k * k / 2;
  const verticalAccelerationG = (amplitude / 1000) * omega * omega / G;
  const surfaceSpeed = compliant.gainMedianMin * (thickness / 2) * amplitude * k * omega;
  const wavesAcrossDiameter = diameter / wavelength;
  const thicknessToWavelength = thickness / wavelength;
  const mappingPass = ideal.gainP05Min >= 0.90 && ideal.leakageP95Max <= 0.10 &&
                      compliant.gainMedianMin >= 0.85 && compliant.verticalP05Min >= 0.85;
  const mechanicsPass = slope <= slopeLimit && surfaceStrain <= strainLimit &&
                        verticalAccelerationG <= verticalAccelerationLimitG &&
                        thicknessToWavelength <= thicknessToWavelengthLimit;
  const finiteDiskPass = wavesAcrossDiameter >= minWavesAcross;
  return {
    N, pitch, samplesPerWave, thickness, amplitude, rpm, regStress, wavelength, diameter,
    wavesAcrossDiameter, thicknessToWavelength, slope, surfaceStrain, verticalAccelerationG, surfaceSpeed,
    mappingPass, mechanicsPass, finiteDiskPass, pass: mappingPass && mechanicsPass && finiteDiskPass,
    mapIdeal: ideal, mapCompliant: compliant,
  };
}

if (import.meta.url === new URL(`file://${process.argv[1]}`).href) {
  selfCheck();
  const rows = [];
  for (const N of [8, 12, 16]) for (const pitch of [15, 20, 25, 30, 40, 50])
    for (const samplesPerWave of [3, 4, 6, 8]) for (const thickness of [4, 6, 8, 10, 12])
      for (const amplitude of [1, 2, 3, 4, 5, 6, 8]) for (const rpm of [100, 150, 200, 300, 400]) {
        const r = evaluateDesign({ N, pitch, samplesPerWave, thickness, amplitude, rpm });
        if (r.pass) rows.push(r);
      }

  rows.sort((a, b) => b.surfaceSpeed - a.surfaceSpeed || a.diameter - b.diameter);
  console.log('成立候補（無荷重写像上界＋必要条件、速度順）');
  console.log('N  p  λ/p  D  λ  t  W rpm  v_surface slope strain az/g');
  for (const r of rows.slice(0, 15)) console.log(
    `${String(r.N).padStart(2)} ${String(r.pitch).padStart(2)} ${String(r.samplesPerWave).padStart(4)} ` +
    `${String(r.diameter).padStart(3)} ${String(r.wavelength).padStart(3)} ${String(r.thickness).padStart(2)} ` +
    `${String(r.amplitude).padStart(2)} ${String(r.rpm).padStart(3)} ` +
    `${r.surfaceSpeed.toFixed(1).padStart(9)} ${r.slope.toFixed(3).padStart(5)} ` +
    `${(100 * r.surfaceStrain).toFixed(1).padStart(5)}% ${r.verticalAccelerationG.toFixed(2).padStart(4)}`);

  const current = evaluateDesign({ N: 12, pitch: 50, samplesPerWave: 2, thickness: 6, amplitude: 5, rpm: 300 });
  const prototype = evaluateDesign({ N: 12, pitch: 20, samplesPerWave: 6, thickness: 10,
    amplitude: 2, rpm: 400, slopeLimit: 0.2, strainLimit: 0.03, verticalAccelerationLimitG: 0.5 });
  const balanced = evaluateDesign({ N: 12, pitch: 25, samplesPerWave: 6, thickness: 10,
    amplitude: 4, rpm: 400 });
  console.log('\n比較');
  for (const [name, r] of [['現行例', current], ['最小構成の原理試作', prototype], ['20kg級の次段候補', balanced]]) console.log(
    `${name}: D=${r.diameter}mm λ=${r.wavelength}mm λ/p=${r.samplesPerWave}, ` +
    `mapping=${r.mappingPass ? 'pass' : 'FAIL'}, finite-disk=${r.finiteDiskPass ? 'pass' : 'FAIL'}, ` +
    `v≈${r.surfaceSpeed.toFixed(1)}mm/s, slope=${r.slope.toFixed(3)}, ` +
    `strain≈${(100 * r.surfaceStrain).toFixed(1)}%, az=${r.verticalAccelerationG.toFixed(2)}g`);
}
