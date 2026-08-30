import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';

const executablePath = process.env.CHROMIUM_PATH || undefined;
const browser = await chromium.launch({ headless: true, ...(executablePath ? { executablePath } : {}) });

try {
  const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  const html = fileURLToPath(new URL('../sim-mapping.html', import.meta.url));
  await page.goto(pathToFileURL(html).href);
  await page.waitForFunction(() => {
    const d = window.mappingSimDiagnostics?.();
    return d?.modelReady === true && d.frameCount > 0;
  });

  const reference = await page.evaluate(() => window.mappingSimDiagnostics());
  assert.equal(reference.pinCount, 112);
  assert.equal(reference.m, 6);
  assert.equal(reference.finite, true);
  assert.ok(reference.gainP05 > 0.8 && reference.gainP05 < 1.2);
  assert.ok(reference.derived.speed > 20 && reference.derived.speed < 22);

  await page.locator('#playButton').click();
  await page.locator('#waveRatio').selectOption('2');
  await page.waitForFunction(() => {
    const d = window.mappingSimDiagnostics?.();
    return d?.modelReady === true && d.m === 2;
  });
  const undersampled = await page.evaluate(() => window.mappingSimDiagnostics());
  assert.ok(undersampled.gainP05 < 0.2);
  assert.match(await page.locator('#statusTitle').textContent(), /条件外/);

  for (const view of ['top', 'section', 'iso']) await page.locator(`[data-view="${view}"]`).click();
  assert.deepEqual(errors, []);
  console.log(`mapping web smoke: pass; λ/p=6 gain05=${reference.gainP05.toFixed(3)}, ` +
              `λ/p=2 gain05=${undersampled.gainP05.toFixed(3)}`);
} finally {
  await browser.close();
}
