// Chromium の場所を決める。
//
//   1. CHROMIUM_PATH が指定されていればそれ
//   2. playwright が自分で入れたもの（npm install && npx playwright install）
//   3. それも無ければ PLAYWRIGHT_BROWSERS_PATH 以下を探す
//
// 3 が要るのは、環境にプリインストールされた Chromium のビルド番号が
// playwright のピン留めと食い違うことがあるため。

import { chromium } from 'playwright';
import { existsSync, readdirSync } from 'fs';
import { join } from 'path';

function scan() {
  const root = process.env.PLAYWRIGHT_BROWSERS_PATH || '/opt/pw-browsers';
  if (!existsSync(root)) return null;
  const dirs = readdirSync(root).filter(d => d.startsWith('chromium')).sort().reverse();
  for (const d of dirs) {
    for (const rel of ['chrome-linux/chrome', 'chrome-linux/headless_shell',
                       'chrome-mac/Chromium.app/Contents/MacOS/Chromium']) {
      const p = join(root, d, rel);
      if (existsSync(p)) return p;
    }
  }
  return null;
}

export async function launch(opts = {}) {
  if (process.env.CHROMIUM_PATH) return chromium.launch({ ...opts, executablePath: process.env.CHROMIUM_PATH });
  try {
    return await chromium.launch(opts);
  } catch (e) {
    const p = scan();
    if (!p) throw e;
    return chromium.launch({ ...opts, executablePath: p });
  }
}
