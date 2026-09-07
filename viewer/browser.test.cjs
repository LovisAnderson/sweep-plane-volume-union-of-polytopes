const { chromium } = require(process.env.PLAYWRIGHT_MODULE || "playwright");
const assert = require("node:assert/strict");
(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({
    viewport: { width: 1400, height: 1000 },
  });
  let calls = 0;
  const errors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  page.on("request", (r) => {
    if (r.url().endsWith("/api/compute")) calls++;
  });
  await page.goto(process.env.VIEWER_URL || "http://127.0.0.1:8765");
  await page
    .locator("#files")
    .setInputFiles(
      require("node:path").resolve(__dirname, "..") + "/inputs/square.ine",
    );
  await page.locator("#initial").fill("1");
  await page.locator("#compute").click();
  await page.waitForFunction(() => !document.getElementById("slider").disabled);
  assert.equal(await page.locator("#minimum").textContent(), "0");
  assert.equal(await page.locator("#maximum").textContent(), "3");
  assert.equal(await page.locator("#area").textContent(), "0.25 / 1");
  assert.equal(calls, 1);
  const before = await page.locator("#geometry line").getAttribute("y1");
  await page.locator("#slider").evaluate((e) => {
    e.value = 1000;
    e.dispatchEvent(new Event("input"));
  });
  assert.equal(await page.locator("#area").textContent(), "1 / 1");
  assert.equal(await page.locator("#lambda").inputValue(), "3");
  assert.notEqual(
    await page.locator("#geometry line").getAttribute("y1"),
    before,
  );
  await page.locator("#lambda").fill("1/2");
  await page.locator("#lambda").press("Tab");
  assert.equal(await page.locator("#area").textContent(), "0.0625 / 1");
  assert.equal(calls, 1);
  await page.locator("#lambda").fill("-9");
  await page.locator("#lambda").press("Tab");
  assert.equal(await page.locator("#lambda").inputValue(), "0");
  assert.equal(calls, 1);
  await page
    .locator("#files")
    .setInputFiles(
      require("node:path").resolve(__dirname, "..") + "/inputs/union.json",
    );
  assert.equal(await page.locator("#slider").isDisabled(), true);
  await page.locator("#compute").click();
  await page.waitForFunction(() => !document.getElementById("slider").disabled);
  await page.screenshot({ path: "/tmp/nefvol-viewer.png", fullPage: true });
  await page.locator("#direction").fill("-1,0");
  await page.locator("#compute").click();
  await page.waitForFunction(() => !document.getElementById("slider").disabled);
  assert.equal(await page.locator("#maximum").textContent(), "0");
  await page
    .locator("#files")
    .setInputFiles(
      require("node:path").resolve(__dirname, "..") + "/inputs/lshape3d.ine",
    );
  await page.locator("#compute").click();
  await page.waitForFunction(
    () => document.getElementById("status").className === "error",
  );
  assert.match(await page.locator("#status").textContent(), /only 2D/);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.screenshot({
    path: "/tmp/nefvol-viewer-mobile.png",
    fullPage: true,
  });
  assert.deepEqual(errors, []);
  console.log(
    "Browser checks passed: loading, initial lambda, endpoints, fractions, markers, cached calls, direction changes, 3D rejection, no console errors",
  );
  await browser.close();
})().catch((e) => {
  console.error(e);
  process.exit(1);
});
