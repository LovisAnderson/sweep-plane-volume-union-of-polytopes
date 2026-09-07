import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
const source = await readFile(new URL("./math.js", import.meta.url), "utf8");
const { rational, exact, prepare, evaluate, interpolate } = await import(
  `data:text/javascript;base64,${Buffer.from(source).toString("base64")}`
);
assert.equal(exact(rational("-1.25e2")), "-125");
assert.equal(exact(rational("1/-2")), "-1/2");
assert.throws(() => rational("1/0"));
assert.throws(() => rational("garbage"));
assert.equal(exact(interpolate(rational("-2"), rational("4"), 500)), "1");
const pieces = prepare({
  pieces: [
    { lo: null, coefficients: ["0"] },
    {
      lo: "10000000000000000",
      coefficients: [
        "100000000000000000000000000000000",
        "-20000000000000000",
        "1",
      ],
    },
  ],
});
assert.equal(exact(evaluate(pieces, rational("10000000000000001"))), "1");
console.log("Rational evaluation and interpolation tests passed");
