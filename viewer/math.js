// Exact rational evaluation avoids cancellation in translated polynomial pieces.
export function rational(s) {
  s = String(s).trim();
  if (/^[+-]?\d+\/[+-]?\d+$/.test(s)) {
    const [n, d] = s.split("/").map(BigInt);
    return norm(n, d);
  }
  const m = s.match(/^([+-]?)(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$/);
  if (!m || !(m[2] || m[3])) throw Error("Enter a number or fraction for λ");
  const exponent = Number(m[4] || 0) - (m[3] || "").length;
  if (Math.abs(exponent) > 1000 || s.length > 2000)
    throw Error("Number is too large");
  const n = BigInt((m[1] === "-" ? "-" : "") + (m[2] || "0") + (m[3] || ""));
  return exponent >= 0
    ? norm(n * 10n ** BigInt(exponent), 1n)
    : norm(n, 10n ** BigInt(-exponent));
}
function norm(n, d) {
  if (!d) throw Error("Zero denominator");
  if (d < 0n) {
    n = -n;
    d = -d;
  }
  let a = n < 0n ? -n : n,
    b = d;
  while (b) {
    [a, b] = [b, a % b];
  }
  return [n / a, d / a];
}
export const add = (a, b) => norm(a[0] * b[1] + b[0] * a[1], a[1] * b[1]);
export const mul = (a, b) => norm(a[0] * b[0], a[1] * b[1]);
export const compare = (a, b) =>
  a[0] * b[1] < b[0] * a[1] ? -1 : a[0] * b[1] > b[0] * a[1] ? 1 : 0;
export const number = (a) => Number(a[0]) / Number(a[1]);
export const exact = (a) => (a[1] === 1n ? String(a[0]) : `${a[0]}/${a[1]}`);
export function prepare(data) {
  return data.pieces.map((p) => ({
    lo: p.lo === null ? null : rational(p.lo),
    coefficients: p.coefficients.map(rational),
  }));
}
export function evaluate(pieces, x) {
  let piece = pieces[0];
  for (const p of pieces) {
    if (p.lo !== null && compare(x, p.lo) < 0) break;
    piece = p;
  }
  return piece.coefficients.reduceRight((v, c) => add(mul(v, x), c), [0n, 1n]);
}
export function interpolate(lo, hi, t) {
  return add(mul(lo, [1000n - BigInt(t), 1000n]), mul(hi, [BigInt(t), 1000n]));
}
