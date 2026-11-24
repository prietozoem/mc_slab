// Run with:  node test/test_engine.js
global.window = global;
require('../js/data.js');
const MC = require('../js/engine.js');
let fails = 0;
function check(name, ok, msg) { console.log((ok ? 'PASS ' : 'FAIL ') + name + (msg ? '  ' + msg : '')); if (!ok) fails++; }
const mat = k => new MC.Material(k, window.XCOM[k]);

// 1. table values / interpolation (water 1 MeV, XCOM: total with coherent = 0.07072 cm2/g)
{
  const m = mat('water').massXS(1);
  const tot = m.reduce((a, b) => a + b, 0);
  check('water 1 MeV total mass XS', Math.abs(tot - 0.07072) < 1e-5, tot.toFixed(5));
  const p = mat('lead').massXS(0.0885);       // just above the K edge (88 keV) -> post-edge photoelectric
  check('lead K-edge jump', p[2] > 7 && mat('lead').massXS(0.0879)[2] < 2, p[2].toFixed(3));
  const q = mat('lead').massXS(1.1);          // pair threshold region, interpolated
  check('lead pair above threshold', q[3] > 0 && q[3] < 0.0004, q[3].toExponential(3));
}

// 2. uncollided transmission against exp(-mu x)
for (const [k, E, T] of [['lead', 1.0, 2.0], ['water', 0.5, 20], ['iron', 0.2, 3], ['concrete', 2.0, 15]]) {
  const sim = MC.makeSim({ material: mat(k), E0: E, T: T, seed: 12345 });
  const N = 200000;
  for (let i = 0; i < N; i++) MC.runHistory(sim);
  const s = MC.summary(sim);
  const p = s.analyticUncollided, sd = Math.sqrt(p * (1 - p) / N);
  const z = (s.uncollided.p - p) / sd;
  check(`uncollided ${k} ${E} MeV ${T} cm`, Math.abs(z) < 4, `MC=${s.uncollided.p.toFixed(5)} analytic=${p.toFixed(5)} z=${z.toFixed(2)}`);
  const t = sim.tally;
  check(`  conservation ${k}`, t.transmitted + t.reflected + t.absorbed === t.N);
}

// 3. Klein-Nishina sampling: mean E'/E and mean cos(theta) against numerical integration
function knMoments(E) {
  const a = E / MC.ME; let n = 0, d = 0, nm = 0;
  const M = 200000;
  for (let i = 0; i < M; i++) {
    const mu = -1 + 2 * (i + 0.5) / M;
    const r = 1 / (1 + a * (1 - mu));
    const f = r * r * (r + 1 / r - (1 - mu * mu));     // dsigma/dmu up to a constant
    d += f; n += f * r; nm += f * mu;
  }
  return { Er: n / d, mu: nm / d };
}
{
  const rng = MC.makeRNG(7);
  for (const E of [0.01, 0.1, 0.662, 2, 10]) {
    let se = 0, sm = 0, n = 400000;
    for (let i = 0; i < n; i++) { const c = MC.sampleCompton(E, rng); se += c.Eout / E; sm += c.mu; }
    const ex = knMoments(E);
    check(`Compton KN E=${E}`, Math.abs(se / n - ex.Er) < 0.003 && Math.abs(sm / n - ex.mu) < 0.004,
      `<E'/E> ${(se / n).toFixed(4)} vs ${ex.Er.toFixed(4)};  <mu> ${(sm / n).toFixed(4)} vs ${ex.mu.toFixed(4)}`);
  }
}

// 4. Rayleigh sampler: forward peaked, valid range
{
  const rng = MC.makeRNG(3);
  let s = 0, n = 100000;
  for (let i = 0; i < n; i++) { const r = MC.sampleRayleigh(0.1, 82, rng); if (r.mu < -1 || r.mu > 1) throw new Error('bad mu'); s += r.mu; }
  check('Rayleigh forward peaked (Pb 100 keV)', s / n > 0.8, '<mu>=' + (s / n).toFixed(3));
  s = 0; for (let i = 0; i < n; i++) s += MC.sampleRayleigh(0.001, 8, rng).mu;
  check('Rayleigh -> Thomson at low E (<mu>~0)', Math.abs(s / n) < 0.01, '<mu>=' + (s / n).toFixed(4));
}

// 5. generator and fast loop give identical tallies for the same seed
{
  const a = MC.makeSim({ material: mat('iron'), E0: 0.5, T: 4, seed: 99 });
  const b = MC.makeSim({ material: mat('iron'), E0: 0.5, T: 4, seed: 99 });
  for (let i = 0; i < 3000; i++) { MC.runHistory(a); for (const ev of MC.historySteps(b)) { } }
  const same = a.tally.transmitted === b.tally.transmitted && a.tally.reflected === b.tally.reflected &&
    a.tally.absorbed === b.tally.absorbed && a.tally.byReaction.join() === b.tally.byReaction.join();
  check('generator == fast loop', same, a.tally.transmitted + ' ' + b.tally.transmitted);
}

// 6. energy range edge cases
{
  const sim = MC.makeSim({ material: mat('lead'), E0: 100, T: 5, seed: 1 });
  for (let i = 0; i < 2000; i++) MC.runHistory(sim);
  check('100 MeV runs', sim.tally.N === 2000);
  const sim2 = MC.makeSim({ material: mat('water'), E0: 0.001, T: 0.01, seed: 1 });
  for (let i = 0; i < 2000; i++) MC.runHistory(sim2);
  check('1 keV runs', sim2.tally.N === 2000);
}
// 7. line source: uniform in z over [-L/2, L/2], same transmission as the point source, generator == fast loop
{
  const L = 6, a = MC.makeSim({ material: mat('lead'), E0: 1, T: 2, seed: 5 }); MC.configureSource(a, L, 'z', 0);
  const b = MC.makeSim({ material: mat('lead'), E0: 1, T: 2, seed: 5 }); MC.configureSource(b, L, 'z', 0);
  let zmin = 1e9, zmax = -1e9, zs = 0, n = 40000;
  a.record = true;
  for (let i = 0; i < n; i++) { const r = MC.runHistory(a); const z0 = r.ph.path[0][3]; zmin = Math.min(zmin, z0); zmax = Math.max(zmax, z0); zs += z0; }
  check('line source z range', zmin >= -L / 2 && zmax <= L / 2 && zmax - zmin > 0.99 * L && Math.abs(zs / n) < 0.05, zmin.toFixed(3) + ' .. ' + zmax.toFixed(3));
  const p = Math.exp(-mat('lead').macro(1).tot * 2), s = MC.summary(a);
  check('line source uncollided = exp(-mu T)', Math.abs(s.uncollided.p - p) < 4 * Math.sqrt(p * (1 - p) / n), s.uncollided.p.toFixed(4) + ' vs ' + p.toFixed(4));
  for (let i = 0; i < 3000; i++) for (const ev of MC.historySteps(b)) { }
  const d = MC.makeSim({ material: mat('lead'), E0: 1, T: 2, seed: 5 }); MC.configureSource(d, L, 'z', 0);
  for (let i = 0; i < 3000; i++) MC.runHistory(d);
  check('line source: generator == fast loop', d.tally.transmitted === b.tally.transmitted && d.tally.reflected === b.tally.reflected && d.tally.byReaction.join() === b.tally.byReaction.join(), d.tally.transmitted + ' ' + b.tally.transmitted);
}
// 8. line source along x (depth range dX inside the slab): uniform depth, analytic uncollided, generator == fast loop
{
  const m = mat('iron'), T = 4, n = 100000;
  const a = MC.makeSim({ material: m, E0: 0.5, T: T, seed: 21 }); MC.configureSource(a, 3, 'x', 0.5);
  check('x-axis source clamped inside slab', a.srcLen === 3 && a.srcX0 === 0.5);
  const a2 = MC.makeSim({ material: m, E0: 0.5, T: T, seed: 21 }); MC.configureSource(a2, 10, 'x', 2.5);
  check('x-axis length clamped to T - x0', a2.srcLen === 1.5, String(a2.srcLen));
  let xmin = 1e9, xmax = -1e9, xs = 0; a.record = true;
  for (let i = 0; i < n; i++) { const r = MC.runHistory(a); const x0 = r.ph.path[0][0]; xmin = Math.min(xmin, x0); xmax = Math.max(xmax, x0); xs += x0; }
  check('x-axis source depth uniform in [x0, x0+dX]', xmin >= 0.5 && xmax <= 3.5 && xmax - xmin > 0.99 * 3 && Math.abs(xs / n - 2.0) < 0.02, xmin.toFixed(3) + ' .. ' + xmax.toFixed(3) + ' mean ' + (xs / n).toFixed(3));
  const s = MC.summary(a), p = s.analyticUncollided;
  check('x-axis source: uncollided = analytic average', Math.abs(s.uncollided.p - p) < 4 * Math.sqrt(p * (1 - p) / n), s.uncollided.p.toFixed(4) + ' vs ' + p.toFixed(4));
  const b = MC.makeSim({ material: m, E0: 0.5, T: T, seed: 21 }); MC.configureSource(b, 3, 'x', 0.5);
  const d = MC.makeSim({ material: m, E0: 0.5, T: T, seed: 21 }); MC.configureSource(d, 3, 'x', 0.5);
  for (let i = 0; i < 3000; i++) { MC.runHistory(d); for (const ev of MC.historySteps(b)) { } }
  check('x-axis source: generator == fast loop', d.tally.transmitted === b.tally.transmitted && d.tally.reflected === b.tally.reflected, d.tally.transmitted + ' ' + b.tally.transmitted);
  // y-axis source
  const c = MC.makeSim({ material: m, E0: 0.5, T: T, seed: 3 }); MC.configureSource(c, 5, 'y', 1); c.record = true;
  let ymax = 0; for (let i = 0; i < 20000; i++) { const r = MC.runHistory(c); ymax = Math.max(ymax, Math.abs(r.ph.path[0][1])); if (r.ph.path[0][0] !== 1) throw new Error('x0'); }
  check('y-axis source range', ymax <= 2.5 && ymax > 2.45, ymax.toFixed(3));
}
console.log(fails ? fails + ' FAILED' : 'ALL PASSED');
process.exit(fails ? 1 : 0);
