/* Tests of the neutron part of the engine. Uses test/data_n_mini.js (tiny subset) unless js/data_n.js has all nuclides.
 * Run: node test/test_neutron.js [path/to/data_n.js] */
'use strict';
global.window = global;
const path = require('path');
require(path.resolve(process.argv[2] || path.join(__dirname, 'data_n_mini.js')));
const MC = require('../js/engine.js');
let fails = 0;
function check(name, ok, info) { console.log((ok ? 'PASS ' : 'FAIL ') + name + (info !== undefined ? '  [' + info + ']' : '')); if (!ok) fails++; }
const have = (k) => { try { return MC.makeNeutronMaterial(k, 1, {}, window.NXS); } catch (e) { return null; } };
const mat = (k, rho) => MC.makeNeutronMaterial(k, rho || (k === 'water' ? 1 : 2.699), {}, window.NXS);

// 1. materials and data
const w = mat('water'), al = mat('aluminum');
check('water: N_H = 2 N_O and N_O = rho NA / M', Math.abs(w.nuc[0].N / w.nuc[1].N - 2) < 1e-12 && Math.abs(w.nuc[1].N - 0.6022140857 / 18.01528) < 2e-5, w.nuc[1].N.toFixed(6));
let miss = null; try { MC.makeNeutronMaterial('lead', 11.35, {}, window.NXS); } catch (e) { miss = e.missing; }
check('missing nuclides are reported', Array.isArray(miss) && miss.indexOf('pb208') >= 0, miss && miss.join(','));

// 2. macroscopic cross sections against a direct sum
{
  const E = 1.0, m = w.macro(E);
  let t = 0;
  for (const nu of w.nuc) for (let r = 0; r < 4; r++) t += nu.N * MC.evalCurve(nu.c[r], E);
  check('macro total = sum of N sigma', Math.abs(m.tot - t) < 1e-12 * t, m.tot.toFixed(5));
  check('macro arr sums to total', Math.abs(m.arr.reduce((a, b) => a + b) - m.tot) < 1e-12);
  const lo = w.macro(1e-20), hi = w.macro(1e3);
  check('clamped outside the data range', lo.tot > 0 && hi.tot > 0);
}

// 3. uncollided transmission = exp(-Sigma T)
for (const [name, M, E, T] of [['water', w, 2, 5], ['water', w, 0.1, 3], ['aluminum', al, 1, 8]]) {
  const sim = MC.makeSim({ material: M, E0: E, T, seed: 11 }), n = 60000;
  for (let i = 0; i < n; i++) MC.runHistory(sim);
  const p = Math.exp(-M.macro(E).tot * T), s = MC.summary(sim);
  check(`${name} E=${E} MeV T=${T} cm: uncollided = exp(-Sigma T)`, Math.abs(s.uncollided.p - p) < 4 * Math.sqrt(p * (1 - p) / n), s.uncollided.p.toFixed(4) + ' vs ' + p.toFixed(4));
  const t = sim.tally;
  check(`${name}: T + R + A (+cutoff) = N`, t.transmitted + t.reflected + t.absorbed === t.N);
  check(`${name}: absorbed = capture + other`, t.absorbed === t.absBy[2] + t.absBy[3]);
}

// 4. generator == fast loop
for (const M of [w, al]) {
  const a = MC.makeSim({ material: M, E0: 3, T: 4, seed: 5 }), b = MC.makeSim({ material: M, E0: 3, T: 4, seed: 5 });
  for (let i = 0; i < 1500; i++) { MC.runHistory(a); for (const ev of MC.historySteps(b)) { } }
  check(`${M.key}: generator == fast loop`, a.tally.transmitted === b.tally.transmitted && a.tally.reflected === b.tally.reflected &&
    a.tally.byReaction.join() === b.tally.byReaction.join() && a.tally.spec.join() === b.tally.spec.join(), a.tally.byReaction.join('/'));
}

// 5. step events for a multi-nuclide material include the nuclide sampling
{
  const sim = MC.makeSim({ material: w, E0: 1, T: 50, seed: 3 }), types = new Set();
  for (const ev of MC.historySteps(sim)) types.add(ev.type);
  check('water history has a nuclide_sample step', types.has('nuclide_sample') && types.has('reaction_sample'), [...types].join(','));
  const s2 = MC.makeSim({ material: al, E0: 1, T: 50, seed: 3 }), t2 = new Set();
  for (const ev of MC.historySteps(s2)) t2.add(ev.type);
  check('single-nuclide material has no nuclide_sample step', !t2.has('nuclide_sample'));
}

// 6. elastic kinematics: E' uniform in [alpha E, E]
for (const A of [1, 15.86, 26.75]) {
  const rng = MC.makeRNG(9), E = 2, n = 200000, al_ = ((A - 1) / (A + 1)) ** 2;
  let s = 0, mn = 1e9, mx = 0, smu = 0;
  for (let i = 0; i < n; i++) { const o = MC.sampleElastic(E, A, rng); s += o.Eout; smu += o.mu; mn = Math.min(mn, o.Eout); mx = Math.max(mx, o.Eout); }
  check(`elastic A=${A}: mean E' = E(1+alpha)/2`, Math.abs(s / n - E * (1 + al_) / 2) < 0.005 * E, (s / n).toFixed(4) + ' vs ' + (E * (1 + al_) / 2).toFixed(4));
  check(`elastic A=${A}: E' in [alpha E, E]`, mn >= E * al_ * (1 - 1e-9) - 1e-12 && mx <= E * (1 + 1e-9), mn.toFixed(4) + '..' + mx.toFixed(4));
  check(`elastic A=${A}: mean lab cosine = 2/(3A)`, Math.abs(smu / n - 2 / (3 * A)) < 0.01, (smu / n).toFixed(4) + ' vs ' + (2 / (3 * A)).toFixed(4));
}

// 7. inelastic sampler: below E - Q, isotropic
{
  const nu = al.nuc[0], rng = MC.makeRNG(4), E = 6, n = 100000;
  let mx = 0, smu = 0;
  for (let i = 0; i < n; i++) { const o = MC.sampleInelastic(E, nu, rng); mx = Math.max(mx, o.Eout); smu += o.mu; }
  check('inelastic E\' < E - Q', mx < E - nu.Qinel + 1e-9, 'Q=' + nu.Qinel.toFixed(3) + ' max=' + mx.toFixed(3));
  check('inelastic emission isotropic', Math.abs(smu / n) < 0.01, (smu / n).toFixed(4));
}

// 8. thermalisation floor and thermal source
{
  const sim = MC.makeSim({ material: w, E0: 1e-6, T: 4, seed: 2 });   // 1 eV source
  let minE = 1e9; sim.record = false;
  for (let i = 0; i < 20000; i++) MC.runHistory(sim);
  const t = sim.tally; let lowest = -1;
  for (let b = 0; b < t.nbins; b++) if (t.spec[b] > 0) { lowest = b; break; }
  const eLow = Math.exp(t.lmin + (t.lmax - t.lmin) * lowest / t.nbins);
  check('no neutron leaves below ~thermal energy', eLow >= 2.0e-8, 'lowest bin edge ' + (eLow * 1e6).toExponential(2) + ' eV');
  check('thermal source: absorbed + reflected + transmitted = N', t.transmitted + t.reflected + t.absorbed === t.N);
}

// 9. thick water: neutrons do not get through
{
  const sim = MC.makeSim({ material: w, E0: 2, T: 80, seed: 8 });
  for (let i = 0; i < 4000; i++) MC.runHistory(sim);
  check('80 cm of water stops 2 MeV neutrons', sim.tally.transmitted === 0 && sim.tally.absorbed + sim.tally.reflected === sim.tally.N, 'absorbed ' + sim.tally.absorbed + ', reflected ' + sim.tally.reflected);
}

// 10. line source and spectrum bins with neutrons
{
  const sim = MC.makeSim({ material: w, E0: 1, T: 5, seed: 7 }); MC.configureSource(sim, 4, 'y', 0);
  for (let i = 0; i < 5000; i++) MC.runHistory(sim);
  check('line source works for neutrons', sim.tally.N === 5000);
  const s = MC.summary(sim);
  check('summary has analytic uncollided', Math.abs(s.analyticUncollided - Math.exp(-s.sigT0 * 5)) < 1e-12);
}

// 11. full library (only when all nuclides are present): every recipe builds, Pb mean free path at 1 MeV
for (const k of Object.keys(MC.N_RECIPES)) {
  const m = have(k) ? mat(k, { water: 1, aluminum: 2.699, iron: 7.874, tungsten: 19.3, lead: 11.35, concrete: 2.3 }[k]) : null;
  if (m) { const t = m.macro(1).tot; check(`${k}: Sigma_t(1 MeV) finite and positive`, t > 0 && isFinite(t), t.toFixed(4) + ' /cm'); }
}
console.log(fails ? fails + ' FAILED' : 'ALL PASSED');
process.exit(fails ? 1 : 0);
