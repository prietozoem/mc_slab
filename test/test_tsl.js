/* Tests of the S(alpha,beta) thermal scattering code. Without an argument a synthetic FREE-GAS S(alpha,beta) of H (analytic) is used;
 * with an argument (path to js/data_tsl.js) the real ENDF/B-VII.1 H(H2O) table is used and the exact free-gas checks are skipped.
 * Run: node test/test_tsl.js [js/data_tsl.js]   (needs the neutron data: test/data_n_mini.js is used for the transport tests) */
'use strict';
global.window = global;
const path = require('path');
require(path.resolve(path.join(__dirname, 'data_n_mini.js')));
const MC = require('../js/engine.js');
let fails = 0;
function check(name, ok, info) { console.log((ok ? 'PASS ' : 'FAIL ') + name + (info !== undefined ? '  [' + info + ']' : '')); if (!ok) fails++; }

function freeGasData() {
  const beta = [0], alpha = [];
  for (let b = 0.01; b < 158.2; b *= 1.03) beta.push(+b.toPrecision(5));
  for (let a = 0.001; a < 633; a *= 1.06) alpha.push(+a.toPrecision(5));
  const A0 = 0.99917, nb = beta.length, na = alpha.length, q = new Uint16Array(nb * na);
  for (let k = 0; k < nb; k++) for (let j = 0; j < na; j++) {
    const a = alpha[j], b = beta[k], S = Math.exp(-(a * a + b * b) / (4 * a)) / Math.sqrt(4 * Math.PI * a);
    q[k * na + j] = S > 0 ? Math.max(1, Math.min(65535, Math.round((Math.log(S) + 235) / 0.004))) : 0;
  }
  return { T: 293.6, kT: 0.0253, A0, sigFree: 20.4363, Emax: 4.0, lnOff: -235, lnStep: 0.004, nb, na, beta, alpha, q: Buffer.from(q.buffer).toString('base64') };
}
let real = false, data;
if (process.argv[2]) { require(path.resolve(process.argv[2])); data = window.TSL_H2O; real = true; } else data = freeGasData();
const T = MC.makeThermal(data);
const erf = (x) => { const t = 1 / (1 + 0.3275911 * x); return 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); };
const fg = (E) => { const x = Math.sqrt(data.A0 * E / data.kT); return data.sigFree * ((1 + 1 / (2 * x * x)) * erf(x) + Math.exp(-x * x) / (x * Math.sqrt(Math.PI))); };

// 1. cross section
if (!real) for (const E of [0.001, 0.0253, 0.1, 1]) {
  const s = T.sigma(E * 1e-6), ex = fg(E);
  check(`free-gas sigma(${E} eV)`, Math.abs(s / ex - 1) < 0.03, s.toFixed(2) + ' vs ' + ex.toFixed(2));
} else {
  const s = T.sigma(0.0253e-6);
  check('H in H2O sigma(0.0253 eV) in the expected range (45-58 b)', s > 45 && s < 58, s.toFixed(2) + ' b');
  check('H in H2O sigma(4 eV) close to the free-atom value 20.4 b', Math.abs(T.sigma(3.99e-6) / 20.44 - 1) < 0.05, T.sigma(3.99e-6).toFixed(2));
}
check('sigma decreases with energy (1/v then flat)', T.sigma(1e-8 * 1e-0) > T.sigma(1e-7) && T.sigma(1e-7) > T.sigma(1e-6) * 0.99);

// 2. kinematics: E' and mu consistent with alpha, beta; E' > 0
{
  const rng = MC.makeRNG(3); let ok = true, sm = 0, n = 20000, up = 0;
  for (let i = 0; i < n; i++) {
    const E = 0.01e-6 * Math.pow(300, rng()) * 1;      // 0.01 .. 3 eV
    const o = T.sample(E, rng);
    ok = ok && o.Eout > 0 && o.mu >= -1 && o.mu <= 1 && isFinite(o.mu);
    sm += o.mu; if (o.Eout > E) up++;
  }
  check('samples valid (E\'>0, |mu|<=1)', ok);
  check('both up- and down-scattering occur', up > 0.05 * n && up < 0.95 * n, 'up fraction ' + (up / n).toFixed(2));
}

// 3. detailed balance: the stationary distribution of the scattering chain is  pi(E) ~ sigma(E) E exp(-E/kT)
{
  const rng = MC.makeRNG(11), kT = data.kT * 1e-6; let E = 0.5e-6, sInv = 0, sE = 0, n = 300000, burn = 2000, cnt = 0;
  for (let i = 0; i < n; i++) {
    const o = T.sample(E, rng); E = o.Eout;
    if (E >= T.EmaxMeV) E = T.EmaxMeV * 0.999;        // above Emax the table does not apply: keep the chain inside
    if (i > burn) { sInv += 1 / T.sigma(E); cnt++; }
  }
  // expectation of 1/sigma under pi:  int E e^{-E/kT} dE / int sigma E e^{-E/kT} dE  (numerical, log grid)
  let num = 0, den = 0; const m = 4000, lo = Math.log(1e-11), hi = Math.log(T.EmaxMeV);
  for (let i = 0; i < m; i++) { const e = Math.exp(lo + (hi - lo) * (i + 0.5) / m), w = e * Math.exp(-e / kT) * e; num += w; den += w * T.sigma(e); }
  const exact = num / den, mc = sInv / cnt;
  check('detailed balance (equilibrium <1/sigma>)', Math.abs(mc / exact - 1) < 0.08, mc.toFixed(5) + ' vs ' + exact.toFixed(5));
}

// 4. water slab with S(alpha,beta): transport conservation, uncollided = exp(-Sigma T), generator == fast loop, no 0.0253 eV floor
{
  const rho = 1, mat = MC.makeNeutronMaterial('water', rho, {}, window.NXS, { tsl: data }), matNo = MC.makeNeutronMaterial('water', rho, {}, window.NXS);
  const sH = mat.nuc[0];
  check('H-1 nuclide carries the thermal law', !!sH.tsl && !matNo.nuc[0].tsl);
  const m0 = mat.macro(1e-8), m1 = matNo.macro(1e-8);
  check('thermal law changes the macroscopic elastic cross section at 0.01 eV', Math.abs(m0.el / m1.el - 1) > 0.02, m0.el.toFixed(3) + ' vs ' + m1.el.toFixed(3));
  const m4 = mat.macro(5e-6), m4b = matNo.macro(5e-6);
  check('above Emax both materials agree', m4.tot === m4b.tot);
  const Tc = 3, E0 = 0.05e-6, n = 40000, sim = MC.makeSim({ material: mat, E0, T: Tc, seed: 21 });
  for (let i = 0; i < n; i++) MC.runHistory(sim);
  const t = sim.tally, p = Math.exp(-mat.macro(E0).tot * Tc), s = MC.summary(sim);
  check('uncollided = exp(-Sigma T)', Math.abs(s.uncollided.p - p) < 4 * Math.sqrt(p * (1 - p) / n), s.uncollided.p.toFixed(4) + ' vs ' + p.toFixed(4));
  check('conservation', t.transmitted + t.reflected + t.absorbed === t.N);
  const a = MC.makeSim({ material: mat, E0: 0.02e-6, T: 2, seed: 5 }), b = MC.makeSim({ material: mat, E0: 0.02e-6, T: 2, seed: 5 });
  for (let i = 0; i < 800; i++) { MC.runHistory(a); for (const ev of MC.historySteps(b)) { } }
  check('generator == fast loop (S(alpha,beta))', a.tally.transmitted === b.tally.transmitted && a.tally.spec.join() === b.tally.spec.join() && a.tally.byReaction.join() === b.tally.byReaction.join());
  // thermalised neutrons: reflected spectrum peaks near kT, and some are below 0.0253 eV
  const sim2 = MC.makeSim({ material: mat, E0: 1e-6, T: 30, seed: 8 });
  for (let i = 0; i < 6000; i++) MC.runHistory(sim2);
  const ta = sim2.tally; let lowBin = ta.nbins;
  for (let k = 0; k < ta.nbins; k++) if (ta.spec[k] > 0) { lowBin = k; break; }
  check('thick water (30 cm): almost nothing is transmitted; T + R + A = N', ta.transmitted < 0.005 * ta.N && ta.transmitted + ta.reflected + ta.absorbed === ta.N, 'T ' + ta.transmitted + ' R ' + ta.reflected + ' A ' + ta.absorbed);
}
console.log(fails ? fails + ' FAILED' : 'ALL PASSED');
process.exit(fails ? 1 : 0);
