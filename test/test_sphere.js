/* Tests of the sphere example (isotropic neutron point source at the centre of a sphere).
 * Uses test/data_n_mini.js (tiny subset). Run: node test/test_sphere.js [path/to/data_n.js] */
'use strict';
global.window = global;
const path = require('path');
require(path.resolve(process.argv[2] || path.join(__dirname, 'data_n_mini.js')));
const MC = require('../js/engine.js');
let fails = 0;
function check(name, ok, info) { console.log((ok ? 'PASS ' : 'FAIL ') + name + (info !== undefined ? '  [' + info + ']' : '')); if (!ok) fails++; }
const water = (rho, tsl) => MC.makeNeutronMaterial('water', rho === undefined ? 1 : rho, {}, window.NXS, tsl ? { tsl } : undefined);
const sphere = (mat, E0, R, seed, n) => { const sim = MC.makeSim({ material: mat, E0, R, geom: 'sphere', seed, nbins: 200 }); for (let i = 0; i < n; i++) MC.runHistory(sim); return sim; };

// 1. geometry: distance to the surface from a point inside
{
  const d = MC.distSphere({ x: 0, y: 0, z: 0, u: 0, v: 1, w: 0 }, 5), d2 = MC.distSphere({ x: 3, y: 0, z: 0, u: 1, v: 0, w: 0 }, 5), d3 = MC.distSphere({ x: 3, y: 0, z: 0, u: -1, v: 0, w: 0 }, 5);
  check('distance to the surface (centre / towards / away)', Math.abs(d - 5) < 1e-12 && Math.abs(d2 - 2) < 1e-12 && Math.abs(d3 - 8) < 1e-12, [d, d2, d3].join(','));
}

// 2. vacuum-like sphere: every neutron flies a straight line of length R -> volume-averaged flux = R/V = 3/(4 pi R^2), all in the E0 bin
{
  const R = 10, E0 = 1.5, n = 4000, sim = sphere(water(1e-7), E0, R, 3, n), s = MC.summary(sim), t = sim.tally;
  const ex = 3 / (4 * Math.PI * R * R);
  check('vacuum: mean track length = R', Math.abs(t.trk / n - R) < 1e-3 * R, (t.trk / n).toFixed(4));
  check('vacuum: volume-averaged flux = 3/(4 pi R^2)', Math.abs(s.flux / ex - 1) < 1e-3, s.flux.toExponential(4) + ' vs ' + ex.toExponential(4));
  check('vacuum: everything leaks uncollided', t.leaked === n && t.uncollided === n && t.absorbed === 0);
  const sp = MC.fluxSpectrum(sim), tot = sp.reduce((a, b) => a + b.y * Math.log(b.e1 / b.e0), 0), hit = sp.filter(b => b.y > 0);
  check('vacuum: flux spectrum sits in one bin containing E0 and integrates to the total', hit.length === 1 && hit[0].e0 <= E0 && E0 < hit[0].e1 && Math.abs(tot / s.flux - 1) < 1e-9, hit.length + ' bin(s)');
}

// 3. water sphere: conservation, uncollided leakage = exp(-Sigma R), absorbed = capture + other
for (const [E, R] of [[2, 5], [0.05, 10], [1e-7, 3]]) {
  const M = water(1), n = 40000, sim = sphere(M, E, R, 11, n), s = MC.summary(sim), t = sim.tally, p = Math.exp(-M.macro(E).tot * R);
  check(`E=${E} MeV, R=${R} cm: leaked + absorbed = N`, t.leaked + t.absorbed === t.N);
  check(`E=${E} MeV, R=${R} cm: uncollided leakage = exp(-Sigma R)`, Math.abs(s.uncollided.p - p) < 4 * Math.sqrt(p * (1 - p) / n) + 1e-9, s.uncollided.p.toFixed(4) + ' vs ' + p.toFixed(4));
  check(`E=${E} MeV, R=${R} cm: absorbed = capture + other`, t.absorbed === t.absBy[2] + t.absBy[3]);
  const spec = Array.from(t.spec).reduce((a, b) => a + b, 0);
  check(`E=${E} MeV, R=${R} cm: leaked spectrum counts collided leaks`, spec === t.leaked - t.uncollided);
}

// 4. track-length estimator = collision estimator:  mean number of collisions = Sum_i (track in bin i) * Sigma_t(E_i)  (bin-resolved check with a fine, nearly flat Sigma_t)
{
  const M = water(1), E = 14, R = 8, n = 30000, sim = sphere(M, E, R, 5, n), t = sim.tally;
  let est = 0; const dl = (t.lmax - t.lmin) / t.nbins;
  for (let i = 0; i < t.nbins; i++) est += (t.fl[i] / n) * M.macro(Math.exp(t.lmin + (i + 0.5) * dl)).tot;
  const coll = t.coll / n;
  check('track-length x Sigma_t (binned) agrees with the collision count', Math.abs(est / coll - 1) < 0.06, est.toFixed(3) + ' vs ' + coll.toFixed(3));
}

// 5. generator and fast loop agree for the same seed (sphere)
{
  const M = water(1), a = MC.makeSim({ material: M, E0: 2, R: 6, geom: 'sphere', seed: 9, nbins: 200 }), b = MC.makeSim({ material: M, E0: 2, R: 6, geom: 'sphere', seed: 9, nbins: 200 });
  let events = 0, birthIso = true;
  for (let i = 0; i < 600; i++) { MC.runHistory(a); for (const ev of MC.historySteps(b)) { events++; if (ev.type === 'birth' && !(ev.iso && ev.mu >= -1 && ev.mu <= 1)) birthIso = false; } }
  const ta = a.tally, tb = b.tally;
  check('generator == fast loop (sphere)', ta.leaked === tb.leaked && ta.absorbed === tb.absorbed && ta.byReaction.join() === tb.byReaction.join() && ta.spec.join() === tb.spec.join() && ta.fl.join() === tb.fl.join() && ta.trk === tb.trk, events + ' events');
  check('birth events carry the isotropic direction', birthIso);
}

// 6. isotropy of the source: <cos theta> = 0, <cos^2 theta> = 1/3 for the first flight directions
{
  const sim = MC.makeSim({ material: water(1e-7), E0: 1, R: 5, geom: 'sphere', seed: 4, nbins: 50 }); let s1 = 0, s2 = 0, n = 60000;
  sim.record = true;
  for (let i = 0; i < n; i++) { const r = MC.runHistory(sim); s1 += r.ph.mu0; s2 += r.ph.mu0 * r.ph.mu0; }
  check('isotropic source: <mu> = 0 and <mu^2> = 1/3', Math.abs(s1 / n) < 0.012 && Math.abs(s2 / n - 1 / 3) < 0.01, (s1 / n).toFixed(4) + ', ' + (s2 / n).toFixed(4));
}

// 7. slab example unaffected by the new fields
{
  const sim = MC.makeSim({ material: water(1), E0: 2, T: 4, seed: 2 }); for (let i = 0; i < 2000; i++) MC.runHistory(sim);
  const t = sim.tally; check('slab: T + R + A = N and no sphere tallies', !sim.sphere && t.transmitted + t.reflected + t.absorbed === t.N && t.leaked === 0 && t.trk === 0);
}
console.log(fails ? fails + ' FAILED' : 'ALL PASSED');
process.exit(fails ? 1 : 0);
