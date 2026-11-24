/* Monte Carlo photon transport in an infinite slab -- teaching engine.
 *
 * Geometry : slab 0 <= x <= T (infinite in y, z). Parallel monoenergetic beam
 *            at normal incidence on the face x = 0, travelling along +x.
 * Physics  : coherent (Rayleigh), incoherent (Compton), photoelectric, pair production.
 *            Cross sections from NIST XCOM, log-log interpolation.
 *            - Photoelectric: photon absorbed (no fluorescence, no electron transport).
 *            - Pair production: photon absorbed (annihilation photons NOT tracked).
 *            - Incoherent: Klein-Nishina on a free electron at rest (Kahn's rejection).
 *            - Coherent: Thomson x (simple Gaussian form factor) -> forward peaked.
 *
 * Example 2 (neutrons only): isotropic point source at the centre of a homogeneous sphere of radius R in
 * vacuum (makeSim with geom:'sphere'). Tallies: leakage, absorption and the volume-averaged flux spectrum in
 * the sphere (track-length estimator, with per-bin statistical errors from the history-to-history scatter).
 *
 * The same physics functions are used by the fast loop (runHistory) and by the
 * step-by-step generator (historySteps). With the same seed both produce identical results.
 */
(function (root) {
  'use strict';

  var ME = 0.51099895;                 // electron rest energy [MeV]
  var REACTIONS = ['coherent', 'incoherent', 'photoelectric', 'pair'];
  var ZEFF = { water: 7.4, aluminum: 13, iron: 26, tungsten: 74, lead: 82, concrete: 11 };
  var E_CUTOFF = 0.001;                // MeV: lowest energy of the XCOM tables used

  /* ---------- random numbers (mulberry32) ---------- */
  function makeRNG(seed) {
    var a = (seed >>> 0) || 1;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      var t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;   // in [0,1)
    };
  }

  /* ---------- material with XCOM tables ---------- */
  function Material(key, data) {
    this.key = key;
    this.rho = data.rho;
    this.rows = data.d;                // [E, coh, inc, pe, pairN, pairE] in MeV, cm2/g
    this.zeff = ZEFF[key] || 10;
    this.name = data.name;
  }

  function interp(a, b, E, k) {
    var y1 = a[k], y2 = b[k];
    if (y1 > 0 && y2 > 0) {
      var l1 = Math.log(a[0]), l2 = Math.log(b[0]);
      return Math.exp(Math.log(y1) + (Math.log(y2) - Math.log(y1)) * (Math.log(E) - l1) / (l2 - l1));
    }
    return y1 + (y2 - y1) * (E - a[0]) / (b[0] - a[0]);   // linear when a value is zero (pair threshold)
  }

  /* mass cross sections [coh, inc, pe, pair] in cm2/g */
  Material.prototype.massXS = function (E) {
    var rows = this.rows, n = rows.length, out = [0, 0, 0, 0];
    if (E <= rows[0][0]) E = rows[0][0];
    if (E >= rows[n - 1][0]) {
      var r = rows[n - 1];
      return [r[1], r[2], r[3], r[4] + r[5]];
    }
    var lo = 0, hi = n - 1;
    while (hi - lo > 1) {                       // rows[lo][0] <= E < rows[hi][0]; at an edge, lo = post-edge row
      var mid = (lo + hi) >> 1;
      if (rows[mid][0] <= E) lo = mid; else hi = mid;
    }
    var a = rows[lo], b = rows[hi];
    out[0] = interp(a, b, E, 1);
    out[1] = interp(a, b, E, 2);
    out[2] = interp(a, b, E, 3);
    out[3] = interp(a, b, E, 4) + interp(a, b, E, 5);
    return out;
  };

  /* macroscopic cross sections [1/cm] */
  Material.prototype.macro = function (E) {
    var m = this.massXS(E), rho = this.rho;
    var c = m[0] * rho, i = m[1] * rho, p = m[2] * rho, q = m[3] * rho;
    return { coh: c, inc: i, pe: p, pair: q, tot: c + i + p + q, arr: [c, i, p, q] };
  };

  /* curves of the table for plotting (macroscopic, 1/cm) */
  Material.prototype.curves = function () {
    var rho = this.rho, E = [], c = [], i = [], p = [], q = [], t = [];
    for (var k = 0; k < this.rows.length; k++) {
      var r = this.rows[k];
      E.push(r[0]); c.push(r[1] * rho); i.push(r[2] * rho); p.push(r[3] * rho);
      q.push((r[4] + r[5]) * rho);
      t.push((r[1] + r[2] + r[3] + r[4] + r[5]) * rho);
    }
    return { E: E, coh: c, inc: i, pe: p, pair: q, tot: t };
  };

  /* ---------- sampling routines ---------- */
  function sampleFlight(xi, sigT) { return -Math.log(1 - xi) / sigT; }   // inverse of F(s) = 1 - exp(-Sigma s)

  function sampleReaction(xi, sig) {
    var x = xi * sig.tot;
    if (x < sig.coh) return 0;
    if (x < sig.coh + sig.inc) return 1;
    if (x < sig.coh + sig.inc + sig.pe) return 2;
    return 3;
  }

  /* Klein-Nishina, free electron at rest, Kahn's method. Returns outgoing energy and cos(theta). */
  function sampleCompton(E, rng) {
    var a = E / ME, tries = 0, p1 = (1 + 2 * a) / (9 + 2 * a);
    for (;;) {
      tries++;
      var r1 = rng(), r2 = rng(), r3 = rng(), x, ok;
      if (r1 <= p1) {
        x = 1 + 2 * a * r2;                                   // x = E / E'
        ok = r3 <= 4 * (1 / x - 1 / (x * x));
      } else {
        x = (1 + 2 * a) / (1 + 2 * a * r2);
        var mu = 1 - (x - 1) / a;
        ok = r3 <= 0.5 * (mu * mu + 1 / x);
      }
      if (ok) return { Eout: E / x, mu: 1 - (x - 1) / a, tries: tries };
    }
  }

  /* Coherent: Thomson (1+mu^2)/2 times F^2 = exp(-(q/q0)^2), q = (E/hc) sin(theta/2) [1/Angstrom],
     q0 = 0.63 Z^(1/3). Crude, but forward peaked and energy/Z dependent. */
  function sampleRayleigh(E, zeff, rng) {
    var q0 = 0.63 * Math.cbrt(zeff);
    var k = E * 1000 / 12.3984;
    var b = k * k / (2 * q0 * q0);                 // F^2 = exp(-b (1 - mu))
    var tries = 0;
    for (;;) {
      tries++;
      var xi = rng(), y;
      y = b < 1e-9 ? 2 * xi : -Math.log(1 - xi * (1 - Math.exp(-2 * b))) / b;
      if (y > 2) y = 2;
      var mu = 1 - y;
      if (rng() <= 0.5 * (1 + mu * mu)) return { mu: mu, tries: tries };
    }
  }

  /* unnormalised angular distributions (used for the teaching plots) */
  function kleinNishinaPDF(E, mu) {
    var a = E / ME, r = 1 / (1 + a * (1 - mu));
    return r * r * (r + 1 / r - (1 - mu * mu));
  }
  function rayleighPDF(E, zeff, mu) {
    var q0 = 0.63 * Math.cbrt(zeff), k = E * 1000 / 12.3984, b = k * k / (2 * q0 * q0);
    return 0.5 * (1 + mu * mu) * Math.exp(-b * (1 - mu));
  }

  function rotate(u, v, w, mu, phi) {
    var s = Math.sqrt(Math.max(0, 1 - mu * mu)), cp = Math.cos(phi), sp = Math.sin(phi), un, vn, wn;
    if (Math.abs(w) > 0.99999) {
      var sg = w > 0 ? 1 : -1;
      un = s * cp; vn = sg * s * sp; wn = sg * mu;
    } else {
      var d = Math.sqrt(1 - w * w);
      un = mu * u + s * (u * w * cp - v * sp) / d;
      vn = mu * v + s * (v * w * cp + u * sp) / d;
      wn = mu * w - s * d * cp;
    }
    var n = Math.sqrt(un * un + vn * vn + wn * wn);
    return [un / n, vn / n, wn / n];
  }

  /* ---------- tallies ---------- */
  var NE_MIN = 1e-11, NE_MAX = 20;      // neutron energy range of the ENDF data [MeV] (1e-5 eV .. 20 MeV)
  function makeTally(E0, nbins, logbins) {
    nbins = nbins || 100;
    return {
      log: !!logbins, lmin: Math.log(NE_MIN), lmax: Math.log(NE_MAX),
      E0: E0, N: 0,
      transmitted: 0, uncollided: 0, reflected: 0, absorbed: 0,
      absBy: [0, 0, 0, 0],                // absorbed by reaction (index 2 = photoelectric, 3 = pair)
      byReaction: [0, 0, 0, 0],           // interactions of each type
      cutoff: 0,
      flux: 0, flux2: 0,                  // sum of 1/u at the exit face (scalar flux estimator)
      leaked: 0, coll: 0,                 // sphere example: neutrons that left the sphere, total number of collisions
      fl: new Float64Array(nbins), fl2: new Float64Array(nbins), hist: new Float64Array(nbins),   // sphere: track length per bin (sum, sum of squares over histories, current history)
      trk: 0, trk2: 0, ht: 0,             // sphere: total track length (sum, sum of squares, current history)
      nbins: nbins, spec: new Float64Array(nbins)   // collided transmitted particles: linear bins in [0, E0] (photons) or log bins in [1e-11, 20] MeV (neutrons)
    };
  }

  /* ---------- simulation object ---------- */
  function makeSim(opts) {
    var mat = opts.material, sph = opts.geom === 'sphere';
    return {
      sphere: sph, R: sph ? opts.R : 0,
      mat: mat, E0: opts.E0, T: sph ? opts.R : opts.T, rng: makeRNG(opts.seed || 1),
      srcLen: 0, srcAxis: 'z', srcX0: 0,   // line source: length [cm] (0 = pencil beam), axis, position x0 (see configureSource)
      tally: makeTally(opts.E0, opts.nbins, opts.material.particle === 'neutron'),
      neutron: opts.material.particle === 'neutron', record: false, history: 0
    };
  }
  function configureSource(sim, len, axis, x0) {
    var T = sim.T;
    sim.srcAxis = axis === 'x' || axis === 'y' ? axis : 'z';
    sim.srcX0 = Math.min(Math.max(x0 || 0, 0), T);
    len = len > 0 ? len : 0;
    if (sim.srcAxis === 'x') len = Math.min(len, T - sim.srcX0);   // the segment must lie inside the slab
    sim.srcLen = len;
    return sim;
  }

  /* Source, monodirectional along +x. srcLen = 0: point at the entrance face (pencil beam).
     srcLen > 0: emission point uniform on a line segment (one random number):
       axis 'x': x0 <= x <= x0 + L  (y = z = 0)   -> distributed in depth (dX)
       axis 'y': x = x0, -L/2 <= y <= L/2, z = 0
       axis 'z': x = x0, y = 0, -L/2 <= z <= L/2 */
  function newPhoton(E0, sim) {
    var ph = { x: 0, y: 0, z: 0, u: 1, v: 0, w: 0, E: E0, nc: 0, path: null };
    if (sim && sim.sphere) {                          // isotropic point source at the centre: cos(theta) = 2 xi1 - 1, phi = 2 pi xi2 (u is the displayed z axis)
      var m0 = 2 * sim.rng() - 1, f0 = 2 * Math.PI * sim.rng(), s0 = Math.sqrt(Math.max(0, 1 - m0 * m0));
      ph.u = m0; ph.v = s0 * Math.cos(f0); ph.w = s0 * Math.sin(f0); ph.mu0 = m0; ph.phi0 = f0;
      return ph;
    }
    if (sim && sim.srcLen > 0) {
      var r = sim.rng(), L = sim.srcLen;
      ph.x = sim.srcX0;
      if (sim.srcAxis === 'x') ph.x = sim.srcX0 + r * L;
      else if (sim.srcAxis === 'y') ph.y = (r - 0.5) * L;
      else ph.z = (r - 0.5) * L;
    }
    return ph;
  }

  /* move the photon a distance s; returns 'T', 'R' or null (stays inside) */
  function transport(ph, s, T) {
    var xn = ph.x + ph.u * s, f;
    if (xn >= T) {
      f = (T - ph.x) / ph.u;
      ph.y += ph.v * f; ph.z += ph.w * f; ph.x = T;
      if (ph.path) ph.path.push([ph.x, ph.y, 'T', ph.z]);
      return 'T';
    }
    if (xn < 0) {
      f = (0 - ph.x) / ph.u;
      ph.y += ph.v * f; ph.z += ph.w * f; ph.x = 0;
      if (ph.path) ph.path.push([ph.x, ph.y, 'R', ph.z]);
      return 'R';
    }
    ph.x = xn; ph.y += ph.v * s; ph.z += ph.w * s;
    return null;
  }

  /* ---- sphere geometry: distance to the surface along the flight direction, flight with track-length scoring ---- */
  function distSphere(ph, R) {
    var b = ph.x * ph.u + ph.y * ph.v + ph.z * ph.w, c = ph.x * ph.x + ph.y * ph.y + ph.z * ph.z - R * R;
    return -b + Math.sqrt(Math.max(0, b * b - c));
  }
  function ebin(t, E) {
    var b = Math.floor((Math.log(Math.max(E, NE_MIN)) - t.lmin) / (t.lmax - t.lmin) * t.nbins);
    return b < 0 ? 0 : b >= t.nbins ? t.nbins - 1 : b;
  }
  function transportS(sim, ph, s) {                   // returns 'L' (leaked) or null; scores the track length at the flight energy
    var db = distSphere(ph, sim.R), out = null, len = s, t = sim.tally;
    if (s >= db) { len = db; out = 'L'; }
    t.hist[ebin(t, ph.E)] += len; t.ht += len;
    ph.x += ph.u * len; ph.y += ph.v * len; ph.z += ph.w * len;
    if (out && ph.path) ph.path.push([ph.x, ph.y, 'L', ph.z]);
    return out;
  }

  /* perform interaction r on photon; returns info {absorbed, ...} */
  function interact(sim, ph, r) {
    var rng = sim.rng, info = { r: r, absorbed: false, Ein: ph.E };
    ph.nc++;
    sim.tally.byReaction[r]++;
    if (r === 2 || r === 3) {                         // photoelectric / pair: absorbed
      info.absorbed = true;
      if (ph.path) ph.path.push([ph.x, ph.y, r, ph.z]);
      return info;
    }
    var mu, tries;
    if (r === 1) {
      var c = sampleCompton(ph.E, rng);
      mu = c.mu; tries = c.tries; info.Eout = c.Eout;
    } else {
      var q = sampleRayleigh(ph.E, sim.mat.zeff, rng);
      mu = q.mu; tries = q.tries; info.Eout = ph.E;
    }
    var phi = 2 * Math.PI * rng();
    var d = rotate(ph.u, ph.v, ph.w, mu, phi);
    ph.u = d[0]; ph.v = d[1]; ph.w = d[2]; ph.E = info.Eout;
    info.mu = mu; info.phi = phi; info.tries = tries; info.dir = d;
    if (ph.path) ph.path.push([ph.x, ph.y, r, ph.z]);
    return info;
  }

  function finish(sim, ph, out) {
    var t = sim.tally, i;
    if (sim.sphere) {
      t.N++; sim.history++; t.coll += ph.nc;
      for (i = 0; i < t.nbins; i++) { var h = t.hist[i]; if (h !== 0) { t.fl[i] += h; t.fl2[i] += h * h; t.hist[i] = 0; } }
      t.trk += t.ht; t.trk2 += t.ht * t.ht; t.ht = 0;
      if (out === 'L') {
        t.leaked++;
        if (ph.nc === 0) t.uncollided++; else t.spec[ebin(t, ph.E)]++;
      } else { t.absorbed++; if (out === 'C') t.cutoff++; }
      return;
    }
    t.N++; sim.history++;
    if (out === 'T') {
      t.transmitted++;
      var inv = 1 / ph.u;
      t.flux += inv; t.flux2 += inv * inv;
      if (ph.nc === 0) t.uncollided++;
      else {
        var b = t.log ? Math.floor((Math.log(Math.max(ph.E, NE_MIN)) - t.lmin) / (t.lmax - t.lmin) * t.nbins)
                      : Math.floor(ph.E / t.E0 * t.nbins);
        t.spec[Math.max(0, Math.min(t.nbins - 1, b))]++;
      }
    } else if (out === 'R') t.reflected++;
    else { t.absorbed++; if (out === 'C') t.cutoff++; }
  }

  /* ---------- fast history ---------- */
  function runHistory(sim) {
    if (sim.neutron) return runHistoryN(sim);
    var ph = newPhoton(sim.E0, sim), mat = sim.mat, rng = sim.rng, T = sim.T;
    if (sim.record) ph.path = [[ph.x, ph.y, 'S', ph.z]];
    for (;;) {
      if (ph.E < E_CUTOFF) { finish(sim, ph, 'C'); return { out: 'C', ph: ph }; }
      var sig = mat.macro(ph.E);
      var s = sampleFlight(rng(), sig.tot);
      var out = transport(ph, s, T);
      if (out) { finish(sim, ph, out); return { out: out, ph: ph }; }
      var r = sampleReaction(rng(), sig);
      var info = interact(sim, ph, r);
      if (info.absorbed) { sim.tally.absBy[r]++; finish(sim, ph, 'A'); return { out: 'A', ph: ph }; }
    }
  }

  /* ---------- step-by-step history (generator) ---------- */
  function* historySteps(sim) {
    if (sim.neutron) { yield* historyStepsN(sim); return; }
    var ph = newPhoton(sim.E0, sim), mat = sim.mat, rng = sim.rng, T = sim.T;
    ph.path = [[ph.x, ph.y, 'S', ph.z]];
    yield { type: 'birth', ph: ph, n: sim.history + 1, srcLen: sim.srcLen, srcAxis: sim.srcAxis };
    for (;;) {
      if (ph.E < E_CUTOFF) {
        finish(sim, ph, 'C');
        yield { type: 'end', out: 'C', ph: ph };
        return;
      }
      var sig = mat.macro(ph.E);
      yield { type: 'xs_total', ph: ph, sig: sig };
      var xi = rng(), s = sampleFlight(xi, sig.tot);
      var db = ph.u > 0 ? (T - ph.x) / ph.u : (0 - ph.x) / ph.u;   // distance to the slab boundary along the flight direction
      yield { type: 'flight', ph: ph, sig: sig, xi: xi, s: s, db: db };
      var from = [ph.x, ph.y];
      var out = transport(ph, s, T);
      yield { type: 'move', ph: ph, from: from, s: s, db: db, out: out };
      if (out) {
        finish(sim, ph, out);
        yield { type: 'end', out: out, ph: ph };
        return;
      }
      yield { type: 'xs_partial', ph: ph, sig: sig };
      var xr = rng(), r = sampleReaction(xr, sig);
      yield { type: 'reaction_sample', ph: ph, sig: sig, xi: xr, r: r };
      var info = interact(sim, ph, r);
      yield { type: 'interaction', ph: ph, info: info };
      if (info.absorbed) {
        sim.tally.absBy[r]++;
        finish(sim, ph, 'A');
        yield { type: 'end', out: 'A', r: r, ph: ph };
        return;
      }
    }
  }


  /* =====================================================================================
   *  NEUTRONS  (ENDF/B-VII.1 pointwise cross sections, see js/data_n.js)
   *  Per nuclide: elastic (MT2), inelastic (MT4), radiative capture (MT102) and "other"
   *  absorption (MT1 - MT2 - MT4 - MT102: (n,2n), (n,p), (n,alpha)... treated as absorbed).
   *  Simple physics: target at rest, elastic isotropic in the centre of mass, inelastic
   *  isotropic with an evaporation spectrum, pseudo-thermalisation at 0.0253 eV.
   * ===================================================================================== */
  var N_REACTIONS = ['elastic', 'inelastic', 'capture', 'other'];
  var E_THERMAL = 2.53e-8;             // MeV (0.0253 eV)
  var M_NEUTRON = 1.00866491588;       // amu
  var NA_B = 0.602214076;              // Avogadro x 1e-24: atoms per (mol b cm)^-1 ... N [atoms/(b cm)] = rho * NA_B * n / sum(n M)

  /* material recipes: relative atom counts n (after lumping minor isotopes) and mean atomic mass M [amu] of the element */
  function concretePart(id, w, M, f) { return { id: id, n: w / M * (f || 1), M: M }; }
  var N_RECIPES = {
    water: [{ id: 'h1', n: 2, M: 1.00794 }, { id: 'o16', n: 1, M: 15.9994 }],
    aluminum: [{ id: 'al27', n: 1, M: 26.9815 }],
    iron: [{ id: 'fe54', n: 0.05845, M: 55.845 }, { id: 'fe56', n: 0.91754, M: 55.845 },
           { id: 'fe57', n: 0.02119, M: 55.845 }, { id: 'fe58', n: 0.00282, M: 55.845 }],
    tungsten: [{ id: 'w182', n: 0.2662, M: 183.84 }, { id: 'w183', n: 0.1431, M: 183.84 },
               { id: 'w184', n: 0.3064, M: 183.84 }, { id: 'w186', n: 0.2843, M: 183.84 }],          // W-180 (0.12 %) lumped into W-182
    lead: [{ id: 'pb204', n: 0.014, M: 207.2 }, { id: 'pb206', n: 0.241, M: 207.2 },
           { id: 'pb207', n: 0.221, M: 207.2 }, { id: 'pb208', n: 0.524, M: 207.2 }],
    concrete: [                                                   // NBS ordinary concrete, mass fractions; minor isotopes lumped
      concretePart('h1', 0.0221, 1.00794), concretePart('c0', 0.002484, 12.011), concretePart('o16', 0.57493, 15.9994),
      concretePart('na23', 0.015208, 22.98977), concretePart('mg24', 0.001266, 24.305), concretePart('al27', 0.019953, 26.9815),
      concretePart('si28', 0.304627, 28.0855, 0.92223), concretePart('si29', 0.304627, 28.0855, 0.04685), concretePart('si30', 0.304627, 28.0855, 0.03092),
      concretePart('k39', 0.010045, 39.0983, 0.93258), concretePart('k41', 0.010045, 39.0983, 0.06730),
      concretePart('ca40', 0.042951, 40.078), concretePart('fe56', 0.006435, 55.845)
    ]
  };
  var N_LABEL = { h1: 'H-1', c0: 'C', o16: 'O-16', al27: 'Al-27', fe54: 'Fe-54', fe56: 'Fe-56', fe57: 'Fe-57', fe58: 'Fe-58',
    w182: 'W-182', w183: 'W-183', w184: 'W-184', w186: 'W-186', pb204: 'Pb-204', pb206: 'Pb-206', pb207: 'Pb-207', pb208: 'Pb-208',
    na23: 'Na-23', mg24: 'Mg-24', si28: 'Si-28', si29: 'Si-29', si30: 'Si-30', k39: 'K-39', k41: 'K-41', ca40: 'Ca-40' };

  function prepCurve(cv) {                                        // cv = {e: [MeV], s: [b]}  ->  typed arrays + logs
    if (!cv) return null;
    var n = cv.e.length, le = new Float64Array(n), ls = new Float64Array(n), e = new Float64Array(n), s = new Float64Array(n);
    for (var i = 0; i < n; i++) { e[i] = cv.e[i]; s[i] = cv.s[i]; le[i] = Math.log(e[i]); ls[i] = s[i] > 0 ? Math.log(s[i]) : NaN; }
    return { e: e, s: s, le: le, ls: ls, n: n };
  }
  /* cross section (b) of a prepared curve at energy E [MeV]: log-log between points, linear if a value is zero */
  function evalCurve(c, E) {
    if (!c) return 0;
    var n = c.n, e = c.e;
    if (E <= e[0]) return c.s[0];
    if (E >= e[n - 1]) return c.s[n - 1];
    var lo = 0, hi = n - 1;
    while (hi - lo > 1) { var m = (lo + hi) >> 1; if (e[m] <= E) lo = m; else hi = m; }
    var s0 = c.s[lo], s1 = c.s[hi];
    if (s0 > 0 && s1 > 0) return Math.exp(c.ls[lo] + (c.ls[hi] - c.ls[lo]) * (Math.log(E) - c.le[lo]) / (c.le[hi] - c.le[lo]));
    return s0 + (s1 - s0) * (E - e[lo]) / (e[hi] - e[lo]);
  }


  /* ---------- thermal scattering law S(alpha,beta) (ENDF-6 MF7/MT4 layout, symmetric form, LASYM=0) ----------
   * sigma(E -> E', mu) = sigma_b/(2kT) sqrt(E'/E) exp(-beta/2) S(alpha,|beta|),
   *   beta = (E'-E)/kT,  alpha = (E + E' - 2 mu sqrt(E E')) / (A0 kT).
   * Integrating over mu (-> alpha):  dsigma/dE' = sigma_b A0/(4E) exp(-beta/2) * int S(alpha,|beta|) dalpha.
   * ln S is interpolated linearly in alpha (ENDF INT=4). Energies inside are in eV, outside (API) in MeV. */
  function b64bytes(str) {
    if (typeof atob === 'function') { var b = atob(str), n = b.length, u = new Uint8Array(n); for (var i = 0; i < n; i++) u[i] = b.charCodeAt(i); return u; }
    return new Uint8Array(Buffer.from(str, 'base64'));
  }
  var thermalCache = null;
  function makeThermal(d) {
    if (thermalCache && thermalCache.src === d) return thermalCache.obj;
    var nb = d.nb, na = d.na, beta = d.beta, al = d.alpha, kT = d.kT, A0 = d.A0, k, j;
    var q = b64bytes(d.q), lnS = new Float64Array(nb * na);
    for (k = 0; k < nb * na; k++) { var qq = q[2 * k] | (q[2 * k + 1] << 8); lnS[k] = qq === 0 ? -300 : qq * d.lnStep + d.lnOff; }
    var sigb = d.sigFree * Math.pow((A0 + 1) / A0, 2), amax = al[na - 1];
    var pw = new Float64Array(nb);                                   // power-law exponent below the first alpha point
    for (k = 0; k < nb; k++) pw[k] = Math.max(-1, Math.min(0, (lnS[k * na + 1] - lnS[k * na]) / Math.log(al[1] / al[0])));

    function rowInteg(kb, a1, a2) {
      if (!(a2 > a1) || a1 >= amax) return 0;
      if (a2 > amax) a2 = amax;
      var o = kb * na, tot = 0, p = pw[kb];
      if (a1 < al[0]) {
        var b = Math.min(a2, al[0]), s0 = Math.exp(lnS[o]);
        tot += p === -1 ? s0 * al[0] * Math.log(b / a1) : s0 * Math.pow(al[0], -p) * (Math.pow(b, p + 1) - Math.pow(a1, p + 1)) / (p + 1);
        a1 = b; if (a2 <= a1) return tot;
      }
      var i = 0; while (i < na - 2 && al[i + 1] <= a1) i++;
      for (; i < na - 1; i++) {
        var x0 = al[i], x1 = al[i + 1]; if (x0 >= a2) break;
        var lo = Math.max(a1, x0), hi = Math.min(a2, x1); if (hi <= lo) continue;
        var l0 = lnS[o + i], m = (lnS[o + i + 1] - l0) / (x1 - x0);
        tot += Math.abs(m) < 1e-12 ? Math.exp(l0 + m * (lo - x0)) * (hi - lo) : (Math.exp(l0 + m * (hi - x0)) - Math.exp(l0 + m * (lo - x0))) / m;
      }
      return tot;
    }
    function rowSample(kb, a1, a2, u) {                              // alpha with density S(alpha,beta_k) on [a1,a2]; null if empty
      var I = rowInteg(kb, a1, a2); if (!(I > 0)) return null;
      var target = u * I, o = kb * na, p = pw[kb];
      if (a2 > amax) a2 = amax;
      if (a1 < al[0]) {
        var b = Math.min(a2, al[0]), s0 = Math.exp(lnS[o]);
        var part = p === -1 ? s0 * al[0] * Math.log(b / a1) : s0 * Math.pow(al[0], -p) * (Math.pow(b, p + 1) - Math.pow(a1, p + 1)) / (p + 1);
        if (target <= part) {
          if (p === -1) return a1 * Math.exp(target / (s0 * al[0]));
          return Math.pow(Math.pow(a1, p + 1) + target * (p + 1) / (s0 * Math.pow(al[0], -p)), 1 / (p + 1));
        }
        target -= part; a1 = b; if (a2 <= a1) return a2;
      }
      var i = 0; while (i < na - 2 && al[i + 1] <= a1) i++;
      for (; i < na - 1; i++) {
        var x0 = al[i], x1 = al[i + 1]; if (x0 >= a2) break;
        var lo = Math.max(a1, x0), hi = Math.min(a2, x1); if (hi <= lo) continue;
        var l0 = lnS[o + i], m = (lnS[o + i + 1] - l0) / (x1 - x0), f0 = Math.exp(l0 + m * (lo - x0));
        var seg = Math.abs(m) < 1e-12 ? f0 * (hi - lo) : (Math.exp(l0 + m * (hi - x0)) - f0) / m;
        if (target <= seg) return Math.abs(m) < 1e-12 ? lo + target / f0 : lo + Math.log(1 + m * target / f0) / m;
        target -= seg;
      }
      return a2;
    }

    /* tables of the outgoing-energy distribution on a log grid of incident energies */
    var Emax = d.Emax, E0g = 1e-5, nE = 110, grid = [], i;
    for (i = 0; i < nE; i++) grid.push(E0g * Math.pow(Emax / E0g, i / (nE - 1)));
    var dists = grid.map(function (E) {
      var list = [];
      for (var kb = 0; kb < nb; kb++) {
        var b = beta[kb], signs = b === 0 ? [1] : [1, -1];
        for (var si = 0; si < signs.length; si++) {
          var bt = signs[si] * b, e2 = E + bt * kT; if (e2 <= 0) continue;
          var sq = Math.sqrt(E * e2), I = rowInteg(kb, (E + e2 - 2 * sq) / (A0 * kT), (E + e2 + 2 * sq) / (A0 * kT));
          list.push({ e: e2, k: kb, w: I * Math.exp(-bt / 2) });
        }
      }
      list.sort(function (a, b2) { return a.e - b2.e; });
      var n = list.length, cum = new Float64Array(n), ek = new Float64Array(n), kk = new Int16Array(n), lo = new Float64Array(n), hi = new Float64Array(n), c = 0;
      for (var jj = 0; jj < n; jj++) {
        lo[jj] = jj === 0 ? list[0].e : 0.5 * (list[jj - 1].e + list[jj].e);
        hi[jj] = jj === n - 1 ? list[n - 1].e : 0.5 * (list[jj].e + list[jj + 1].e);
        c += list[jj].w * (hi[jj] - lo[jj]); cum[jj] = c; ek[jj] = list[jj].e; kk[jj] = list[jj].k;
      }
      return { cum: cum, lo: lo, hi: hi, kk: kk, n: n, sigma: sigb * A0 / (4 * E) * c };
    });
    var lg = grid.map(Math.log), ls = dists.map(function (x) { return Math.log(x.sigma); }), step = lg[1] - lg[0];

    var obj = {
      EmaxMeV: Emax * 1e-6, kT: kT, A0: A0, sigb: sigb,
      /* thermal (inelastic incoherent) cross section per H atom [b] at E [MeV], log-log interpolation on the grid */
      sigma: function (EMeV) {
        var E = EMeV * 1e6; if (E <= grid[0]) return dists[0].sigma * Math.sqrt(grid[0] / E);   // ~1/v below the grid
        if (E >= Emax) return dists[nE - 1].sigma;
        var x = (Math.log(E) - lg[0]) / step, ii = Math.min(nE - 2, Math.floor(x)), f = x - ii;
        return Math.exp(ls[ii] * (1 - f) + ls[ii + 1] * f);
      },
      sample: function (EMeV, rng) {
        var E = EMeV * 1e6, x = E <= grid[0] ? 0 : (Math.log(E) - lg[0]) / step, ii = Math.min(nE - 1, Math.floor(x));
        if (ii < nE - 1 && rng() < x - ii) ii++;                   // pick one of the two neighbouring grid energies
        var D = dists[ii], e2 = 0, kb = 0, tries = 0;
        for (;;) {
          var u = rng() * D.cum[D.n - 1], a = 0, b = D.n - 1;
          while (a < b) { var m = (a + b) >> 1; if (D.cum[m] < u) a = m + 1; else b = m; }
          kb = D.kk[a];
          var Ep = D.lo[a] + (D.hi[a] - D.lo[a]) * rng();          // energy transfer taken from the grid energy: beta = (E'-E_i)/kT
          e2 = E + (Ep - grid[ii]);
          if (e2 > 1e-9 || ++tries > 50) break;
        }
        if (!(e2 > 1e-9)) e2 = E;
        var sq = Math.sqrt(E * e2), a1 = (E + e2 - 2 * sq) / (A0 * kT), a2 = (E + e2 + 2 * sq) / (A0 * kT);
        var al_ = rowSample(kb, a1, a2, rng()), mu;
        if (al_ === null) mu = 2 * rng() - 1;
        else mu = Math.max(-1, Math.min(1, (E + e2 - al_ * A0 * kT) / (2 * sq)));
        return { Eout: e2 * 1e-6, mu: mu, beta: (e2 - E) / kT, alpha: al_ === null ? NaN : al_ };
      }
    };
    thermalCache = { src: d, obj: obj };
    return obj;
  }

  function NeutronMaterial(key, rho, name, nxs, opts) {
    this.particle = 'neutron';
    this.key = key; this.rho = rho; this.name = name;
    this.tsl = opts && opts.tsl && key === 'water' ? makeThermal(opts.tsl) : null;
    var recipe = N_RECIPES[key];
    if (!recipe) throw new Error('no neutron recipe for ' + key);
    var den = 0, i, missing = [];
    for (i = 0; i < recipe.length; i++) {
      if (!nxs || !nxs.nuclides || !nxs.nuclides[recipe[i].id]) missing.push(recipe[i].id);
      den += recipe[i].n * recipe[i].M;
    }
    if (missing.length) { var er = new Error('missing neutron data: ' + missing.join(', ')); er.missing = missing; throw er; }
    this.nuc = [];
    for (i = 0; i < recipe.length; i++) {
      var d = nxs.nuclides[recipe[i].id], A = d.m / M_NEUTRON, thr = 0, c;
      var inel = prepCurve(d.inel);
      if (inel) {                                                 // inelastic threshold = last zero before the first positive value
        for (var k = 0; k < inel.n; k++) if (inel.s[k] > 0) { thr = k > 0 ? inel.e[k - 1] : inel.e[0]; break; }
      }
      this.nuc.push({
        id: recipe[i].id, label: N_LABEL[recipe[i].id] || recipe[i].id, A: A,
        N: rho * NA_B * recipe[i].n / den,                        // atoms / (b cm)
        c: [prepCurve(d.el), inel, prepCurve(d.cap), prepCurve(d.oth)],
        Qinel: thr * A / (A + 1),                                 // lab-energy loss of the first excited level (approx.)
        tsl: this.tsl && recipe[i].id === 'h1' ? this.tsl : null  // S(alpha,beta) below its Emax replaces the elastic cross section
      });
    }
    this.single = this.nuc.length === 1;
  }
  /* macroscopic cross sections [1/cm]: totals by reaction + the same for every nuclide */
  NeutronMaterial.prototype.macro = function (E) {
    var nn = this.nuc.length, arr = [0, 0, 0, 0], nt = new Array(nn), na = new Array(nn), i, r;
    for (i = 0; i < nn; i++) {
      var nu = this.nuc[i], a = [0, 0, 0, 0], t = 0;
      for (r = 0; r < 4; r++) {
        a[r] = (r === 0 && nu.tsl && E < nu.tsl.EmaxMeV ? nu.tsl.sigma(E) : evalCurve(nu.c[r], E)) * nu.N;
        t += a[r]; arr[r] += a[r];
      }
      nt[i] = t; na[i] = a;
    }
    return { arr: arr, tot: arr[0] + arr[1] + arr[2] + arr[3], el: arr[0], inel: arr[1], cap: arr[2], oth: arr[3], nucTot: nt, nucArr: na };
  };
  /* log-spaced curves for plotting (macroscopic, 1/cm) */
  NeutronMaterial.prototype.curves = function () {
    if (this._curves) return this._curves;
    var o = { E: [], el: [], inel: [], cap: [], oth: [], tot: [] }, n = 700;
    for (var k = 0; k <= n; k++) {
      var E = NE_MIN * Math.pow(NE_MAX / NE_MIN, k / n), m = this.macro(E);
      o.E.push(E); o.el.push(m.arr[0]); o.inel.push(m.arr[1]); o.cap.push(m.arr[2]); o.oth.push(m.arr[3]); o.tot.push(m.tot);
    }
    this._curves = o; return o;
  };

  function pickIndex(xi, vals, total) {                           // cumulative-vector sampling
    var x = xi * total, c = 0, i;
    for (i = 0; i < vals.length - 1; i++) { c += vals[i]; if (x < c) return i; }
    return vals.length - 1;
  }

  /* elastic scattering on a nucleus of mass ratio A at rest; isotropic in the centre of mass */
  function sampleElastic(E, A, rng, noFloor) {
    var muc = 2 * rng() - 1, D = A * A + 1 + 2 * A * muc;
    var Eout = E * D / ((A + 1) * (A + 1)), mu = (1 + A * muc) / Math.sqrt(D);
    var floored = false;
    if (!noFloor && Eout < E_THERMAL) { Eout = E_THERMAL; floored = true; } // pseudo-thermalisation: no energy below 0.0253 eV
    return { Eout: Eout, mu: mu, muc: muc, floored: floored };
  }
  /* inelastic: evaporation spectrum f(E') ~ E' exp(-E'/theta), E' < E - Q, isotropic emission */
  function sampleInelastic(E, nu, rng) {
    var avail = Math.max(E - nu.Qinel, 1e-9), theta = Math.max(0.02, Math.sqrt(8 * avail / nu.A)), Eout = 0, tries = 0;
    for (;;) {
      tries++;
      Eout = -theta * Math.log(Math.max(rng() * rng(), 1e-300));
      if (Eout < avail) break;
      if (tries > 200) { Eout = avail * rng(); break; }
    }
    return { Eout: Math.max(Eout, E_THERMAL), mu: 2 * rng() - 1, theta: theta, avail: avail, tries: tries };
  }

  function interactN(sim, ph, sig, ni, r) {
    var rng = sim.rng, nu = sim.mat.nuc[ni], info = { r: r, nuc: ni, label: nu.label, absorbed: false, Ein: ph.E };
    ph.nc++;
    sim.tally.byReaction[r]++;
    if (r >= 2) {                                                 // capture / other: absorbed
      info.absorbed = true;
      if (ph.path) ph.path.push([ph.x, ph.y, r, ph.z]);
      return info;
    }
    var th = r === 0 && nu.tsl && ph.E < nu.tsl.EmaxMeV, o = th ? nu.tsl.sample(ph.E, rng) : r === 0 ? sampleElastic(ph.E, nu.A, rng, !!sim.mat.tsl) : sampleInelastic(ph.E, nu, rng);
    info.thermal = th;
    var phi = 2 * Math.PI * rng(), d = rotate(ph.u, ph.v, ph.w, o.mu, phi);
    ph.u = d[0]; ph.v = d[1]; ph.w = d[2]; ph.E = o.Eout;
    info.Eout = o.Eout; info.mu = o.mu; info.phi = phi; info.dir = d; info.detail = o;
    if (ph.path) ph.path.push([ph.x, ph.y, r, ph.z]);
    return info;
  }

  var N_MAX_COLLISIONS = 200000;
  function runHistoryN(sim) {
    var ph = newPhoton(sim.E0, sim), mat = sim.mat, rng = sim.rng, T = sim.T, sph = sim.sphere;
    if (sim.record) ph.path = [[ph.x, ph.y, 'S', ph.z]];
    for (;;) {
      if (ph.nc > N_MAX_COLLISIONS) { finish(sim, ph, 'C'); return { out: 'C', ph: ph }; }
      var sig = mat.macro(ph.E);
      var s = sampleFlight(rng(), sig.tot);
      var out = sph ? transportS(sim, ph, s) : transport(ph, s, T);
      if (out) { finish(sim, ph, out); return { out: out, ph: ph }; }
      var ni = mat.single ? 0 : pickIndex(rng(), sig.nucTot, sig.tot);
      var r = pickIndex(rng(), sig.nucArr[ni], sig.nucTot[ni]);
      var info = interactN(sim, ph, sig, ni, r);
      if (info.absorbed) { sim.tally.absBy[r]++; finish(sim, ph, 'A'); return { out: 'A', ph: ph }; }
    }
  }
  function* historyStepsN(sim) {
    var ph = newPhoton(sim.E0, sim), mat = sim.mat, rng = sim.rng, T = sim.T, sph = sim.sphere;
    ph.path = [[ph.x, ph.y, 'S', ph.z]];
    yield { type: 'birth', ph: ph, n: sim.history + 1, srcLen: sim.srcLen, srcAxis: sim.srcAxis, iso: sph, mu: ph.mu0, phi: ph.phi0 };
    for (;;) {
      if (ph.nc > N_MAX_COLLISIONS) { finish(sim, ph, 'C'); yield { type: 'end', out: 'C', ph: ph }; return; }
      var sig = mat.macro(ph.E);
      yield { type: 'xs_total', ph: ph, sig: sig };
      var xi = rng(), s = sampleFlight(xi, sig.tot);
      var db = sph ? distSphere(ph, sim.R) : ph.u > 0 ? (T - ph.x) / ph.u : (0 - ph.x) / ph.u;
      yield { type: 'flight', ph: ph, sig: sig, xi: xi, s: s, db: db };
      var from = [ph.x, ph.y];
      var out = sph ? transportS(sim, ph, s) : transport(ph, s, T);
      yield { type: 'move', ph: ph, from: from, s: s, db: db, out: out };
      if (out) { finish(sim, ph, out); yield { type: 'end', out: out, ph: ph }; return; }
      yield { type: 'xs_partial', ph: ph, sig: sig };
      var ni = 0;
      if (!mat.single) {
        var xn = rng(); ni = pickIndex(xn, sig.nucTot, sig.tot);
        yield { type: 'nuclide_sample', ph: ph, sig: sig, xi: xn, ni: ni };
      }
      var xr = rng(), r = pickIndex(xr, sig.nucArr[ni], sig.nucTot[ni]);
      yield { type: 'reaction_sample', ph: ph, sig: sig, xi: xr, r: r, ni: ni };
      var info = interactN(sim, ph, sig, ni, r);
      yield { type: 'interaction', ph: ph, info: info };
      if (info.absorbed) { sim.tally.absBy[r]++; finish(sim, ph, 'A'); yield { type: 'end', out: 'A', r: r, ph: ph }; return; }
    }
  }

  function makeNeutronMaterial(key, rho, name, nxs, opts) { return new NeutronMaterial(key, rho, name, nxs, opts); }

  /* ---------- derived results ---------- */
  /* probability of crossing the rest of the slab without interaction, averaged over the source */
  function analyticUncollided(sim, mu) {
    var T = sim.T, x0 = sim.srcX0 || 0, L = sim.srcLen || 0;
    if (!(L > 0)) return Math.exp(-mu * T);
    if (sim.srcAxis === 'x') return (Math.exp(-mu * (T - x0 - L)) - Math.exp(-mu * (T - x0))) / (mu * L);
    return Math.exp(-mu * (T - x0));
  }

  /* volume-averaged flux spectrum in the sphere per unit lethargy, E phi(E) [per source neutron and cm^2], with 1-sigma errors */
  function fluxSpectrum(sim) {
    var t = sim.tally, N = t.N, nb = t.nbins, V = 4 / 3 * Math.PI * Math.pow(sim.R, 3), dl = (t.lmax - t.lmin) / nb, out = [], i;
    for (i = 0; i < nb; i++) {
      var m = N ? t.fl[i] / N : 0, v = N > 1 ? Math.max(0, t.fl2[i] / N - m * m) / (N - 1) : 0;
      out.push({ e0: Math.exp(t.lmin + i * dl), e1: Math.exp(t.lmin + (i + 1) * dl), y: m / (V * dl), err: Math.sqrt(v) / (V * dl) });
    }
    return out;
  }
  function summarySphere(sim, mu) {
    var t = sim.tally, N = t.N, V = 4 / 3 * Math.PI * Math.pow(sim.R, 3);
    function frac(n) { var p = N ? n / N : 0; return { p: p, err: N && p > 0 ? Math.sqrt((1 - p) / (N * p)) : 0, n: n }; }
    var fm = N ? t.trk / N : 0, fv = N > 1 ? Math.max(0, t.trk2 / N - fm * fm) / (N - 1) : 0;
    return {
      N: N, sigT0: mu, mfp: 1 / mu, thickness_mfp: sim.R * mu, V: V,
      leaked: frac(t.leaked), uncollided: frac(t.uncollided), scattered: frac(t.leaked - t.uncollided), absorbed: frac(t.absorbed),
      analyticUncollided: Math.exp(-mu * sim.R), collMean: N ? t.coll / N : 0,
      flux: fm / V, fluxErr: Math.sqrt(fv) / V
    };
  }

  function summary(sim) {
    var t = sim.tally, N = t.N, mu = sim.mat.macro(sim.E0).tot;
    if (sim.sphere) return summarySphere(sim, mu);
    function frac(n) {
      var p = N ? n / N : 0;
      return { p: p, err: N && p > 0 ? Math.sqrt((1 - p) / (N * p)) : 0, n: n };   // err = relative error
    }
    var fm = N ? t.flux / N : 0;
    var fv = N > 1 ? Math.max(0, t.flux2 / N - fm * fm) / (N - 1) : 0;
    return {
      N: N, sigT0: mu, mfp: 1 / mu, thickness_mfp: sim.T * mu,
      transmitted: frac(t.transmitted), uncollided: frac(t.uncollided),
      scattered: frac(t.transmitted - t.uncollided), reflected: frac(t.reflected),
      absorbed: frac(t.absorbed),
      analyticUncollided: analyticUncollided(sim, mu),
      flux: fm, fluxErr: Math.sqrt(fv)
    };
  }

  var api = {
    ME: ME, REACTIONS: REACTIONS, ZEFF: ZEFF, E_CUTOFF: E_CUTOFF,
    makeRNG: makeRNG, Material: Material, makeSim: makeSim, makeTally: makeTally, configureSource: configureSource,
    sampleFlight: sampleFlight, sampleReaction: sampleReaction,
    sampleCompton: sampleCompton, sampleRayleigh: sampleRayleigh, rotate: rotate,
    kleinNishinaPDF: kleinNishinaPDF, rayleighPDF: rayleighPDF,
    runHistory: runHistory, historySteps: historySteps, summary: summary, fluxSpectrum: fluxSpectrum, distSphere: distSphere,
    makeNeutronMaterial: makeNeutronMaterial, N_RECIPES: N_RECIPES, N_REACTIONS: N_REACTIONS, NE_MIN: NE_MIN, NE_MAX: NE_MAX,
    makeThermal: makeThermal, E_THERMAL: E_THERMAL, sampleElastic: sampleElastic, sampleInelastic: sampleInelastic, evalCurve: evalCurve, prepCurve: prepCurve
  };
  root.MCEngine = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
