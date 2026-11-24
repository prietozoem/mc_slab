/* UI for the Monte Carlo transport examples: 1 = slab shielding (photons / neutrons), 2 = neutron point source in a sphere. */
(function () {
  'use strict';
  var MC = window.MCEngine, I18 = window.I18N;
  function t(k, v) {                                                    // sphere variant (.s), else neutron variant (.n), else the base string
    var sk = k + '.s', nk = k + '.n';
    if (S.example === 'sphere' && I18.has(sk)) return I18.t(sk, v);
    return I18.t(S.particle === 'neutron' && I18.has(nk) ? nk : k, v);
  }
  var $ = function (id) { return document.getElementById(id); };

  var PH_COL = ['#7e57c2', '#1e88e5', '#e53935', '#00838f'], N_COL = ['#1e88e5', '#8e24aa', '#e53935', '#00838f'];
  var NUC_COL = ['#1e88e5', '#e53935', '#43a047', '#fb8c00', '#8e24aa', '#00838f', '#6d4c41', '#546e7a', '#c0ca33', '#d81b60', '#3949ab', '#00897b', '#7cb342'];
  var RX_COL = PH_COL, RX_KEY = ['coh', 'inc', 'pe', 'pair'];
  var RX_SHAPE = ['circle', 'square', 'triangle', 'diamond'];
  function setTables() {
    var n = S.particle === 'neutron';
    RX_COL = n ? N_COL : PH_COL; RX_KEY = n ? ['el', 'inel', 'cap', 'oth'] : ['coh', 'inc', 'pe', 'pair'];
    $('geomLegend').dataset.lang = '';
  }
  var OUT_COL = { T: '#2e7d32', R: '#ef6c00', L: '#2e7d32', A: '#6b7280', C: '#6b7280' };
  var MAT_ORDER = ['water', 'aluminum', 'iron', 'concrete', 'tungsten', 'lead'];

  var S = {
    example: 'slab', fluxLog: false, params: { slab: null, sphere: null },
    particle: 'photon', mat: null, sim: null, E0: 1, T: 2, nHist: 20000, seed: 12345, mode: 'auto', view: '2d',
    running: false, raf: 0, timer: 0, gen: null, ev: null, current: null, finished: [],
    winRender: null
  };
  var winPlot = new Plot($('winCanvas'));
  var specPlot = new Plot($('spec'));
  var fluxPlot = new Plot($('fluxSpec'));

  /* ---------------- helpers ---------------- */
  function fmt(v, p) {
    if (v === 0) return '0';
    var a = Math.abs(v); p = p || 4;
    if (a >= 1e5 || a < 1e-3) return v.toExponential(p - 1).replace('e+', 'e');
    return String(Number(v.toPrecision(p)));
  }
  function pct(p) { return (100 * p).toFixed(p < 0.001 && p > 0 ? 4 : 2) + ' %'; }
  function lang() { return window.I18N.getLang(); }
  function matName() { return S.mat.name[lang()] || S.mat.name.en; }
  function clamp(v, lo, hi) { return Math.min(hi, Math.max(lo, v)); }

  function sizeCanvas(cv) {
    var dpr = window.devicePixelRatio || 1, W = cv.clientWidth, H = cv.clientHeight;
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    var g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
    return { g: g, W: W, H: H };
  }

  /* ---------------- setup / reset ---------------- */
  function readParams() {
    var neu = S.particle === 'neutron', sph = S.example === 'sphere';
    S.E0 = neu ? clamp(parseFloat($('energy').value) || 2, 1e-11, 20) : clamp(parseFloat($('energy').value) || 1, 0.001, 100);
    S.T = Math.max(0.001, parseFloat($('thickness').value) || 1);
    S.nHist = Math.max(1, Math.floor(parseFloat($('nhist').value) || 1000));
    S.seed = Math.max(1, Math.floor(parseFloat($('seed').value) || 1));
    S.srcType = sph ? 'beam' : $('srctype').value;
    S.srcAxis = 'y';                                   // the line lies on the entrance face; code axis y is displayed as X
    S.srcX0 = 0;
    S.srcLen = S.srcType === 'line' ? Math.max(0.001, parseFloat($('srclen').value) || 1) : 0;
    ['srclenRow'].forEach(function (id) { $(id).hidden = S.srcType !== 'line'; });
    $('energy').value = S.E0; $('thickness').value = S.T;
    var key = $('material').value;
    $('tslRow').hidden = !(neu && key === 'water');
    if (neu) {
      var wantTsl = key === 'water' && $('tslChk').checked;
      if (wantTsl && !window.TSL_H2O) loadTsl();
      try { S.mat = MC.makeNeutronMaterial(key, window.XCOM[key].rho, window.XCOM[key].name, window.NXS, { tsl: wantTsl ? window.TSL_H2O : null }); if (!S.tslLoading) showMsg(''); }
      catch (e) {
        showMsg(t('msg.missing', { m: (e.missing || [e.message]).join(', ') }));
        S.particle = 'photon'; S.example = 'slab'; S.E0 = 1; $('energy').value = 1; syncExampleUI(); syncParticleUI(); setTables();
        S.mat = new MC.Material(key, window.XCOM[key]);
      }
    } else S.mat = new MC.Material(key, window.XCOM[key]);
  }
  function loadTsl() {                                                   // lazy-load the S(alpha,beta) table, then rebuild the simulation
    if (S.tslLoading) return;
    S.tslLoading = true;
    var sc = document.createElement('script'); sc.src = 'js/data_tsl.js';
    sc.onload = function () { S.tslLoading = false; if (window.TSL_H2O) resetSim(); else tslFail(); };
    sc.onerror = tslFail; document.head.appendChild(sc);
  }
  function tslFail() { S.tslLoading = false; $('tslChk').checked = false; showMsg(t('msg.notsl')); resetSim(); }
  function showMsg(h) { $('dataMsg').hidden = !h; $('dataMsg').innerHTML = h || ''; }
  function syncParticleUI() {
    var neu = S.particle === 'neutron';
    document.querySelectorAll('#particleToggle button').forEach(function (b) { b.classList.toggle('on', b.dataset.particle === S.particle); });
    $('energy').min = neu ? '1e-11' : '0.001'; $('energy').max = neu ? '20' : '100';
    document.title = t('app.title');
    document.querySelectorAll('[data-i18n]').forEach(function (el) { el.innerHTML = t(el.getAttribute('data-i18n')); });
  }
  function loadNeutronData(done) {
    if (window.NXS) return done(true);
    var sc = document.createElement('script'); sc.src = 'js/data_n.js';
    sc.onload = function () { done(!!window.NXS); }; sc.onerror = function () { done(false); };
    document.head.appendChild(sc);
  }
  function setParticle(p) {
    if (p === S.particle) return;
    function go() {
      S.particle = p; setTables();
      $('energy').value = p === 'neutron' ? 2 : 1;
      showMsg(''); syncParticleUI(); applyLang(); resetSim();
    }
    if (p === 'neutron') loadNeutronData(function (ok) { if (ok) go(); else showMsg(t('msg.nodata')); });
    else go();
  }

  /* ---------------- example selector (1 = slab, 2 = point source in a sphere) ---------------- */
  var EX_DEFAULT = {
    slab: { particle: 'photon', material: 'lead', energy: 1, thickness: 4, nhist: 20000, tsl: true },
    sphere: { particle: 'neutron', material: 'water', energy: 2, thickness: 30, nhist: 5000, tsl: true }
  };
  function grabInputs() {
    return { particle: S.particle, material: $('material').value, energy: $('energy').value, thickness: $('thickness').value, nhist: $('nhist').value, tsl: $('tslChk').checked };
  }
  function applyInputs(p) {
    $('material').value = p.material; $('energy').value = p.energy; $('thickness').value = p.thickness; $('nhist').value = p.nhist; $('tslChk').checked = !!p.tsl;
  }
  function syncExampleUI() {
    var sph = S.example === 'sphere';
    document.querySelectorAll('#exToggle button').forEach(function (b) { b.classList.toggle('on', b.dataset.ex === S.example); });
    $('particleLabel').hidden = $('particleToggle').hidden = sph;      // the sphere example is neutrons only
    $('srctypeRow').hidden = sph;                                      // the source is always an isotropic point at the centre
    $('fluxPanel').hidden = !sph;
    $('specTitle').hidden = $('spec').hidden = sph;                    // example 2 has no slab-style output spectrum
    $('geomLegend').dataset.lang = '';
  }
  function setExample(ex) {
    if (ex === S.example) return;
    var target = S.params[ex] || EX_DEFAULT[ex];
    function go() {
      stopAll();
      S.params[S.example] = grabInputs();
      S.example = ex; S.sim = null; S.particle = target.particle; setTables();
      cam2.zoom = 1; cam2.px = 0; cam2.py = 0;
      showMsg(''); syncExampleUI(); syncParticleUI(); applyLang(); applyInputs(target); resetSim();
      try { history.replaceState(null, '', ex === 'sphere' ? '#sphere' : location.pathname + location.search); } catch (e) { /* file:// or sandboxed */ }
    }
    if (target.particle === 'neutron') loadNeutronData(function (ok) { if (ok) go(); else showMsg(t('msg.nodata')); });
    else go();
  }

  function stopAll() {
    S.running = false; cancelAnimationFrame(S.raf);
    clearInterval(S.timer); S.timer = 0;
    $('btnCont').classList.remove('active'); $('btnCont').textContent = t('btn.cont');
    $('btnRun').disabled = false; $('btnPause').disabled = true;
  }

  function newSim() {                                              // a fresh simulation with the current parameters (same seed -> same random sequence)
    var sph = S.example === 'sphere';
    var sim = MC.makeSim({ material: S.mat, E0: S.E0, T: S.T, R: S.T, geom: sph ? 'sphere' : 'slab', seed: S.seed, nbins: S.particle === 'neutron' ? 200 : 100 });
    if (!sph) MC.configureSource(sim, S.srcLen, S.srcAxis, S.srcX0);
    return sim;
  }
  function resetSim() {
    stopAll(); readParams();
    S.sim = newSim();
    if (S.example !== 'sphere') { S.srcLen = S.sim.srcLen; S.srcX0 = S.sim.srcX0; }   // clamped values
    
    S.gen = null; S.ev = null; S.current = null; S.finished = []; S.winRender = null; S.actions = [];
    updatePrevBtn();
    if (S.mode === 'step') showReady(); else $('win').hidden = true;
    updateInfo(); drawGeom(); updateTallies();
  }

  function updateInfo() {
    var sig = S.mat.macro(S.E0).tot;
    $('mfpInfo').innerHTML = t('ctl.mfp', { v: fmt(sig * S.T, 3), mu: fmt(sig, 3) });
  }

  /* ---------------- geometry drawing ---------------- */
  function srcSegment() {                                         // end points [x,y,z] of the line source
    var x0 = S.srcX0, h = S.srcLen / 2;
    if (S.srcAxis === 'x') return { a: [x0, 0, 0], b: [x0 + S.srcLen, 0, 0] };
    if (S.srcAxis === 'y') return { a: [x0, -h, 0], b: [x0, h, 0] };
    return { a: [x0, 0, -h], b: [x0, 0, h] };
  }

  function view2(W, H) {                                          // scale [px/cm] at zoom 1 and screen x of the origin
    return S.example === 'sphere' ? { base: 0.41 * Math.min(W, H) / S.T, x0: W / 2 } : { base: W * 0.38 / S.T, x0: W * 0.31 };
  }
  function drawGeom2D(proposal) {
    var cv = $('geom'), c = sizeCanvas(cv), g = c.g, W = c.W, H = c.H;
    if (!W || !H || !S.sim) return;
    var T = S.T, sph = S.example === 'sphere', vp = view2(W, H), sx = vp.base * cam2.zoom, x0 = vp.x0 + cam2.px, zv = S.view === '2dz';       // zv: vertical axis is z instead of y
    function X(x) { return x0 + x * sx; }
    function Y(y) { return H / 2 + cam2.py - y * sx; }
    function PV(pt) { return zv ? (pt[3] || 0) : pt[1]; }
    if (sph) sphereBackdrop(); else {
    // slab
    g.fillStyle = '#e8edf3'; g.fillRect(X(0), 0, T * sx, H);
    g.strokeStyle = '#8a97a8'; g.lineWidth = 1.5; g.setLineDash([]);
    g.beginPath(); g.moveTo(X(0), 0); g.lineTo(X(0), H); g.moveTo(X(T), 0); g.lineTo(X(T), H); g.stroke();
    g.fillStyle = '#4b5563'; g.font = '12px system-ui, sans-serif'; g.textAlign = 'center';
    g.fillText(t('geom.slab') + ': ' + matName(), X(T / 2), 16);
    g.fillText(t('geom.scale', { T: fmt(T, 3) }), X(T / 2), H - 8);
    g.textAlign = 'right'; g.fillText(t('geom.entrance'), X(0) - 6, H - 8);
    g.textAlign = 'left'; g.fillText(t('geom.exit'), X(T) + 6, H - 8);
    // beam arrow
    function arrowAt(yy) {
      g.strokeStyle = '#111'; g.fillStyle = '#111'; g.lineWidth = 2;
      g.beginPath(); g.moveTo(6, yy); g.lineTo(X(0) - 4, yy); g.stroke();
      g.beginPath(); g.moveTo(X(0), yy); g.lineTo(X(0) - 10, yy - 5); g.lineTo(X(0) - 10, yy + 5); g.closePath(); g.fill();
    }
    g.fillStyle = '#8a97a8'; g.textAlign = 'right'; g.fillText(t('geom.hint2d'), W - 8, 34); g.textAlign = 'left';
    if (S.srcLen > 0) {
      var sg = srcSegment(), xa = sg.a[0], xb = sg.b[0], va = zv ? sg.a[2] : sg.a[1], vb = zv ? sg.b[2] : sg.b[1], k0;
      for (k0 = 0; k0 < 5; k0++) {                              // emission arrows from sample points of the segment
        var f = (k0 + 0.5) / 5, ax = xa + (xb - xa) * f, av = va + (vb - va) * f;
        if (S.srcAxis !== 'x' && S.srcX0 === 0) arrowAt(Y(av));
        else {
          g.strokeStyle = '#111'; g.fillStyle = '#111'; g.lineWidth = 2;
          g.beginPath(); g.moveTo(X(ax), Y(av)); g.lineTo(X(ax) + 26, Y(av)); g.stroke();
          g.beginPath(); g.moveTo(X(ax) + 34, Y(av)); g.lineTo(X(ax) + 24, Y(av) - 5); g.lineTo(X(ax) + 24, Y(av) + 5); g.closePath(); g.fill();
        }
      }
      g.strokeStyle = '#BA7517'; g.fillStyle = '#BA7517'; g.lineWidth = 6; g.lineCap = 'round';
      g.beginPath(); g.moveTo(X(xa), Y(va)); g.lineTo(X(xb), Y(vb)); g.stroke(); g.lineCap = 'butt';
      if (Math.abs(xb - xa) * sx < 1 && Math.abs(vb - va) * sx < 1) { g.beginPath(); g.arc(X(xa), Y(va), 5, 0, 6.2832); g.fill(); }
      g.fillStyle = '#8a5a00';
      g.fillText(t('geom.srcline', { A: 'X', L: fmt(S.srcLen, 3) }) + ', ' + fmt(S.E0, 3) + ' MeV', 8, 18);
    } else {
      arrowAt(Y(0));
      g.fillStyle = '#111';
      g.fillText(t('geom.beam') + ' ' + fmt(S.E0, 3) + ' MeV', 8, Y(0) - 8);
    }

    }
    function sphereBackdrop() {
      var cx = X(0), cy = Y(0), r = T * sx;
      g.fillStyle = '#e8edf3'; g.strokeStyle = '#8a97a8'; g.lineWidth = 1.5; g.setLineDash([]);
      g.beginPath(); g.arc(cx, cy, r, 0, 6.2832); g.fill(); g.stroke();
      g.fillStyle = '#4b5563'; g.font = '12px system-ui, sans-serif'; g.textAlign = 'center';
      g.fillText(t('geom.sphere') + ': ' + matName(), W / 2, 16);
      g.fillText(t('geom.scale', { T: fmt(T, 3) }), W / 2, H - 8);
      g.setLineDash([4, 3]); g.lineWidth = 1; g.beginPath(); g.moveTo(cx, cy); g.lineTo(cx + r, cy); g.stroke(); g.setLineDash([]);
      g.fillText('R', cx + r / 2, cy - 5);
      g.fillStyle = '#8a97a8'; g.textAlign = 'right'; g.fillText(t('geom.hint2d'), W - 8, 34); g.textAlign = 'left';
    }
    function sphereSource() {
      g.fillStyle = '#BA7517'; g.strokeStyle = '#fff'; g.lineWidth = 1.5;
      g.beginPath(); g.arc(X(0), Y(0), 6, 0, 6.2832); g.fill(); g.stroke();
      g.fillStyle = '#8a5a00'; g.font = '12px system-ui, sans-serif'; g.textAlign = 'left';
      g.fillText(t('geom.src.s') + ', ' + fmt(S.E0, 3) + ' MeV', 8, 18);
    }
    function drawPath(p, out, ux, uy, alpha, markers, width) {
      var pts = p.path; if (!pts || pts.length < 1) return;
      g.globalAlpha = alpha; g.strokeStyle = OUT_COL[out] || '#0b3d91'; g.lineWidth = width;
      g.beginPath();
      for (var i = 0; i < pts.length; i++) { if (i) g.lineTo(X(pts[i][0]), Y(PV(pts[i]))); else g.moveTo(X(pts[i][0]), Y(PV(pts[i]))); }
      g.stroke();
      var last = pts[pts.length - 1];
      if ((out === 'T' || out === 'R' || out === 'L') && ux !== undefined) {           // outgoing ray after leaving the slab / sphere
        var n = Math.hypot(ux, uy);
        if (n > 1e-6) {
          var L = (sph ? 0.2 : 0.35) * T * sx;
          g.setLineDash([3, 3]); g.beginPath(); g.moveTo(X(last[0]), Y(PV(last)));
          g.lineTo(X(last[0]) + ux / n * L, Y(PV(last)) - uy / n * L); g.stroke(); g.setLineDash([]);
        }
      }
      if (markers) for (var k = 1; k < pts.length; k++) {
        var tag = pts[k][2];
        if (typeof tag === 'number') { g.fillStyle = RX_COL[tag]; g.beginPath(); g.arc(X(pts[k][0]), Y(PV(pts[k])), 4, 0, 6.2832); g.fill(); }
      }
      g.globalAlpha = 1;
    }
    var nf = S.finished.length;
    S.finished.forEach(function (p, i) {
      var a = S.mode === 'auto' ? 0.18 + 0.4 * (i / Math.max(1, nf)) : 0.2 + 0.35 * (i / Math.max(1, nf));
      drawPath(p, p.out, p.u, zv ? p.w : p.v, a, S.mode === 'step', 1.2);
    });
    var cur = S.current;
    if (cur) {
      var ph = cur.ph, pts = ph.path, last = pts[pts.length - 1];
      var pyv = zv ? ph.z : ph.y, pvv = zv ? ph.w : ph.v;
      var drawn = { path: pts.slice() };
      if (!cur.out && (Math.abs(last[0] - ph.x) > 1e-12 || Math.abs(PV(last) - pyv) > 1e-12)) drawn.path.push([ph.x, ph.y, 'p', ph.z]);
      drawPath(drawn, cur.out || 'N', ph.u, pvv, 1, true, 2.4);
      if (!cur.out) {                                           // the photon itself
        var px = X(ph.x), py = Y(pyv), n = Math.hypot(ph.u, pvv);
        g.fillStyle = '#ffb300'; g.strokeStyle = '#111'; g.lineWidth = 1.5;
        g.beginPath(); g.arc(px, py, 6, 0, 6.2832); g.fill(); g.stroke();
        if (n > 1e-6) {
          g.strokeStyle = '#111'; g.lineWidth = 2; g.beginPath(); g.moveTo(px, py);
          g.lineTo(px + ph.u / n * 22, py - pvv / n * 22); g.stroke();
        }
      }
      if (proposal && proposal.type === 'flight') {            // proposed flight (before the move)
        var xe = ph.x + ph.u * proposal.s, ye = pyv + pvv * proposal.s;
        g.strokeStyle = '#ffb300'; g.lineWidth = 2.2; g.setLineDash([6, 4]);
        g.beginPath(); g.moveTo(X(ph.x), Y(pyv)); g.lineTo(X(xe), Y(ye)); g.stroke(); g.setLineDash([]);
        g.fillStyle = '#ffb300'; g.beginPath(); g.arc(X(xe), Y(ye), 4, 0, 6.2832); g.fill();
      }
    }
    if (sph) sphereSource();
  }

  function updateLegend() {
    var L = $('geomLegend');
    if (!L.dataset.lang || L.dataset.lang !== lang()) {
      L.dataset.lang = lang();
      L.innerHTML = (S.example === 'sphere' ? ['L', 'A'] : ['T', 'R', 'A']).map(function (k) { return '<span><i style="background:' + OUT_COL[k] + '"></i>' + t('leg.' + k) + '</span>'; }).join('') +
        RX_KEY.map(function (k, i) { return '<span><b style="background:' + RX_COL[i] + '"></b>' + t('rx.' + k + '.s') + '</span>'; }).join('');
    }
  }

  function drawGeom(proposal) {
    if (S.view === '3d') drawGeom3D(proposal); else drawGeom2D(proposal);
    updateLegend();
  }
  function redraw() { drawGeom(S.ev && S.ev.type === 'flight' ? S.ev : undefined); }

  /* ---------------- 3D view (orbit camera, own projection) ---------------- */
  var cam2 = { zoom: 1, px: 0, py: 0 };
  var cam = { yaw: -0.6, pitch: 0.38, zoom: 1 };

  function drawGeom3D(proposal) {
    var cv = $('geom'), c = sizeCanvas(cv), g = c.g, W = c.W, H = c.H;
    if (!W || !H || !S.sim) return;
    var T = S.T, L = T * 2.5, k, sph = S.example === 'sphere', cx0 = sph ? 0 : T / 2;
    function scan(p) {
      var pts = p.path || []; for (var i = 0; i < pts.length; i++) L = Math.max(L, Math.abs(pts[i][1]) * 1.15, Math.abs(pts[i][3] || 0) * 1.15);
    }
    if (sph) L = T; else {
      S.finished.forEach(scan); if (S.current) scan(S.current.ph);
      L = Math.max(L, S.srcLen / 2 * 1.15); L = Math.min(L, Math.max(T * 12, S.srcLen));
    }                                  // half-size of the drawn slab in y, z
    var rad = sph ? 1.25 * T : Math.hypot(1.15 * T, L, L), scale = 0.46 * Math.min(W, H) / rad * cam.zoom, D = rad * 5;
    var cyw = Math.cos(cam.yaw), syw = Math.sin(cam.yaw), cp = Math.cos(cam.pitch), sp = Math.sin(cam.pitch);
    function P(x, y, z) {
      var X = x - cx0, x1 = X * cyw + z * syw, z1 = -X * syw + z * cyw;
      var y2 = y * cp - z1 * sp, z2 = y * sp + z1 * cp, f = D / (D - z2);
      return [W / 2 + x1 * scale * f, H / 2 - y2 * scale * f, z2];
    }
    if (sph) sphere3(); else {
    // slab: 6 faces sorted far -> near
    var C = [];
    [0, T].forEach(function (x) { [-L, L].forEach(function (y) { [-L, L].forEach(function (z) { C.push(P(x, y, z)); }); }); });
    // corner index = ix*4 + iy*2 + iz
    var faces = [
      { q: [0, 1, 3, 2], col: 'rgba(30,136,229,.13)', lab: t('geom.entrance') },      // x = 0
      { q: [4, 5, 7, 6], col: 'rgba(46,125,50,.13)', lab: t('geom.exit') },           // x = T
      { q: [0, 1, 5, 4], col: 'rgba(140,155,175,.14)' }, { q: [2, 3, 7, 6], col: 'rgba(140,155,175,.14)' },
      { q: [0, 2, 6, 4], col: 'rgba(140,155,175,.14)' }, { q: [1, 3, 7, 5], col: 'rgba(140,155,175,.14)' }
    ];
    faces.forEach(function (f) { f.d = (C[f.q[0]][2] + C[f.q[1]][2] + C[f.q[2]][2] + C[f.q[3]][2]) / 4; });
    faces.sort(function (a, b) { return a.d - b.d; });
    g.lineJoin = 'round'; g.font = '12px system-ui, sans-serif';
    faces.forEach(function (f) {
      g.fillStyle = f.col; g.strokeStyle = '#8a97a8'; g.lineWidth = 1.2;
      g.beginPath(); f.q.forEach(function (i, n) { if (n) g.lineTo(C[i][0], C[i][1]); else g.moveTo(C[i][0], C[i][1]); });
      g.closePath(); g.fill(); g.stroke();
    });
    g.fillStyle = '#4b5563'; g.textAlign = 'center';
    g.fillText(t('geom.slab') + ': ' + matName() + '  (T = ' + fmt(T, 3) + ' cm)', W / 2, H - 8);
    g.textAlign = 'left';                                      // face colour key
    g.fillStyle = 'rgba(30,136,229,.45)'; g.fillRect(W - 200, H - 36, 10, 10);
    g.fillStyle = 'rgba(46,125,50,.5)'; g.fillRect(W - 200, H - 20, 10, 10);
    g.fillStyle = '#4b5563'; g.fillText(t('geom.entrance'), W - 185, H - 27); g.fillText(t('geom.exit'), W - 185, H - 11);

    // beam
    function beam3(zz) {
      var b0 = P(-0.9 * T, 0, zz), b1 = P(0, 0, zz);
      g.strokeStyle = '#111'; g.fillStyle = '#111'; g.lineWidth = 2.4;
      g.beginPath(); g.moveTo(b0[0], b0[1]); g.lineTo(b1[0], b1[1]); g.stroke();
      var ang = Math.atan2(b1[1] - b0[1], b1[0] - b0[0]);
      g.beginPath(); g.moveTo(b1[0], b1[1]); g.lineTo(b1[0] - 11 * Math.cos(ang - 0.4), b1[1] - 11 * Math.sin(ang - 0.4));
      g.lineTo(b1[0] - 11 * Math.cos(ang + 0.4), b1[1] - 11 * Math.sin(ang + 0.4)); g.closePath(); g.fill();
      return b0;
    }
    g.textAlign = 'left';
    if (S.srcLen > 0) {
      var sg = srcSegment(), kk, bl;
      function arrow3(p0, p1) {
        var q0 = P(p0[0], p0[1], p0[2]), q1 = P(p1[0], p1[1], p1[2]), an = Math.atan2(q1[1] - q0[1], q1[0] - q0[0]);
        g.strokeStyle = '#111'; g.fillStyle = '#111'; g.lineWidth = 2.4;
        g.beginPath(); g.moveTo(q0[0], q0[1]); g.lineTo(q1[0], q1[1]); g.stroke();
        g.beginPath(); g.moveTo(q1[0], q1[1]); g.lineTo(q1[0] - 10 * Math.cos(an - 0.4), q1[1] - 10 * Math.sin(an - 0.4));
        g.lineTo(q1[0] - 10 * Math.cos(an + 0.4), q1[1] - 10 * Math.sin(an + 0.4)); g.closePath(); g.fill();
      }
      for (kk = 0; kk < 5; kk++) {
        var f = (kk + 0.5) / 5, c = [0, 1, 2].map(function (i) { return sg.a[i] + (sg.b[i] - sg.a[i]) * f; });
        var xs0 = (S.srcAxis !== 'x' && S.srcX0 === 0) ? -0.6 * T : c[0];
        arrow3([xs0, c[1], c[2]], [xs0 === c[0] ? c[0] + 0.3 * T : 0, c[1], c[2]]);
      }
      var s0 = P(sg.a[0], sg.a[1], sg.a[2]), s1 = P(sg.b[0], sg.b[1], sg.b[2]);
      g.strokeStyle = '#BA7517'; g.lineWidth = 6; g.lineCap = 'round';
      g.beginPath(); g.moveTo(s0[0], s0[1]); g.lineTo(s1[0], s1[1]); g.stroke(); g.lineCap = 'butt';
      g.fillStyle = '#8a5a00'; g.fillText(t('geom.srcline', { A: 'X', L: fmt(S.srcLen, 3) }) + ', ' + fmt(S.E0, 3) + ' MeV', 12, 22);
    } else {
      var b0 = beam3(0); g.fillStyle = '#111';
      g.fillText(t('geom.beam') + ' ' + fmt(S.E0, 3) + ' MeV', b0[0], b0[1] - 8);
    }

    }
    function sphere3() {
      var c0 = P(0, 0, 0), rs = T * scale * D / Math.sqrt(D * D - T * T), n = 72, a, k2;
      g.fillStyle = 'rgba(30,136,229,.10)'; g.strokeStyle = '#8a97a8'; g.lineWidth = 1.5;
      g.beginPath(); g.arc(c0[0], c0[1], rs, 0, 6.2832); g.fill(); g.stroke();
      g.strokeStyle = 'rgba(100,116,139,.55)'; g.lineWidth = 1;
      for (k2 = 0; k2 < 3; k2++) {                              // three orthogonal great circles
        g.beginPath();
        for (a = 0; a <= n; a++) {
          var c = T * Math.cos(2 * Math.PI * a / n), sn = T * Math.sin(2 * Math.PI * a / n);
          var q = k2 === 0 ? P(0, c, sn) : k2 === 1 ? P(c, 0, sn) : P(c, sn, 0);
          if (a) g.lineTo(q[0], q[1]); else g.moveTo(q[0], q[1]);
        }
        g.stroke();
      }
      g.font = '12px system-ui, sans-serif'; g.fillStyle = '#4b5563'; g.textAlign = 'center';
      g.fillText(t('geom.sphere') + ': ' + matName() + '  (R = ' + fmt(T, 3) + ' cm)', W / 2, H - 8);
      g.textAlign = 'left';
    }
    function sphereSrc3() {
      var c0 = P(0, 0, 0);
      g.fillStyle = '#BA7517'; g.strokeStyle = '#fff'; g.lineWidth = 1.5;
      g.beginPath(); g.arc(c0[0], c0[1], 6, 0, 6.2832); g.fill(); g.stroke();
      g.fillStyle = '#8a5a00'; g.textAlign = 'left'; g.fillText(t('geom.src.s') + ', ' + fmt(S.E0, 3) + ' MeV', 12, 22);
    }
    function path3(p, out, dir, alpha, markers, width) {
      var pts = p.path; if (!pts || !pts.length) return;
      g.globalAlpha = alpha; g.strokeStyle = OUT_COL[out] || '#0b3d91'; g.lineWidth = width;
      g.beginPath();
      for (var i = 0; i < pts.length; i++) { var q = P(pts[i][0], pts[i][1], pts[i][3] || 0); if (i) g.lineTo(q[0], q[1]); else g.moveTo(q[0], q[1]); }
      g.stroke();
      var last = pts[pts.length - 1];
      if ((out === 'T' || out === 'R' || out === 'L') && dir) {
        var a = P(last[0], last[1], last[3] || 0), e = (sph ? 0.25 : 0.7) * T;
        var d = P(last[0] + dir[0] * e, last[1] + dir[1] * e, (last[3] || 0) + dir[2] * e);
        g.setLineDash([3, 3]); g.beginPath(); g.moveTo(a[0], a[1]); g.lineTo(d[0], d[1]); g.stroke(); g.setLineDash([]);
      }
      if (markers) for (var m = 1; m < pts.length; m++) {
        var tag = pts[m][2];
        if (typeof tag === 'number') { var mp = P(pts[m][0], pts[m][1], pts[m][3] || 0); g.fillStyle = RX_COL[tag]; g.beginPath(); g.arc(mp[0], mp[1], 4, 0, 6.2832); g.fill(); }
      }
      g.globalAlpha = 1;
    }
    var nf = S.finished.length;
    S.finished.forEach(function (p, i) {
      var a = S.mode === 'auto' ? 0.2 + 0.4 * (i / Math.max(1, nf)) : 0.25 + 0.35 * (i / Math.max(1, nf));
      path3(p, p.out, [p.u, p.v, p.w], a, S.mode === 'step', 1.2);
    });
    var cur = S.current;
    if (cur) {
      var ph = cur.ph, pts = ph.path, last = pts[pts.length - 1], drawn = { path: pts.slice() };
      if (!cur.out && (Math.abs(last[0] - ph.x) > 1e-12 || Math.abs(last[1] - ph.y) > 1e-12 || Math.abs((last[3] || 0) - ph.z) > 1e-12)) drawn.path.push([ph.x, ph.y, 'p', ph.z]);
      path3(drawn, cur.out || 'N', [ph.u, ph.v, ph.w], 1, true, 2.6);
      if (!cur.out) {
        var pp = P(ph.x, ph.y, ph.z), e = 0.16 * rad, pe = P(ph.x + ph.u * e, ph.y + ph.v * e, ph.z + ph.w * e);
        g.fillStyle = '#ffb300'; g.strokeStyle = '#111'; g.lineWidth = 1.5;
        g.beginPath(); g.arc(pp[0], pp[1], 6, 0, 6.2832); g.fill(); g.stroke();
        g.beginPath(); g.moveTo(pp[0], pp[1]); g.lineTo(pe[0], pe[1]); g.stroke();
        if (proposal && proposal.type === 'flight') {
          var pq = P(ph.x + ph.u * proposal.s, ph.y + ph.v * proposal.s, ph.z + ph.w * proposal.s);
          g.strokeStyle = '#ffb300'; g.lineWidth = 2.4; g.setLineDash([6, 4]); g.beginPath(); g.moveTo(pp[0], pp[1]); g.lineTo(pq[0], pq[1]); g.stroke(); g.setLineDash([]);
          g.fillStyle = '#ffb300'; g.beginPath(); g.arc(pq[0], pq[1], 4, 0, 6.2832); g.fill();
        }
      }
    }
    if (sph) sphereSrc3();
    // axes gizmo
    var ax = 36, ay = H - 36, axes = [['z', 1, 0, 0, '#111'], ['x', 0, 1, 0, '#2e7d32'], ['y', 0, 0, 1, '#1565c0']];
    axes.forEach(function (a) {
      var X = a[1], Y = a[2], Z = a[3], x1 = X * cyw + Z * syw, z1 = -X * syw + Z * cyw, y2 = Y * cp - z1 * sp;
      g.strokeStyle = a[4]; g.fillStyle = a[4]; g.lineWidth = 2;
      g.beginPath(); g.moveTo(ax, ay); g.lineTo(ax + x1 * 24, ay - y2 * 24); g.stroke();
      g.textAlign = 'center'; g.fillText(a[0], ax + x1 * 32, ay - y2 * 32 + 4);
    });
    g.fillStyle = '#8a97a8'; g.textAlign = 'right'; g.fillText(t('geom.hint3d'), W - 8, 16);
  }

  /* ---------------- tallies ---------------- */
  function updateTallies() {
    if (S.example === 'sphere') return updateTalliesSphere();
    var sim = S.sim, ta = sim.tally, s = MC.summary(sim);
    function row(label, f, cls) {
      return '<tr class="' + (cls || '') + '"><td>' + label + '</td><td>' + pct(f.p) + (f.n ? ' ± ' + (100 * f.err).toFixed(1) + ' %' : '') + '</td></tr>';
    }
    var h = '<table class="tally">';
    h += '<tr><td>' + t('tally.N') + '</td><td>' + ta.N + '</td></tr>';
    h += row('<b>' + t('tally.trans') + '</b>', s.transmitted);
    h += row(t('tally.unc'), s.uncollided, 'sub');
    h += '<tr class="sub"><td>' + t('tally.analytic') + '</td><td>' + pct(s.analyticUncollided) + '</td></tr>';
    h += row(t('tally.scat'), s.scattered, 'sub');
    h += row(t('tally.refl'), s.reflected);
    h += row(t('tally.abs'), s.absorbed);
    h += '<tr class="sub"><td>' + t('rx.' + RX_KEY[2] + '.s') + ' / ' + t('rx.' + RX_KEY[3] + '.s') + '</td><td>' + ta.absBy[2] + ' / ' + ta.absBy[3] + '</td></tr>';
    h += '<tr><td>' + t('tally.flux') + '<br><span style="color:#5b6573;font-size:11px">' + t('tally.flux.u') + '</span></td><td>' +
      (ta.N ? fmt(s.flux, 4) + ' ± ' + fmt(s.fluxErr, 2) : '–') + '</td></tr>';
    h += '<tr><td>' + t('tally.int') + '</td><td>' + RX_KEY.map(function (k, i) {
      return '<span class="sw" style="background:' + RX_COL[i] + '"></span>' + ta.byReaction[i];
    }).join(' &nbsp;') + '</td></tr>';
    h += '</table>';
    $('tallyTable').innerHTML = h;
    drawSpectrum(s);
  }

  function updateTalliesSphere() {
    var sim = S.sim, ta = sim.tally, s = MC.summary(sim);
    function row(label, f, cls) {
      return '<tr class="' + (cls || '') + '"><td>' + label + '</td><td>' + pct(f.p) + (f.n ? ' ± ' + (100 * f.err).toFixed(1) + ' %' : '') + '</td></tr>';
    }
    var h = '<table class="tally">';
    h += '<tr><td>' + t('tally.N') + '</td><td>' + ta.N + '</td></tr>';
    h += row('<b>' + t('sph.leak') + '</b>', s.leaked);
    h += row(t('sph.unc'), s.uncollided, 'sub');
    h += '<tr class="sub"><td>' + t('sph.analytic') + '</td><td>' + pct(s.analyticUncollided) + '</td></tr>';
    h += row(t('sph.scat'), s.scattered, 'sub');
    h += row(t('sph.abs'), s.absorbed);
    h += '<tr class="sub"><td>' + t('rx.cap.s') + ' / ' + t('rx.oth.s') + '</td><td>' + ta.absBy[2] + ' / ' + ta.absBy[3] + '</td></tr>';
    h += '<tr><td>' + t('sph.coll') + '</td><td>' + (ta.N ? fmt(s.collMean, 4) : '–') + '</td></tr>';
    h += '<tr><td>' + t('sph.flux') + '<br><span style="color:#5b6573;font-size:11px">' + t('sph.flux.u') + '</span></td><td>' +
      (ta.N ? fmt(s.flux, 4) + ' ± ' + fmt(s.fluxErr, 2) : '–') + '</td></tr>';
    h += '<tr><td>' + t('tally.int') + '</td><td>' + RX_KEY.map(function (k, i) {
      return '<span class="sw" style="background:' + RX_COL[i] + '"></span>' + ta.byReaction[i];
    }).join(' &nbsp;') + '</td></tr>';
    h += '</table>';
    $('tallyTable').innerHTML = h;
    drawFluxSpectrum();
  }

  /* bottom plot: flux spectrum in the sphere, volume average, per unit lethargy: E phi(E) with +-1 sigma bars */
  function drawFluxSpectrum() {
    if (S.example !== 'sphere' || !S.sim || !S.sim.sphere) return;
    var sp = MC.fluxSpectrum(S.sim), N = S.sim.tally.N, bars = [], errs = [], ymax = 0, ymin = Infinity, log = S.fluxLog;
    sp.forEach(function (b) {
      if (!(b.y > 0)) return;
      bars.push({ x0: b.e0, x1: b.e1, y: b.y, color: '#1e88e5' });
      if (b.y > ymax) ymax = b.y; if (b.y < ymin) ymin = b.y;
      if (N > 1 && b.err > 0) errs.push({ x: Math.sqrt(b.e0 * b.e1), y0: Math.max(b.y - b.err, b.y * 1e-3), y1: b.y + b.err });
    });
    if (ymax <= 0) { ymax = 1; ymin = 1e-3; }
    var yl = log ? [Math.pow(10, Math.floor(Math.log10(Math.max(ymin, ymax * 1e-8)))), Math.pow(10, Math.ceil(Math.log10(ymax * 1.2)))] : [0, ymax * 1.12];
    var ytop = log ? yl[1] / 1.7 : ymax * 1.12 * 0.93, vl = [{ x: 2.53e-8, color: '#6b7280', dash: [4, 4], width: 1.2 }], texts = [{ x: 2.53e-8, y: ytop, text: t('flux.th'), color: '#6b7280', dx: 6, align: 'left' }];
    if (N) {
      vl.push({ x: S.E0, color: '#c62828', width: 2, dash: [] });
      texts.push({ x: S.E0, y: ytop, text: t('flux.src', { E: fmt(S.E0, 3) }), color: '#c62828', dx: -6, align: 'right' });
    }
    fluxPlot.draw({ ml: 88, xlog: true, ylog: log, xlim: [1e-9, 30], ylim: yl, xlabel: t('flux.x'), ylabel: t('flux.y'), bars: bars, errbars: errs, vlines: vl, texts: texts });
  }

  function drawSpectrumN(s) {
    var ta = S.sim.tally, nb = ta.nbins, N = Math.max(1, ta.N), dl = (ta.lmax - ta.lmin) / nb, bars = [], ymax = 0, i;
    for (i = 0; i < nb; i++) {
      var y = ta.spec[i] / (N * dl);
      if (y > 0) bars.push({ x0: Math.exp(ta.lmin + i * dl), x1: Math.exp(ta.lmin + (i + 1) * dl), y: y, color: '#1e88e5' });
      if (y > ymax) ymax = y;
    }
    if (ymax <= 0) ymax = 1;
    var texts = [], vl = [];
    if (ta.N) {
      vl.push({ x: S.E0, color: '#c62828', width: 3, dash: [] });
      texts.push({ x: S.E0, y: ymax * 1.12 * 0.93, text: t('spec.unc', { p: pct(s.uncollided.p) }), color: '#c62828', align: S.E0 > 1e-4 ? 'right' : 'left', dx: S.E0 > 1e-4 ? -6 : 6 });
    }
    specPlot.draw({
      xlog: true, xlim: [1e-11, 30], ylim: [0, ymax * 1.12], xlabel: t('spec.x'), ylabel: t('spec.y'), bars: bars, vlines: vl, texts: texts,
      legend: [{ color: '#1e88e5', label: t('spec.scat') }], legendPos: 'left'
    });
  }
  function drawSpectrum(s) {
    if (S.particle === 'neutron') return drawSpectrumN(s);
    var ta = S.sim.tally, nb = ta.nbins, E0 = S.E0, bw = E0 / nb, N = Math.max(1, ta.N);
    var bars = [], ymax = 0, i;
    for (i = 0; i < nb; i++) { var y = ta.spec[i] / (N * bw); bars.push({ x0: i * bw, x1: (i + 1) * bw, y: y, color: '#1e88e5' }); if (y > ymax) ymax = y; }
    if (ymax <= 0) ymax = 1;
    var xmax = E0 * 1.12;
    var texts = [], vl = [];
    if (ta.N) {
      vl.push({ x: E0, color: '#c62828', width: 3, dash: [] });
      texts.push({ x: E0, y: ymax * 1.12 * 0.93, text: t('spec.unc', { p: pct(s.uncollided.p) }), color: '#c62828', align: 'right', dx: -6 });
    }
    specPlot.draw({
      xlim: [0, xmax], ylim: [0, ymax * 1.12], xlabel: t('spec.x'), ylabel: t('spec.y'), bars: bars, vlines: vl, texts: texts,
      legend: [{ color: '#1e88e5', label: t('spec.scat') }], legendPos: 'left'
    });
  }

  /* ---------------- automatic mode ---------------- */
  function autoLoop() {
    if (!S.running) return;
    var sim = S.sim, t0 = performance.now(), n = 0;
    while (sim.tally.N < S.nHist && performance.now() - t0 < 14) {
      sim.record = n < (S.example === 'sphere' ? 6 : 20); n++;
      var r = MC.runHistory(sim);
      if (sim.record) {
        S.finished.push({ path: r.ph.path, out: r.out, u: r.ph.u, v: r.ph.v, w: r.ph.w });
        if (S.finished.length > (S.example === 'sphere' ? 90 : 250)) S.finished.shift();
      }
    }
    sim.record = false;
    drawGeom(); updateTallies();
    if (sim.tally.N >= S.nHist) { stopAll(); return; }
    S.raf = requestAnimationFrame(autoLoop);
  }
  function runAuto() {
    if (S.sim.tally.N >= S.nHist) return;
    S.running = true; $('btnRun').disabled = true; $('btnPause').disabled = false;
    S.raf = requestAnimationFrame(autoLoop);
  }

  /* ---------------- step-by-step mode ---------------- */
  function archiveCurrent() {
    if (S.current && S.current.out) {
      S.finished.push({ path: S.current.ph.path, out: S.current.out, u: S.current.ph.u, v: S.current.ph.v, w: S.current.ph.w });
      if (S.finished.length > 25) S.finished.shift();
    }
    S.current = null;
  }
  function startHistory() { archiveCurrent(); S.gen = MC.historySteps(S.sim); }

  /* One user action ('n' = next step, 'p' = next particle) advances the simulation and returns the event to show.
     The actions are logged: "previous step" rebuilds the state by replaying the log minus its last action from a fresh
     simulation with the same seed (the random sequence is deterministic, so the rebuilt state is identical). */
  function advance(kind) {
    var ev, r;
    if (kind === 'p') {
      if (S.gen) {
        r = S.gen.next();
        while (!r.done) { if (r.value.type === 'end') S.current.out = r.value.out; r = S.gen.next(); }
      }
      startHistory();
      ev = S.gen.next().value;
    } else {
      r = S.gen ? S.gen.next() : { done: true };
      if (r.done) { startHistory(); r = S.gen.next(); }
      ev = r.value;
    }
    S.ev = ev;
    if (ev.type === 'birth') S.current = { ph: ev.ph, out: null };
    if (ev.type === 'end') S.current.out = ev.out;
    return ev;
  }
  function showStep(ev) { showWindow(ev); drawGeom(ev); updateTallies(); }
  function updatePrevBtn() { $('btnPrev').disabled = !S.actions || S.actions.length === 0; }
  function stepNext() { var ev = advance('n'); S.actions.push('n'); updatePrevBtn(); showStep(ev); }
  function stepNextParticle() { var ev = advance('p'); S.actions.push('p'); updatePrevBtn(); showStep(ev); }
  function stepPrev() {
    if (!S.actions.length) return;
    if (S.timer) stopAll();                                        // stop the continuous run
    var acts = S.actions.slice(0, -1);
    S.sim = newSim(); S.gen = null; S.ev = null; S.current = null; S.finished = []; S.actions = [];
    acts.forEach(function (a) { advance(a); S.actions.push(a); });
    updatePrevBtn();
    if (S.ev) showStep(S.ev); else { showReady(); drawGeom(); updateTallies(); }
  }
  function toggleContinuous() {
    if (S.timer) { stopAll(); return; }
    var d = parseInt($('delay').value, 10);
    S.timer = setInterval(stepNext, d);
    $('btnCont').classList.add('active'); $('btnCont').textContent = t('btn.pause');
  }

  /* ---------------- floating window ---------------- */
  function showReady() {
    $('win').hidden = false;
    $('winTitle').textContent = t('win.title'); $('winHist').textContent = '';
    $('winText').innerHTML = t('win.ready');
    $('winCanvas').classList.add('off'); $('winText').classList.add('big');
    S.winRender = null;
  }

  function showWindow(ev) {
    var win = $('win'); win.hidden = false;
    var ph = ev.ph, E = ph.E, hist = S.sim.history + (ev.type === 'end' ? 0 : 1);
    $('winHist').textContent = t('win.hist', { n: ev.type === 'birth' ? ev.n : hist });
    var title = '', text = '', render = null;
    var sig = ev.sig, mat = S.mat;
    switch (ev.type) {
      case 'birth':
        title = t('step.birth.title');
        text = ev.iso ? t('step.birth.text', { E: fmt(E), mu: fmt(ev.mu, 4), th: fmt(Math.acos(ev.mu) * 180 / Math.PI, 4), phi: fmt(ev.phi * 180 / Math.PI, 4) })
                      : t(ev.srcLen > 0 ? 'step.birth.text.line' : 'step.birth.text', { E: fmt(E) });
        if (ev.srcLen > 0) text += t('step.birth.line', { A: 'x', L: fmt(ev.srcLen, 3), z: fmt(ph.y, 4) });
        break;
      case 'xs_total':
        title = t('step.xs_total.title');
        text = t('step.xs_total.text', { mat: matName(), rho: mat.rho, E: fmt(E), sig: fmt(sig.tot), mfp: fmt(1 / sig.tot), tsl: mat.tsl && E < mat.tsl.EmaxMeV ? t('note.tsl') : '' });
        render = function () { plotXsTotal(ev); }; break;
      case 'flight':
        title = t('step.flight.title');
        text = t('step.flight.text', { xi: fmt(ev.xi), s: fmt(ev.s), db: fmt(ev.db) });
        render = function () { plotFlight(ev); }; break;
      case 'move':
        title = t('step.move.title');
        text = t(ev.out ? 'step.move.' + ev.out : 'step.move.in',
          { s: fmt(ev.s), db: fmt(ev.db), x: fmt(ph.x), y: fmt(ph.y), z: fmt(ph.z), r: fmt(Math.sqrt(ph.x * ph.x + ph.y * ph.y + ph.z * ph.z)), u: fmt(ph.u), E: fmt(ph.E) });
        break;
      case 'xs_partial':
        title = t('step.xs_partial.title');
        text = t('step.xs_partial.text', { E: fmt(E) }) + '<ul>' + RX_KEY.map(function (k, i) {
          var v = sig.arr[i];
          return '<li><span class="sw" style="background:' + RX_COL[i] + '"></span><span class="nm">' + t('rx.' + k) + '</span><span class="v">' +
            (v > 0 ? fmt(v) + ' cm⁻¹ (' + pct(v / sig.tot) + ')' : t('plot.below')) + '</span></li>';
        }).join('') + '<li><span class="sw" style="background:#111"></span><span class="nm"><b>Σ<sub>t</sub></b></span><span class="v"><b>' + fmt(sig.tot) + ' cm⁻¹</b></span></li></ul>';
        render = function () { plotXsPartial(ev); }; break;
      case 'nuclide_sample':
        title = t('step.nuc.title');
        text = t('step.nuc.text', { xi: fmt(ev.xi), x: fmt(ev.xi * sig.tot), nuc: mat.nuc[ev.ni].label }) + '<ul>' + mat.nuc.map(function (nu, i) {
          return '<li><span class="sw" style="background:' + NUC_COL[i % NUC_COL.length] + '"></span><span class="nm">' + nu.label + '</span><span class="v">' +
            fmt(sig.nucTot[i]) + ' cm⁻¹ (' + pct(sig.nucTot[i] / sig.tot) + ')</span></li>';
        }).join('') + '</ul>';
        render = function () { drawVector(ev); }; break;
      case 'reaction_sample':
        title = t('step.reaction.title');
        var nt = S.particle === 'neutron' ? sig.nucTot[ev.ni] : sig.tot;
        text = t('step.reaction.text', { sig: fmt(nt), xi: fmt(ev.xi), x: fmt(ev.xi * nt), rx: t('rx.' + RX_KEY[ev.r]), nuc: S.particle === 'neutron' ? mat.nuc[ev.ni].label : '' });
        render = function () { drawVector(ev); }; break;
      case 'interaction':
        title = t('step.interaction.title');
        var inf = ev.info;
        text = t(inf.thermal ? 'step.int.t' : 'step.int.' + inf.r, {
          beta: inf.thermal ? fmt(inf.detail.beta, 4) : '', alpha: inf.thermal ? fmt(inf.detail.alpha, 4) : '',
          nuc: inf.label || '', A: inf.nuc !== undefined && S.particle === 'neutron' ? fmt(S.mat.nuc[inf.nuc].A, 4) : '',
          avail: inf.detail && inf.detail.avail !== undefined ? fmt(inf.detail.avail) : '', floor: inf.detail && inf.detail.floored ? t('note.floor') : '',
          E: fmt(inf.Ein), Eout: fmt(inf.Eout || 0), mu: inf.mu !== undefined ? fmt(inf.mu) : '', tries: inf.tries,
          th: inf.mu !== undefined ? fmt(Math.acos(inf.mu) * 180 / Math.PI, 4) : '', phi: inf.phi !== undefined ? fmt(inf.phi * 180 / Math.PI, 4) : ''
        });
        if (inf.r <= 1 && !inf.thermal) { text += '<br><small>' + t('step.polar') + '</small>'; render = function () { drawPolar(ev); }; }
        break;
      case 'end':
        title = t('step.end.title');
        text = t('step.end.' + ev.out, { nc: ph.nc, E: fmt(ph.E), rx: ev.r !== undefined ? t('rx.' + RX_KEY[ev.r]) : '' });
        break;
    }
    $('winTitle').textContent = title;
    $('winText').innerHTML = text;
    var cvEl = $('winCanvas');
    cvEl.classList.toggle('off', !render);
    $('winText').classList.toggle('big', !render);
    S.winRender = render;
    if (render) requestAnimationFrame(render);
  }

  function logLim(vals) {
    var lo = Infinity, hi = 0;
    vals.forEach(function (v) { if (v > 0) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
    return [Math.pow(10, Math.floor(Math.log10(lo))), Math.pow(10, Math.ceil(Math.log10(hi)))];
  }
  function xlim() { return S.particle === 'neutron' ? [1e-11, 30] : [0.0008, 130]; }

  function plotXsTotal(ev) {
    var cur = S.mat.curves(), E = ev.ph.E, st = ev.sig.tot, yl = logLim(cur.tot);
    winPlot.draw({
      xlog: true, ylog: true, xlim: xlim(), ylim: yl, xlabel: t('plot.E'), ylabel: t('plot.sigma'),
      title: 'Σt(E): ' + matName(),
      lines: [{ x: cur.E, y: cur.tot, color: '#111', width: 2.2 }],
      vlines: [{ x: E, y1: st, color: '#c62828' }], hlines: [{ y: st, x1: E, color: '#c62828' }],
      points: [{ x: E, y: st, color: '#c62828', r: 6 }],
      texts: [{ x: E, y: st, text: 'E = ' + fmt(E, 4) + ' MeV,  Σt = ' + fmt(st, 4) + ' cm⁻¹', dx: E > (S.particle === 'neutron' ? 1e-3 : 1) ? -10 : 10, dy: -14, align: E > (S.particle === 'neutron' ? 1e-3 : 1) ? 'right' : 'left', color: '#c62828' }]
    });
  }

  function plotXsPartial(ev) {
    var cur = S.mat.curves(), E = ev.ph.E, sig = ev.sig;
    var cols = S.particle === 'neutron' ? [cur.el, cur.inel, cur.cap, cur.oth] : [cur.coh, cur.inc, cur.pe, cur.pair];
    var hi = logLim(cur.tot)[1];
    var lowVals = [];
    sig.arr.forEach(function (v) { if (v > 0) lowVals.push(v); });
    var lo = Math.pow(10, Math.floor(Math.log10(Math.min.apply(null, lowVals))) - 1);
    lo = Math.max(lo, hi * 1e-9);
    var lines = cols.map(function (c, i) { return { x: cur.E, y: c, color: RX_COL[i], width: 1.8 }; });
    lines.push({ x: cur.E, y: cur.tot, color: '#111', width: 1.6, dash: [6, 4] });
    var pts = sig.arr.map(function (v, i) { return { x: E, y: v, color: RX_COL[i], r: 6, shape: RX_SHAPE[i] }; });
    pts.push({ x: E, y: sig.tot, color: '#111', r: 5 });
    winPlot.draw({
      xlog: true, ylog: true, xlim: xlim(), ylim: [lo, hi], xlabel: t('plot.E'), ylabel: t('plot.sigma'),
      title: matName() + ': Σ(E)', lines: lines, points: pts, vlines: [{ x: E, color: '#888', dash: [4, 4], width: 1 }],
      legend: RX_KEY.map(function (k, i) { return { color: RX_COL[i], label: t('rx.' + k) }; }).concat([{ color: '#111', label: t('plot.tot') }]),
      legendPos: 'left'
    });
  }

  function plotFlight(ev) {
    var sig = ev.sig.tot, s = ev.s, db = ev.db, xi = ev.xi;
    var xmax = Math.max(4 / sig, s * 1.25, Math.min(db * 1.12, 12 / sig));
    var xs = [], ys = [], n = 240, i;
    for (i = 0; i <= n; i++) { var x = xmax * i / n; xs.push(x); ys.push(1 - Math.exp(-sig * x)); }
    var vl = [{ x: s, y1: xi, color: '#c62828' }];
    if (db <= xmax) vl.push({ x: db, color: '#ef6c00', dash: [2, 3], width: 2 });
    winPlot.draw({
      xlim: [0, xmax], ylim: [0, 1.02], xlabel: t('plot.s'), ylabel: t('plot.F'), title: 'F(s)',
      lines: [{ x: xs, y: ys, color: '#1e5fd0', width: 2.4 }],
      hlines: [{ y: xi, x1: s, color: '#c62828' }], vlines: vl,
      points: [{ x: s, y: xi, color: '#c62828', r: 6 }],
      texts: [
        { x: 0, y: xi, text: 'ξ = ' + fmt(xi, 4), dx: 6, dy: -9, color: '#c62828' },
        { x: s, y: 0.04, text: 's = ' + fmt(s, 4) + ' cm', dx: s > xmax * 0.6 ? -6 : 6, align: s > xmax * 0.6 ? 'right' : 'left', color: '#c62828' }
      ].concat(db <= xmax ? [{ x: db, y: 0.5, text: t('plot.boundary') + ' (' + fmt(db, 3) + ' cm)', dx: db > xmax * 0.6 ? -6 : 6, align: db > xmax * 0.6 ? 'right' : 'left', color: '#ef6c00' }] : [])
    });
  }

  function vectorCfg(ev) {
    var sig = ev.sig, c = {}, i;
    if (ev.type === 'nuclide_sample') {
      c.vals = sig.nucTot; c.tot = sig.tot; c.sel = ev.ni; c.labels = S.mat.nuc.map(function (n) { return n.label; });
      c.cols = c.vals.map(function (v, k) { return NUC_COL[k % NUC_COL.length]; });
      c.head = 'Σᵢ per nuclide  →  Σt = ' + fmt(sig.tot) + ' cm⁻¹';
    } else if (S.particle === 'neutron') {
      c.vals = sig.nucArr[ev.ni]; c.tot = sig.nucTot[ev.ni]; c.sel = ev.r; c.cols = RX_COL;
      c.labels = RX_KEY.map(function (k) { return t('rx.' + k + '.s'); });
      c.head = S.mat.nuc[ev.ni].label + ':  [Σel, Σinel, Σcap, Σoth]  →  Σ = ' + fmt(c.tot) + ' cm⁻¹';
    } else {
      c.vals = sig.arr; c.tot = sig.tot; c.sel = ev.r; c.cols = RX_COL;
      c.labels = RX_KEY.map(function (k) { return t('rx.' + k + '.s'); });
      c.head = '[Σcoh, Σinc, Σpe, Σpair]  →  Σt = ' + fmt(sig.tot) + ' cm⁻¹';
    }
    return c;
  }
  function drawVector(ev) {
    var cv = $('winCanvas'), c = sizeCanvas(cv), g = c.g, W = c.W, H = c.H, cfg = vectorCfg(ev), sig = { tot: cfg.tot };
    var ml = 24, mr = 24, bw = W - ml - mr, by = H * 0.40, bh = Math.min(64, H * 0.22), x = ml, i, n = cfg.vals.length, lab = 0;
    g.font = '12px system-ui, sans-serif'; g.textBaseline = 'middle';
    g.fillStyle = '#111'; g.textAlign = 'center'; g.font = 'bold 13px system-ui, sans-serif';
    g.fillText(cfg.head, W / 2, 18);
    g.font = '12px system-ui, sans-serif';
    for (i = 0; i < n; i++) {
      var f = cfg.vals[i] / cfg.tot, w = bw * f;
      g.fillStyle = cfg.cols[i]; g.globalAlpha = i === cfg.sel ? 1 : 0.55; g.fillRect(x, by, Math.max(w, 0), bh); g.globalAlpha = 1;
      g.strokeStyle = '#fff'; g.lineWidth = 1; g.strokeRect(x, by, Math.max(w, 0), bh);
      if (n <= 4 || f >= 0.03 || i === cfg.sel) {
        var yy = by + bh + 22 + (lab++ % 2) * 32;                  // staggered labels
        g.fillStyle = cfg.cols[i]; g.textAlign = 'center';
        var cx = Math.min(W - 70, Math.max(70, x + w / 2));
        g.fillText(cfg.labels[i], cx, yy); g.fillText(fmt(cfg.vals[i], 3) + '  (' + pct(f) + ')', cx, yy + 13);
        if (w > 14) { g.strokeStyle = cfg.cols[i]; g.globalAlpha = 0.5; g.beginPath(); g.moveTo(x + w / 2, by + bh); g.lineTo(cx, yy - 8); g.stroke(); g.globalAlpha = 1; }
      }
      x += w;
    }
    // pointer at xi * Sigma_t
    var px = ml + bw * ev.xi;
    g.fillStyle = '#111'; g.beginPath(); g.moveTo(px, by - 2); g.lineTo(px - 7, by - 16); g.lineTo(px + 7, by - 16); g.closePath(); g.fill();
    g.strokeStyle = '#111'; g.lineWidth = 2; g.beginPath(); g.moveTo(px, by - 2); g.lineTo(px, by + bh + 4); g.stroke();
    g.fillStyle = '#111'; g.textAlign = px > W / 2 ? 'right' : 'left';
    g.fillText('ξ = ' + fmt(ev.xi, 4) + '  →  ξ·Σt = ' + fmt(ev.xi * sig.tot, 4) + ' cm⁻¹', px > W / 2 ? px + 6 : px - 6, by - 26);
  }

  var elCache = {};
  function elasticPDF(A) {                                          // lab-frame cosine density of isotropic-CM elastic scattering
    var key = A.toFixed(3);
    if (!elCache[key]) {
      var B = 120, K = 60000, f = new Float64Array(B), k;
      for (k = 0; k < K; k++) {
        var cc = -1 + 2 * (k + 0.5) / K, mu = (1 + A * cc) / Math.sqrt(A * A + 1 + 2 * A * cc);
        f[Math.min(B - 1, Math.max(0, Math.floor((mu + 1) / 2 * B)))] += 1;
      }
      for (k = 0; k < B; k++) f[k] /= K * 2 / B;
      elCache[key] = function (m) { return f[Math.min(B - 1, Math.max(0, Math.floor((m + 1) / 2 * B)))]; };
    }
    return elCache[key];
  }

  function drawPolar(ev) {
    var cv = $('winCanvas'), c = sizeCanvas(cv), g = c.g, W = c.W, H = c.H, inf = ev.info;
    var cx = W * 0.42, cy = H / 2, R = Math.min(W * 0.36, H * 0.42), n = 360, i;
    var pdf;
    if (S.particle === 'neutron') pdf = inf.r === 0 ? elasticPDF(S.mat.nuc[inf.nuc].A) : function () { return 0.5; };
    else pdf = inf.r === 1 ? function (mu) { return MC.kleinNishinaPDF(inf.Ein, mu); } : function (mu) { return MC.rayleighPDF(inf.Ein, S.mat.zeff, mu); };
    var fmax = 1e-300, f = [];
    for (i = 0; i <= n; i++) { var th = Math.PI * i / n; f.push(pdf(Math.cos(th))); if (f[i] > fmax) fmax = f[i]; }
    // rings and axis
    g.strokeStyle = '#e3e6ea'; g.lineWidth = 1;
    [0.5, 1].forEach(function (k) { g.beginPath(); g.arc(cx, cy, R * k, 0, 6.2832); g.stroke(); });
    g.beginPath(); g.moveTo(cx - R * 1.05, cy); g.lineTo(cx + R * 1.1, cy); g.stroke();
    // distribution (symmetrical about the incident axis)
    g.fillStyle = RX_COL[inf.r]; g.strokeStyle = RX_COL[inf.r]; g.lineWidth = 2;
    g.beginPath();
    for (i = 0; i <= n; i++) { var a = Math.PI * i / n, r = R * f[i] / fmax; g.lineTo(cx + r * Math.cos(a), cy - r * Math.sin(a)); }
    for (i = n; i >= 0; i--) { var a2 = Math.PI * i / n, r2 = R * f[i] / fmax; g.lineTo(cx + r2 * Math.cos(a2), cy + r2 * Math.sin(a2)); }
    g.closePath(); g.globalAlpha = 0.25; g.fill(); g.globalAlpha = 1; g.stroke();
    // incident direction
    g.strokeStyle = '#111'; g.fillStyle = '#111'; g.lineWidth = 2;
    g.beginPath(); g.moveTo(cx - R * 1.02, cy); g.lineTo(cx - 4, cy); g.stroke();
    g.beginPath(); g.moveTo(cx, cy); g.lineTo(cx - 11, cy - 5); g.lineTo(cx - 11, cy + 5); g.closePath(); g.fill();
    // sampled direction
    var th0 = Math.acos(inf.mu) * (Math.sin(inf.phi) >= 0 ? 1 : -1);
    var ex = cx + R * Math.cos(th0), ey = cy - R * Math.sin(th0);
    g.strokeStyle = '#c62828'; g.fillStyle = '#c62828'; g.lineWidth = 3;
    g.beginPath(); g.moveTo(cx, cy); g.lineTo(ex, ey); g.stroke();
    g.beginPath(); g.arc(ex, ey, 5, 0, 6.2832); g.fill();
    g.font = '12px system-ui, sans-serif'; g.textBaseline = 'middle'; g.fillStyle = '#111'; g.textAlign = 'left';
    g.fillText('θ = ' + fmt(Math.acos(inf.mu) * 180 / Math.PI, 4) + '°', W * 0.78, cy - 22);
    g.fillText('E = ' + fmt(inf.Ein, 4) + ' MeV', W * 0.78, cy);
    g.fillText("E' = " + fmt(inf.Eout, 4) + ' MeV', W * 0.78, cy + 18);
    g.textAlign = 'center'; g.fillStyle = '#5b6573';
    g.fillText(t('rx.' + RX_KEY[inf.r]), cx, 14);
  }

  /* window drag */
  (function () {
    var bar = $('winbar'), win = $('win'), dx = 0, dy = 0, drag = false;
    bar.addEventListener('pointerdown', function (e) {
      if (e.target.tagName === 'BUTTON') return;
      var r = win.getBoundingClientRect(); dx = e.clientX - r.left; dy = e.clientY - r.top; drag = true; bar.setPointerCapture(e.pointerId);
    });
    bar.addEventListener('pointermove', function (e) {
      if (!drag) return;
      win.style.left = Math.max(0, Math.min(window.innerWidth - 80, e.clientX - dx)) + 'px';
      win.style.top = Math.max(0, Math.min(window.innerHeight - 40, e.clientY - dy)) + 'px'; win.style.right = 'auto';
    });
    bar.addEventListener('pointerup', function () { drag = false; });
    $('winClose').addEventListener('click', function () { win.hidden = true; });
    if (window.ResizeObserver) new ResizeObserver(function () { if (S.winRender && !win.hidden) S.winRender(); }).observe($('winCanvas'));
  })();

  /* ---------------- UI wiring ---------------- */
  function applyLang() {
    document.documentElement.lang = lang(); document.title = t('app.title');
    document.querySelectorAll('[data-i18n]').forEach(function (el) { el.innerHTML = t(el.getAttribute('data-i18n')); });
    document.querySelectorAll('.lang button').forEach(function (b) { b.classList.toggle('on', b.dataset.lang === lang()); });
    var sel = $('material'), cur = sel.value;
    sel.innerHTML = MAT_ORDER.map(function (k) { return '<option value="' + k + '">' + (window.XCOM[k].name[lang()]) + ' (ρ = ' + window.XCOM[k].rho + ' g/cm³)</option>'; }).join('');
    sel.value = cur || 'lead';
    $('geomLegend').dataset.lang = '';
    $('btnCont').textContent = S.timer ? t('btn.pause') : t('btn.cont');
    $('hint').innerHTML = S.mode === 'step' ? t('hint.step') : t('hint.auto');
    if (S.sim) { updateInfo(); drawGeom(); updateTallies(); if (!$('win').hidden) { if (S.ev) showWindow(S.ev); else showReady(); } }
  }
  function setMode(m) {
    S.mode = m;
    $('autoBtns').hidden = m !== 'auto'; $('stepBtns').hidden = m !== 'step';
    $('hint').innerHTML = m === 'step' ? t('hint.step') : t('hint.auto');
    resetSim();
  }

  document.querySelectorAll('.lang button').forEach(function (b) {
    b.addEventListener('click', function () { window.I18N.setLang(b.dataset.lang); applyLang(); });
  });
  ['material', 'energy', 'thickness', 'nhist', 'seed', 'srctype', 'srclen', 'tslChk'].forEach(function (id) { $(id).addEventListener('change', resetSim); });
  document.querySelectorAll('input[name=mode]').forEach(function (r) { r.addEventListener('change', function () { if (r.checked) setMode(r.value); }); });
  $('btnRun').addEventListener('click', function () { readParamsKeep(); runAuto(); });
  function readParamsKeep() { S.nHist = Math.max(1, Math.floor(parseFloat($('nhist').value) || 1000)); }
  $('btnPause').addEventListener('click', stopAll);
  $('btnShowWin').addEventListener('click', function () { if (S.ev) showWindow(S.ev); else showReady(); });
  document.querySelectorAll('.viewtoggle button').forEach(function (b) {
    b.addEventListener('click', function () {
      S.view = b.dataset.view;
      document.querySelectorAll('.viewtoggle button').forEach(function (x) { x.classList.toggle('on', x === b); });
      $('geom').classList.toggle('rot', S.view === '3d');
      redraw();
    });
  });
  (function () {                                   // 2D view: wheel = zoom about the cursor, drag = pan, double-click = reset
    var cv = $('geom'), drag = false, lx = 0, ly = 0;
    cv.addEventListener('pointerdown', function (e) { if (S.view === '3d') return; drag = true; lx = e.clientX; ly = e.clientY; cv.setPointerCapture(e.pointerId); cv.style.cursor = 'grabbing'; });
    cv.addEventListener('pointermove', function (e) {
      if (!drag || S.view === '3d') return;
      cam2.px += e.clientX - lx; cam2.py += e.clientY - ly; lx = e.clientX; ly = e.clientY; redraw();
    });
    cv.addEventListener('pointerup', function () { drag = false; cv.style.cursor = ''; });
    cv.addEventListener('dblclick', function () { if (S.view === '3d') return; cam2.zoom = 1; cam2.px = 0; cam2.py = 0; redraw(); });
    cv.addEventListener('wheel', function (e) {
      if (S.view === '3d') return; e.preventDefault();
      var r = cv.getBoundingClientRect(), cx = e.clientX - r.left, cy = e.clientY - r.top, W = r.width, H = r.height;
      var vp = view2(W, H), base = vp.base, sOld = base * cam2.zoom, x0 = vp.x0;
      var wx = (cx - x0 - cam2.px) / sOld, wy = (H / 2 + cam2.py - cy) / sOld;   // world point under the cursor
      cam2.zoom = Math.max(0.2, Math.min(60, cam2.zoom * Math.exp(-e.deltaY * 0.001)));
      var sNew = base * cam2.zoom;
      cam2.px = cx - x0 - wx * sNew; cam2.py = cy - H / 2 + wy * sNew; redraw();
    }, { passive: false });
  })();
  (function () {                                   // orbit camera
    var cv = $('geom'), drag = false, lx = 0, ly = 0;
    cv.addEventListener('pointerdown', function (e) { if (S.view !== '3d') return; drag = true; lx = e.clientX; ly = e.clientY; cv.setPointerCapture(e.pointerId); });
    cv.addEventListener('pointermove', function (e) {
      if (!drag) return;
      cam.yaw += (e.clientX - lx) * 0.01; cam.pitch = Math.max(-1.5, Math.min(1.5, cam.pitch + (e.clientY - ly) * 0.01));
      lx = e.clientX; ly = e.clientY; redraw();
    });
    cv.addEventListener('pointerup', function () { drag = false; });
    cv.addEventListener('wheel', function (e) {
      if (S.view !== '3d') return; e.preventDefault();
      cam.zoom = Math.max(0.3, Math.min(5, cam.zoom * Math.exp(-e.deltaY * 0.001))); redraw();
    }, { passive: false });
  })();
  document.querySelectorAll('#particleToggle button').forEach(function (b) {
    b.addEventListener('click', function () { setParticle(b.dataset.particle); });
  });
  document.querySelectorAll('#exToggle button').forEach(function (b) { b.addEventListener('click', function () { setExample(b.dataset.ex); }); });
  document.querySelectorAll('#fluxScale button').forEach(function (b) {
    b.addEventListener('click', function () {
      S.fluxLog = b.dataset.yscale === 'log';
      document.querySelectorAll('#fluxScale button').forEach(function (x) { x.classList.toggle('on', x === b); });
      drawFluxSpectrum();
    });
  });
  $('btnReset').addEventListener('click', resetSim);
  $('btnPrev').addEventListener('click', stepPrev);
  $('btnNext').addEventListener('click', stepNext);
  $('btnNextP').addEventListener('click', stepNextParticle);
  $('btnCont').addEventListener('click', toggleContinuous);
  $('delay').addEventListener('input', function () {
    $('delayOut').textContent = $('delay').value + ' ms';
    if (S.timer) { clearInterval(S.timer); S.timer = setInterval(stepNext, parseInt($('delay').value, 10)); }
  });
  document.addEventListener('keydown', function (e) {
    if (S.mode !== 'step') return;
    var tag = (e.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'select' || tag === 'textarea') return;
    if (e.code === 'Space' || e.key === 'ArrowRight') { e.preventDefault(); stepNext(); }
    else if (e.key === 'b' || e.key === 'B' || e.key === 'ArrowLeft') { e.preventDefault(); stepPrev(); }
    else if (e.key === 'p' || e.key === 'P') stepNextParticle();
  });
  window.addEventListener('resize', function () { drawGeom(); updateTallies(); if (S.winRender && !$('win').hidden) S.winRender(); });

  window.I18N.setLang((navigator.language || 'en').slice(0, 2) === 'es' ? 'es' : 'en');
  applyLang();
  $('material').value = 'lead';
  resetSim();
  if (location.hash === '#sphere') setExample('sphere');          // direct link to example 2
})();
