/* Minimal canvas plotting (linear / log axes, lines, markers, guide lines, legend). No dependencies. */
(function (root) {
  'use strict';
  var SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
  function pow10(n) { return '10' + String(n).split('').map(function (c) { return SUP[c] || c; }).join(''); }
  function fmt(v) {
    if (v === 0) return '0';
    var a = Math.abs(v);
    if (a >= 1e4 || a < 1e-3) return v.toExponential(2).replace('e+', 'e');
    return String(Number(v.toPrecision(4)));
  }

  function niceTicks(lo, hi, n) {
    var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var err = span / n / step;
    if (err >= 5) step *= 5; else if (err >= 2) step *= 2;
    var t = [], v = Math.ceil(lo / step - 1e-9) * step;
    for (; v <= hi + step * 1e-9; v += step) t.push(Number(v.toPrecision(12)));
    return t;
  }

  function Plot(canvas) { this.c = canvas; }

  Plot.prototype.draw = function (o) {
    var cv = this.c, dpr = window.devicePixelRatio || 1;
    var W = cv.clientWidth, H = cv.clientHeight;
    if (!W || !H) return;
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    var g = cv.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    var m = { l: 62, r: 14, t: o.title ? 26 : 12, b: 44 };
    var pw = W - m.l - m.r, ph = H - m.t - m.b;
    var xl = o.xlim, yl = o.ylim;
    var tx = o.xlog ? Math.log10 : function (v) { return v; };
    var ty = o.ylog ? Math.log10 : function (v) { return v; };
    var x0 = tx(xl[0]), x1 = tx(xl[1]), y0 = ty(yl[0]), y1 = ty(yl[1]);
    function X(v) { return m.l + (tx(v) - x0) / (x1 - x0) * pw; }
    function Y(v) { return m.t + ph - (ty(v) - y0) / (y1 - y0) * ph; }
    this.X = X; this.Y = Y;
    g.font = '12px system-ui, sans-serif';
    g.lineJoin = 'round';

    // grid + ticks
    function ticks(lo, hi, log) {
      var r = [];
      if (log) {
        var st = Math.max(1, Math.ceil((hi - lo) / 7)); for (var e = Math.ceil(lo - 1e-9); e <= hi + 1e-9; e += st) r.push({ v: Math.pow(10, e), label: pow10(e), major: true });
        if (hi - lo < 4) for (var e2 = Math.floor(lo); e2 < hi; e2++) for (var k = 2; k < 10; k++) { var v = k * Math.pow(10, e2); if (Math.log10(v) > lo && Math.log10(v) < hi) r.push({ v: v, major: false }); }
      } else niceTicks(lo, hi, 6).forEach(function (v) { r.push({ v: v, label: fmt(v), major: true }); });
      return r;
    }
    g.strokeStyle = '#e3e6ea'; g.fillStyle = '#4b5563'; g.lineWidth = 1;
    g.textAlign = 'center'; g.textBaseline = 'top';
    ticks(x0, x1, o.xlog).forEach(function (t) {
      var px = X(t.v); if (px < m.l - 1 || px > m.l + pw + 1) return;
      g.globalAlpha = t.major ? 1 : 0.4; g.beginPath(); g.moveTo(px, m.t); g.lineTo(px, m.t + ph); g.stroke(); g.globalAlpha = 1;
      if (t.label) g.fillText(t.label, px, m.t + ph + 5);
    });
    g.textAlign = 'right'; g.textBaseline = 'middle';
    ticks(y0, y1, o.ylog).forEach(function (t) {
      var py = Y(t.v); if (py < m.t - 1 || py > m.t + ph + 1) return;
      g.globalAlpha = t.major ? 1 : 0.4; g.beginPath(); g.moveTo(m.l, py); g.lineTo(m.l + pw, py); g.stroke(); g.globalAlpha = 1;
      if (t.label) g.fillText(t.label, m.l - 6, py);
    });
    g.strokeStyle = '#9aa3ad'; g.strokeRect(m.l, m.t, pw, ph);

    // clip to plot area
    g.save(); g.beginPath(); g.rect(m.l, m.t, pw, ph); g.clip();

    (o.bars || []).forEach(function (b) {          // {x0,x1,y,color}
      var a = X(b.x0), c = X(b.x1), top = Y(b.y), bot = Y(yl[0]);
      g.fillStyle = b.color; g.fillRect(a, top, Math.max(1, c - a - 0.5), bot - top);
    });
    (o.lines || []).forEach(function (L) {
      g.strokeStyle = L.color; g.lineWidth = L.width || 1.6; g.setLineDash(L.dash || []);
      g.beginPath();
      var pen = false;
      for (var i = 0; i < L.x.length; i++) {
        var yv = L.y[i];
        if (!(yv > 0) && o.ylog) { pen = false; continue; }
        var px = X(L.x[i]), py = Y(yv);
        if (!pen) { g.moveTo(px, py); pen = true; } else g.lineTo(px, py);
      }
      g.stroke(); g.setLineDash([]);
    });
    (o.vlines || []).concat(o.hlines || []).forEach(function (L) {
      g.strokeStyle = L.color || '#111'; g.lineWidth = L.width || 1.2; g.setLineDash(L.dash || [5, 4]);
      g.beginPath();
      if (L.x !== undefined) { g.moveTo(X(L.x), m.t + (L.y1 !== undefined ? Y(L.y1) - m.t : 0)); g.lineTo(X(L.x), L.y0 !== undefined ? Y(L.y0) : m.t + ph); }
      else { g.moveTo(L.x0 !== undefined ? X(L.x0) : m.l, Y(L.y)); g.lineTo(L.x1 !== undefined ? X(L.x1) : m.l + pw, Y(L.y)); }
      g.stroke(); g.setLineDash([]);
    });
    g.restore();

    // markers (outside clip so they are never cut)
    (o.points || []).forEach(function (P) {
      if (P.y <= 0 && o.ylog) return;
      var px = X(P.x), py = Y(P.y), r = P.r || 5;
      g.fillStyle = P.color; g.strokeStyle = '#fff'; g.lineWidth = 1.5;
      g.beginPath();
      if (P.shape === 'square') g.rect(px - r, py - r, 2 * r, 2 * r);
      else if (P.shape === 'diamond') { g.moveTo(px, py - r * 1.3); g.lineTo(px + r * 1.3, py); g.lineTo(px, py + r * 1.3); g.lineTo(px - r * 1.3, py); g.closePath(); }
      else if (P.shape === 'triangle') { g.moveTo(px, py - r * 1.2); g.lineTo(px + r * 1.2, py + r); g.lineTo(px - r * 1.2, py + r); g.closePath(); }
      else g.arc(px, py, r, 0, 2 * Math.PI);
      g.fill(); g.stroke();
    });

    // axis labels / title
    g.fillStyle = '#111'; g.textAlign = 'center'; g.textBaseline = 'alphabetic';
    if (o.xlabel) g.fillText(o.xlabel, m.l + pw / 2, H - 6);
    if (o.ylabel) { g.save(); g.translate(14, m.t + ph / 2); g.rotate(-Math.PI / 2); g.fillText(o.ylabel, 0, 0); g.restore(); }
    if (o.title) { g.font = 'bold 13px system-ui, sans-serif'; g.fillText(o.title, m.l + pw / 2, 17); g.font = '12px system-ui, sans-serif'; }

    // text annotations
    (o.texts || []).forEach(function (T) {
      var px = X(T.x) + (T.dx || 0), py = Y(T.y) + (T.dy || 0);
      g.textAlign = T.align || 'left'; g.textBaseline = 'middle';
      g.lineWidth = 3; g.strokeStyle = 'rgba(255,255,255,.9)'; g.strokeText(T.text, px, py);
      g.fillStyle = T.color || '#111'; g.fillText(T.text, px, py);
    });

    // legend
    if (o.legend) {
      var items = o.legend, lw = 0, lh = items.length * 16 + 8;
      items.forEach(function (it) { lw = Math.max(lw, g.measureText(it.label).width); });
      lw += 36;
      var lx = o.legendPos === 'left' ? m.l + 8 : m.l + pw - lw - 6, ly = m.t + 6;
      g.fillStyle = 'rgba(255,255,255,.88)'; g.strokeStyle = '#c9ced4'; g.lineWidth = 1;
      g.fillRect(lx, ly, lw, lh); g.strokeRect(lx, ly, lw, lh);
      g.textAlign = 'left'; g.textBaseline = 'middle';
      items.forEach(function (it, i) {
        var yy = ly + 12 + i * 16;
        g.strokeStyle = it.color; g.lineWidth = 2.5; g.beginPath(); g.moveTo(lx + 6, yy); g.lineTo(lx + 24, yy); g.stroke();
        g.fillStyle = '#111'; g.fillText(it.label, lx + 29, yy);
      });
    }
  };

  root.Plot = Plot;
})(window);
