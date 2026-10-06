/* Sky Sketch — shape matcher
 *
 * Compares a hand-drawn star path against constellation stick figures.
 *
 * Both shapes are lifted out of screen space first: every figure is
 * gnomonically projected onto the tangent plane at its own centroid, which is
 * the shape you'd see looking straight at it. That makes the comparison
 * independent of where the drawing happened to sit on screen.
 *
 * The comparison itself is full Procrustes on arc-length-resampled paths,
 * written in complex form: with both point sets centred and scaled to unit
 * Frobenius norm, |sum z_k * conj(w_k)| is the best achievable correlation over
 * all rotations and scales. 1.0 means the two outlines are the same shape.
 */
(function (root) {
  'use strict';

  var N = 64;           // Procrustes resample resolution
  var TURN_N = 96;      // turning-function resolution
  var TURN_SCALE = 0.85; // radians of RMS turning error that halves the score
  var PRO_POWER = 2.2;  // how hard the global-form term bites
  var REFLECT_PENALTY = 0.96;

  // ---------- geometry ----------

  function tangentFrame(vecs) {
    var cx = 0, cy = 0, cz = 0, i;
    for (i = 0; i < vecs.length; i++) { cx += vecs[i][0]; cy += vecs[i][1]; cz += vecs[i][2]; }
    var n = Math.sqrt(cx * cx + cy * cy + cz * cz) || 1;
    var c = [cx / n, cy / n, cz / n];
    var e = [-c[1], c[0], 0];
    var en = Math.sqrt(e[0] * e[0] + e[1] * e[1]);
    if (en < 1e-8) { e = [1, 0, 0]; en = 1; }
    e = [e[0] / en, e[1] / en, 0];
    var up = [c[1] * e[2] - c[2] * e[1], c[2] * e[0] - c[0] * e[2], c[0] * e[1] - c[1] * e[0]];
    return { c: c, e: e, n: up };
  }

  function planeProject(v, tp) {
    var z = v[0] * tp.c[0] + v[1] * tp.c[1] + v[2] * tp.c[2];
    if (z < 0.08) return null;
    return [
      (v[0] * tp.e[0] + v[1] * tp.e[1] + v[2] * tp.e[2]) / z,
      (v[0] * tp.n[0] + v[1] * tp.n[1] + v[2] * tp.n[2]) / z
    ];
  }

  /** Angular radius (deg) of a set of unit vectors about their centroid. */
  function spread(vecs) {
    var tp = tangentFrame(vecs), worst = 0, i, d;
    for (i = 0; i < vecs.length; i++) {
      d = vecs[i][0] * tp.c[0] + vecs[i][1] * tp.c[1] + vecs[i][2] * tp.c[2];
      d = Math.acos(Math.max(-1, Math.min(1, d)));
      if (d > worst) worst = d;
    }
    return { centre: tp.c, radiusDeg: worst * 180 / Math.PI };
  }

  // ---------- path shaping ----------

  function resample(pts, count) {
    if (!pts || pts.length < 2) return null;
    var seg = [], total = 0, i, d;
    for (i = 1; i < pts.length; i++) {
      d = Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]);
      seg.push(d); total += d;
    }
    if (total < 1e-12) return null;
    var out = [pts[0].slice()];
    var step = total / (count - 1);
    var idx = 1, walked = 0;
    for (var k = 1; k < count - 1; k++) {
      var target = k * step;
      while (idx < pts.length - 1 && walked + seg[idx - 1] < target) { walked += seg[idx - 1]; idx++; }
      var len = seg[idx - 1] || 1e-12;
      var t = Math.max(0, Math.min(1, (target - walked) / len));
      out.push([
        pts[idx - 1][0] + (pts[idx][0] - pts[idx - 1][0]) * t,
        pts[idx - 1][1] + (pts[idx][1] - pts[idx - 1][1]) * t
      ]);
    }
    out.push(pts[pts.length - 1].slice());
    return out;
  }

  /** Centre at origin, scale to unit Frobenius norm, flatten to [x0,y0,...]. */
  function normalise(pts) {
    if (!pts) return null;
    var n = pts.length, mx = 0, my = 0, i;
    for (i = 0; i < n; i++) { mx += pts[i][0]; my += pts[i][1]; }
    mx /= n; my /= n;
    var out = new Float64Array(n * 2), ss = 0;
    for (i = 0; i < n; i++) {
      var x = pts[i][0] - mx, y = pts[i][1] - my;
      out[i * 2] = x; out[i * 2 + 1] = y;
      ss += x * x + y * y;
    }
    var s = Math.sqrt(ss);
    if (s < 1e-12) return null;
    for (i = 0; i < out.length; i++) out[i] /= s;
    return out;
  }

  function reflectY(z) {
    var o = new Float64Array(z.length);
    for (var i = 0; i < z.length; i += 2) { o[i] = z[i]; o[i + 1] = -z[i + 1]; }
    return o;
  }

  function reverseOrder(z) {
    var n = z.length / 2, o = new Float64Array(z.length);
    for (var i = 0; i < n; i++) { o[i * 2] = z[(n - 1 - i) * 2]; o[i * 2 + 1] = z[(n - 1 - i) * 2 + 1]; }
    return o;
  }

  /** |<z,w>| over unit-norm shapes = best similarity-transform correlation. */
  function correlate(z, w) {
    var re = 0, im = 0;
    for (var i = 0; i < z.length; i += 2) {
      re += z[i] * w[i] + z[i + 1] * w[i + 1];
      im += z[i + 1] * w[i] - z[i] * w[i + 1];
    }
    return Math.hypot(re, im);
  }

  /**
   * Turning function: cumulative heading change sampled at equal arc length.
   * Invariant to translation and scale; rotation shows up as a constant offset,
   * which turnDistance divides out. Far more discriminating than Procrustes
   * alone, because it reads the corners rather than the overall drift.
   */
  function turningProfile(pts, count) {
    if (!pts || pts.length < 2) return null;
    var segLen = [], segAng = [], i;
    for (i = 1; i < pts.length; i++) {
      var dx = pts[i][0] - pts[i - 1][0], dy = pts[i][1] - pts[i - 1][1];
      var L = Math.hypot(dx, dy);
      if (L < 1e-12) continue;
      segLen.push(L); segAng.push(Math.atan2(dy, dx));
    }
    if (!segLen.length) return null;
    var cum = [segAng[0]];
    for (i = 1; i < segAng.length; i++) {
      var d = segAng[i] - segAng[i - 1];
      while (d > Math.PI) d -= 2 * Math.PI;
      while (d < -Math.PI) d += 2 * Math.PI;
      cum.push(cum[i - 1] + d);
    }
    var total = 0;
    for (i = 0; i < segLen.length; i++) total += segLen[i];
    var bounds = [], run = 0;
    for (i = 0; i < segLen.length; i++) { run += segLen[i]; bounds.push(run / total); }
    var out = new Float64Array(count), k = 0;
    for (var j = 0; j < count; j++) {
      var s = (j + 0.5) / count;
      while (k < bounds.length - 1 && s > bounds[k]) k++;
      out[j] = cum[k];
    }
    return out;
  }

  function negate(p) {
    var o = new Float64Array(p.length);
    for (var i = 0; i < p.length; i++) o[i] = -p[i];
    return o;
  }

  /** RMS turning-function distance, minimised over rotation (a constant offset). */
  function turnDistance(a, b) {
    var i, d = 0;
    for (i = 0; i < a.length; i++) d += a[i] - b[i];
    var m = d / a.length, s = 0;
    for (i = 0; i < a.length; i++) { var e = a[i] - b[i] - m; s += e * e; }
    return Math.sqrt(s / a.length);
  }

  /** Complex rotate+scale taking z onto w, for drawing the overlay. */
  function alignment(z, w) {
    var re = 0, im = 0;
    for (var i = 0; i < z.length; i += 2) {
      re += w[i] * z[i] + w[i + 1] * z[i + 1];
      im += w[i + 1] * z[i] - w[i] * z[i + 1];
    }
    return { re: re, im: im };
  }

  // ---------- turning polyline sets into candidate strokes ----------

  function edgeKey(a, b) { return a < b ? a + ':' + b : b + ':' + a; }

  function adjacency(polys) {
    var adj = new Map(), i, j;
    for (i = 0; i < polys.length; i++) {
      for (j = 1; j < polys[i].length; j++) {
        var a = polys[i][j - 1], b = polys[i][j];
        if (a === b) continue;
        if (!adj.has(a)) adj.set(a, []);
        if (!adj.has(b)) adj.set(b, []);
        if (adj.get(a).indexOf(b) < 0) adj.get(a).push(b);
        if (adj.get(b).indexOf(a) < 0) adj.get(b).push(a);
      }
    }
    return adj;
  }

  /** One pen-stroke through the figure's graph, starting from `start`. */
  function greedyWalk(adj, start) {
    var walk = [start], used = new Set(), node = start, guard = 0, i;
    while (guard++ < 400) {
      var nbrs = adj.get(node) || [], next = null, nd = -1;
      for (i = 0; i < nbrs.length; i++) {
        if (used.has(edgeKey(node, nbrs[i]))) continue;
        var deg = (adj.get(nbrs[i]) || []).length;
        if (next === null || deg > nd) { next = nbrs[i]; nd = deg; }
      }
      if (next === null) break;
      used.add(edgeKey(node, next));
      walk.push(next);
      node = next;
    }
    return walk.length >= 2 ? walk : null;
  }

  /**
   * Candidate strokes: each polyline the figure is stored as, the polylines run
   * together, and single-stroke walks begun from each loose end. People trace a
   * figure in whatever order feels natural, so the catalogue keeps several.
   */
  function strokeSets(lines) {
    var polys = [], i, j;
    for (i = 0; i < lines.length; i++) if (lines[i].p.length >= 2) polys.push(lines[i].p);
    if (!polys.length) return [];
    var out = polys.slice();
    if (polys.length > 1) {
      var cat = [];
      for (i = 0; i < polys.length; i++) {
        var p = polys[i];
        cat = cat.concat(cat.length && cat[cat.length - 1] === p[0] ? p.slice(1) : p);
      }
      out.push(cat);
      var adj = adjacency(polys);
      var ends = [], hubs = [];
      adj.forEach(function (v, k) { (v.length === 1 ? ends : hubs).push([k, v.length]); });
      hubs.sort(function (a, b) { return b[1] - a[1]; });
      var starts = ends.map(function (e) { return e[0]; }).slice(0, 7);
      for (i = 0; i < hubs.length && starts.length < 9; i++) starts.push(hubs[i][0]);
      var seen = {};
      for (i = 0; i < starts.length; i++) {
        var w = greedyWalk(adj, starts[i]);
        if (!w) continue;
        var k2 = w.join(',');
        if (seen[k2]) continue;
        seen[k2] = 1;
        out.push(w);
      }
      out.sort(function (a, b) { return b.length - a.length; });
    }
    var uniq = [], keys = {};
    for (i = 0; i < out.length; i++) {
      var kk = out[i].join(',');
      if (keys[kk]) continue;
      keys[kk] = 1;
      uniq.push(out[i]);
      if (uniq.length >= 12) break;
    }
    return uniq;
  }

  // ---------- catalogue ----------

  /**
   * @param {Array} figures  [{ iau, en, la, lines:[{t,p:[hip...]}], ... }]
   * @param {Function} lookup  hip -> [vx,vy,vz] | null
   */
  function buildCatalogue(figures, lookup) {
    var cat = [];
    for (var i = 0; i < figures.length; i++) {
      var f = figures[i];
      var members = new Set(), k, j;
      for (k = 0; k < f.lines.length; k++) for (j = 0; j < f.lines[k].p.length; j++) members.add(f.lines[k].p[j]);
      var vecs = [];
      members.forEach(function (h) { var v = lookup(h); if (v) vecs.push(v); });
      if (vecs.length < 2) continue;
      var geo = spread(vecs);
      var tp = tangentFrame(vecs);
      var descs = [];
      var strokes = strokeSets(f.lines);
      for (k = 0; k < strokes.length; k++) {
        var pts = [];
        for (j = 0; j < strokes[k].length; j++) {
          var v = lookup(strokes[k][j]);
          if (!v) continue;
          var pp = planeProject(v, tp);
          if (pp) pts.push(pp);
        }
        var z = normalise(resample(pts, N));
        var tf = turningProfile(pts, TURN_N);
        if (!z || !tf) continue;
        descs.push({ z: z, zr: reflectY(z), turn: tf, turnR: negate(tf), hips: strokes[k], pts: pts });
      }
      if (!descs.length) continue;
      cat.push({
        fig: f, members: members, descs: descs,
        centre: geo.centre, radiusDeg: geo.radiusDeg, frame: tp
      });
    }
    return cat;
  }

  /**
   * @param {number[]} hips  stars the user connected, in draw order
   * @param {Array} cat      output of buildCatalogue
   * @param {Function} lookup
   */
  function matchPath(hips, cat, lookup) {
    var vecs = [], i;
    for (i = 0; i < hips.length; i++) { var v = lookup(hips[i]); if (v) vecs.push(v); }
    if (vecs.length < 2) return [];
    var tp = tangentFrame(vecs);
    var pts = [];
    for (i = 0; i < vecs.length; i++) { var pp = planeProject(vecs[i], tp); if (pp) pts.push(pp); }
    var rev = pts.slice().reverse();
    var base = normalise(resample(pts, N));
    var baseTurn = turningProfile(pts, TURN_N);
    if (!base || !baseTurn) return [];
    var forms = [
      { z: base, turn: baseTurn, rev: false },
      { z: normalise(resample(rev, N)), turn: turningProfile(rev, TURN_N), rev: true }
    ].filter(function (f) { return f.z && f.turn; });
    var drawn = new Set(hips);
    var out = [];
    for (i = 0; i < cat.length; i++) {
      var c = cat[i], best = -1, keep = null;
      for (var d = 0; d < c.descs.length; d++) {
        var dc = c.descs[d];
        for (var f = 0; f < forms.length; f++) {
          var fm = forms[f];
          for (var m = 0; m < 2; m++) {
            var pro = correlate(fm.z, m ? dc.zr : dc.z);
            var td = turnDistance(fm.turn, m ? dc.turnR : dc.turn);
            var s = Math.exp(-td / TURN_SCALE) * Math.pow(Math.max(0, pro), PRO_POWER);
            if (m) s *= REFLECT_PENALTY;
            if (s > best) {
              best = s;
              keep = { desc: dc, form: fm, reflected: !!m, pro: pro, turnDist: td };
            }
          }
        }
      }
      if (!keep) continue;
      var shared = 0;
      drawn.forEach(function (h) { if (c.members.has(h)) shared++; });
      out.push({
        cat: c,
        shape: best,
        pro: keep.pro,
        turnDist: keep.turnDist,
        shared: shared,
        score: Math.min(0.999, best * (1 + Math.min(0.06, shared * 0.014))),
        desc: keep.desc,
        reversed: keep.form.rev,
        reflected: keep.reflected,
        userShape: keep.form.z
      });
    }
    out.sort(function (a, b) { return b.score - a.score; });
    return out;
  }

  /** Map a raw shape score onto a 0–1 "how close is this" reading. */
  function confidence(score) {
    var t = (score - 0.20) / 0.72;
    return Math.max(0.01, Math.min(1, Math.pow(Math.max(0, t), 0.8)));
  }

  root.Matcher = {
    N: N,
    tangentFrame: tangentFrame,
    planeProject: planeProject,
    spread: spread,
    resample: resample,
    normalise: normalise,
    reflectY: reflectY,
    reverseOrder: reverseOrder,
    correlate: correlate,
    alignment: alignment,
    turningProfile: turningProfile,
    turnDistance: turnDistance,
    strokeSets: strokeSets,
    buildCatalogue: buildCatalogue,
    matchPath: matchPath,
    confidence: confidence
  };
})(typeof window !== 'undefined' ? window : (typeof globalThis !== 'undefined' ? globalThis : this));
if (typeof module !== 'undefined' && module.exports) module.exports = (typeof globalThis !== 'undefined' ? globalThis : this).Matcher;
