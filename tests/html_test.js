/*
 * Verifies merge_tomo.html against merge_tomo.py.
 *
 * Both <script> blocks are pulled out of the page and run in a vm context with
 * a stub DOM, so the real merge path AND the real drawing path execute.  The
 * fake 2-D context flags any non-finite coordinate: in a browser those draw
 * nothing at all, which is exactly the plotting bug that would otherwise reach
 * the user unnoticed.
 *
 *   python tests/make_fixtures.py     (once, or after changing merge_tomo.py)
 *   node   tests/html_test.js
 */
const fs = require('fs');
const vm = require('vm');
const path = require('path');

const HERE = __dirname;
const FIX  = path.join(HERE, 'fixtures');
// optional argument: validate any built page, e.g. the published docs/index.html
const HTML = process.argv[2]
  ? path.resolve(process.argv[2])
  : path.join(HERE, '..', 'merge_tomo.html');
const fails = [];
const check = (name, cond, extra) => {
  console.log((cond ? '  PASS  ' : '  FAIL  ') + name + (extra ? '   ' + extra : ''));
  if (!cond) fails.push(name);
};
const REF = JSON.parse(fs.readFileSync(path.join(FIX, 'ref.json'), 'utf8'));

/* ------------------------------------------------ pull the scripts out */
const html = fs.readFileSync(HTML, 'utf8');
// executable blocks only: the published page also carries a JSON-LD block for
// search engines, which is data, not code
const blocks = [...html.matchAll(/<script([^>]*)>([\s\S]*?)<\/script>/g)]
  .filter(m => {
    const t = /type\s*=\s*["']([^"']+)/.exec(m[1]);
    return !t || /javascript/i.test(t[1]) || t[1] === 'module';
  })
  .map(m => m[2]);
check('page has exactly two executable script blocks', blocks.length === 2,
      `${path.basename(HTML)}: found ${blocks.length}`);

/* ------------------------------------------------------ stub DOM + canvas */
let drawCalls = 0, badCoord = null;
const ctx2d = new Proxy({}, {
  get(t, k){
    if (k === 'measureText') return () => ({width: 10});
    if (typeof k === 'string' &&
        /^(fillStyle|strokeStyle|lineWidth|font|textAlign|textBaseline|globalAlpha)$/.test(k))
      return t[k];
    return (...a) => {
      drawCalls++;
      for (const v of a)
        if (typeof v === 'number' && !Number.isFinite(v) && !badCoord) badCoord = String(k);
    };
  },
  set(t, k, v){ t[k] = v; return true; }
});
const els = new Map();
function mkEl(id){
  const kid = () => ({style: {}});
  const el = {
    id, value: '', textContent: '', innerHTML: '', className: '', disabled: false,
    min: '0.05', max: '10', step: '0.01', width: 1440, height: 1010, checked: false,
    scrollTop: 0, scrollHeight: 0, style: {}, files: [], dataset: {}, title: '',
    children: [], firstElementChild: kid(), lastElementChild: kid(),
    _cls: new Set(),
    addEventListener(){}, click(){}, remove(){},
    appendChild(c){ this.children.push(c); return c; },
    getContext(){ return ctx2d; }, toBlob(cb){ cb({}); }
  };
  el.classList = {
    add: c => el._cls.add(c), remove: c => el._cls.delete(c),
    toggle: (c, on) => { const want = on === undefined ? !el._cls.has(c) : !!on;
                         want ? el._cls.add(c) : el._cls.delete(c); return want; },
    contains: c => el._cls.has(c)
  };
  return el;
}
const sandbox = {
  console, performance: {now: () => Date.now()},
  setTimeout, clearTimeout, Math, Date, Number, Object, Array, Set, Map,
  Float64Array, isFinite, parseFloat, parseInt, JSON, Error, String, Infinity,
  URL: {createObjectURL: () => '', revokeObjectURL(){}}, Blob: function(){},
  document: {
    getElementById(id){ if (!els.has(id)) els.set(id, mkEl(id)); return els.get(id); },
    createElement(){ return mkEl('tmp'); },
    addEventListener(){}, body: {appendChild(){}}
  }
};
const ctxv = vm.createContext(sandbox);
vm.runInContext(blocks[0], ctxv, {filename: 'core.js'});
vm.runInContext(blocks[1], ctxv, {filename: 'view.js'});
const S = vm.runInContext('S', ctxv);      // top-level const: not a global property
const $ = id => ctxv.document.getElementById(id);
const load = f => ctxv.parseXYZ(fs.readFileSync(path.join(FIX, f), 'utf8'), f);

/* ============================================================== parsing */
console.log('\n--- parsing ---');
const A = load('west.xyz'), B = load('east.xyz');
check('west parsed', A.lon.length === 21291, `${A.lon.length} nodes`);
check('east parsed', B.lon.length === 21291);
check('comments and blank lines skipped',
      ctxv.parseXYZ('# hi\n\n> seg\n% c\n-120 40 3.5\n', 't').lon.length === 1);
check('comma-separated accepted', ctxv.parseXYZ('-120,40,3.5\n', 't').vel[0] === 3.5);
check('Fortran D exponent accepted',
      Math.abs(ctxv.parseXYZ('-120 40 3.5D0\n', 't').vel[0] - 3.5) < 1e-12);
for (const [txt, why] of [['-120 40\n', /expected 3 columns/],
                          ['-120 40 abc\n', /non-numeric/],
                          ['\n\n', /no data rows/]]){
  let ok = false;
  try { ctxv.parseXYZ(txt, 't'); } catch (e){ ok = why.test(e.message); }
  check(`rejects ${JSON.stringify(txt).slice(0,16)}`, ok);
}

/* ====================================================== kernel vs Python */
console.log('\n--- kernel maths matches merge_tomo.py ---');
{
  let worst = 0;
  for (const [lon, L, w] of REF.weights)
    worst = Math.max(worst, Math.abs(ctxv.westWeight(lon, L, -119, -110) - w));
  check('westWeight over 40 (lon,L) pairs', worst < 1e-14, `max diff ${worst.toExponential(2)}`);
}
{
  let worst = 0;
  for (const [lon, lat, L, ax, ay, bx, by, w] of REF.weights2d)
    worst = Math.max(worst, Math.abs(ctxv.blendWeight(lon, lat, L, [ax,ay], [bx,by]) - w));
  check('blendWeight over 2-D splits (E-W, N-S, diagonal)', worst < 1e-14,
        `${REF.weights2d.length} cases, max diff ${worst.toExponential(2)}`);
}
check('scalar and [lon,lat] centres agree on a shared latitude',
      [-130,-119,-114.5,-112,-100].every(x =>
        ctxv.westWeight(x, 4.52, -119, -110) ===
        ctxv.blendWeight(x, 43.5, 4.52, [-119,43.5], [-110,43.5])));
{
  let worst = 0;
  for (const [L, s] of REF.scales)
    worst = Math.max(worst, Math.abs(ctxv.handoverScale(L, -119, -110) - s));
  check('handoverScale', worst < 1e-14, `max diff ${worst.toExponential(2)}`);
}
const op = ctxv.overlapProjection(A, B, -119, -110);
check('overlap measured from shared nodes', op.n === 7701 &&
      Math.abs(op.tLo + 2.5) < 1e-9 && Math.abs(op.tHi - 2.5) < 1e-9,
      `${op.n} nodes, t ${op.tLo} … ${op.tHi}`);
const Ls = ctxv.suggestLProj(-119, -110, op.tLo, op.tHi);
check('suggested L agrees with Python to 1e-12',
      Math.abs(Ls - REF.suggestL) < 1e-12, `js ${Ls.toFixed(9)}  py ${REF.suggestL.toFixed(9)}`);
check('edges land exactly on the 90/10 target at that L',
      Math.abs(ctxv.westWeight(-117, Ls, -119, -110) - 0.9) < 1e-12 &&
      Math.abs(ctxv.westWeight(-112, Ls, -119, -110) - 0.1) < 1e-12);
check('midpoint outside the overlap -> null', ctxv.suggestLProj(-119, -110, 1, 5) === null);
check('weights stay finite and in [0,1] for L from 1e-6 to 1e6',
      [1e-6,1e-3,1,1e3,1e6].every(L => [-200,-119,-114.5,-110,-50].every(x => {
        const w = ctxv.westWeight(x, L, -119, -110);
        return Number.isFinite(w) && w >= 0 && w <= 1;
      })));
check('edgeQuality ok / loose',
      ctxv.edgeQuality(Ls, -119, -110, op.tLo, op.tHi).verdict === 'ok' &&
      ctxv.edgeQuality(Ls*2, -119, -110, op.tLo, op.tHi).verdict === 'loose');

/* ============================== output byte-identical to the Python tool */
console.log('\n--- merged output is byte-identical to merge_tomo.py ---');
function compare(a, b, L, cA, cB, pyFile, label){
  const mm = ctxv.gaussianMerge(a, b, L, cA, cB);
  const js = ctxv.formatXYZ(mm);
  const py = fs.readFileSync(path.join(FIX, pyFile), 'utf8');
  if (js === py){ check(label, true, `${mm.lon.length} lines identical`); return mm; }
  const x = js.split('\n'), y = py.split('\n');
  let first = -1, n = 0;
  for (let i = 0; i < Math.max(x.length, y.length); i++)
    if (x[i] !== y[i]){ n++; if (first < 0) first = i; }
  check(label, false, `${n} differing lines; first at ${first+1}\n    js "${x[first]}"\n    py "${y[first]}"`);
  return mm;
}
const m = compare(A, B, 4.52, -119, -110, 'py_default.xyz', 'east-west, L=4.52');
compare(A, B, 4.0, -118.5, -109.5, 'py_awk.xyz', 'east-west, the awk parameters');
const So = load('south.xyz'), No = load('north.xyz');
compare(So, No, 3.3, [-115,39], [-115,45], 'py_ns.xyz', 'north-south split');
const Dg1 = load('diagA.xyz'), Dg2 = load('diagB.xyz');
compare(Dg1, Dg2, 3.0, [-118,39], [-112,45], 'py_diag.xyz', 'diagonal split');

/* ==================================================== diagnostics match */
console.log('\n--- diagnostics ---');
const od = ctxv.overlapDiff(A, B, m);
const st = ctxv.stats(od.diff);
check('blended count', od.idx.length === 7701, String(od.idx.length));
check('A-only / B-only counts', m.src.filter(s=>s==='W').length === 13590 &&
      m.src.filter(s=>s==='E').length === 13590);
check('disagreement mean/rms', Math.abs(st.mean - REF.diffMean) < 1e-9 &&
      Math.abs(st.rms - REF.diffRms) < 1e-9,
      `mean ${st.mean.toFixed(4)} rms ${st.rms.toFixed(4)}`);
const seam = ctxv.seamReport(m, od, -119, -110);
check('seam rms at both overlap edges matches Python',
      Math.abs(seam.rmsW - REF.seamRmsW) < 1e-6 && Math.abs(seam.rmsE - REF.seamRmsE) < 1e-6,
      `${seam.rmsW.toFixed(4)} / ${seam.rmsE.toFixed(4)}`);
check('grid spacing detected', Math.abs(ctxv.spacing(m.lon) - 0.1) < 1e-9);
check('negative zero keeps its sign like Python', ctxv.fmt(-0,12,6) === '   -0.000000');
check('field widths', ctxv.fmt(-126,12,6) === ' -126.000000' &&
      ctxv.fmt(36,12,7) === '  36.0000000' && ctxv.fmt(4.11528444,14,8) === '    4.11528444');
check('diagnostic file has a header and one row per node',
      ctxv.formatDiag(m).trim().split('\n').length === m.lon.length + 1);

/* ================================================= node match tolerance */
console.log('\n--- node matching ---');
const base = load('base.xyz'), off = load('offset.xyz'), jit = load('jitter.xyz');
ctxv.setMatchTol(1e-4);
check('half-cell-offset grids share nothing at the default tolerance',
      ctxv.overlapProjection(base, off, -115, -114.9).n === 0);
check('nearestGap measures the actual grid offset',
      Math.abs(ctxv.nearestGap(base, off) - 0.125) < 1e-9,
      ctxv.nearestGap(base, off).toFixed(6));
check('same grid, coarser rounding: gap is tiny but above the default tolerance',
      ctxv.nearestGap(base, jit) > 1e-4 && ctxv.nearestGap(base, jit) < 1e-3,
      ctxv.nearestGap(base, jit).toExponential(2));
ctxv.setMatchTol(5e-4);
check('and a slightly wider tolerance pairs every node of it',
      ctxv.overlapProjection(base, jit, -115, -114.9).n === base.lon.length,
      `${ctxv.overlapProjection(base, jit, -115, -114.9).n} of ${base.lon.length}`);
ctxv.setMatchTol(1e-4);
check('tolerance <= 0 is rejected', (() => {
  try { ctxv.setMatchTol(0); return false; } catch { return true; }
})());
/* The bucket-snapping shortcut this replaced would drop pairs that straddle a
   boundary; a true radius must not, whatever the offset. */
{
  const n = 400, lonA = new Float64Array(n), latA = new Float64Array(n),
        lonB = new Float64Array(n), latB = new Float64Array(n);
  let seed = 1, rnd = () => (seed = (seed*16807) % 2147483647)/2147483647;
  for (let i=0;i<n;i++){
    lonA[i] = -120 + i*0.05; latA[i] = 40 + (i % 7)*0.05;
    lonB[i] = lonA[i] + (rnd()-0.5)*1.6e-5;      // well inside 1e-4
    latB[i] = latA[i] + (rnd()-0.5)*1.6e-5;
  }
  const all = Array.from({length:n}, (_, i) => i);
  const p = ctxv.matchNodes(lonA, latA, all, lonB, latB, all, 1e-4);
  check('radius matching finds every jittered pair (bucket snapping would not)',
        Array.from(p).every((v, i) => v === i), `${Array.from(p).filter((v,i)=>v===i).length}/${n}`);
}

/* ======================================================= the view layer */
console.log('\n--- the app, driven end to end ---');
S.A = A; S.B = B;
ctxv.initialiseFromData();
check('centres default to the map centroids on a shared latitude',
      +$('ax').value === -119 && +$('bx').value === -110 && $('ay').value === $('by').value,
      `A(${$('ax').value},${$('ay').value}) B(${$('bx').value},${$('by').value})`);
check('L auto-set to the floored recommendation',
      Math.abs(S.L - Math.floor(Ls*100)/100) < 1e-12, `L = ${S.L}`);
check('colour scale auto-filled', +$('vmin').value < +$('vmax').value,
      `${$('vmin').value} … ${$('vmax').value}`);
check('blend axis reported as east-west', /east–west/.test($('axisHint').textContent),
      $('axisHint').textContent);
check('verdict is the clean one', $('verdict').className === 'verdict ok');

drawCalls = 0; badCoord = null;
ctxv.runMerge();
check('merge + full QC render completes', drawCalls > 20000, `${drawCalls} canvas ops`);
check('no non-finite canvas coordinate', badCoord === null, badCoord || '');
check('app result matches the standalone merge',
      ctxv.formatXYZ(S.merged) === ctxv.formatXYZ(m));
check('downloads enabled', $('dlXYZ').disabled === false && $('dlPNG').disabled === false &&
      $('dlDiag').disabled === false);
console.log('\n--- colour palettes ---');
const PAL = vm.runInContext('PALETTES', ctxv);
const ORDER = vm.runInContext('PAL_ORDER', ctxv);
check('all 19 GMT palettes present, plus the uniform ones',
      ORDER.length >= 21 && ORDER.every(n => PAL[n]), `${Object.keys(PAL).length} palettes`);
{
  let bad = [];
  for (const n of Object.keys(PAL)){
    const a = PAL[n];
    if (a[0][0] !== 0 || Math.abs(a[a.length-1][0] - 1) > 1e-12) bad.push(n + ':range');
    for (let i=1;i<a.length;i++) if (a[i][0] < a[i-1][0]) bad.push(n + ':order');
    for (const s of a) if (!s.slice(1).every(v => Number.isFinite(v) && v >= 0 && v <= 255))
      bad.push(n + ':rgb');
  }
  check('every palette spans 0..1, is monotone, and stays in 0..255',
        bad.length === 0, bad.slice(0,4).join(' '));
}
{
  let bad = [];
  for (const n of Object.keys(PAL)){
    const f = ctxv.cmapFor(n), r = ctxv.cmapFor(n + '_r');
    for (const t of [0, 0.25, 0.5, 0.75, 1]){
      const c = f(t), c2 = r(1 - t);
      if (!c.every(Number.isFinite) || Math.max(...c.map((v,i)=>Math.abs(v-c2[i]))) > 1e-9)
        bad.push(n);
    }
    // out-of-range input must clamp, not produce NaN
    if (!f(-5).every(Number.isFinite) || !f(9).every(Number.isFinite)) bad.push(n + ':clamp');
  }
  check('_r really is the mirror image, and out-of-range input clamps',
        bad.length === 0, [...new Set(bad)].slice(0,4).join(' '));
}
check('unknown palette name falls back instead of throwing',
      ctxv.cmapFor('not_a_palette')(0.5).every(Number.isFinite));
check('css gradient built for every palette',
      Object.keys(PAL).every(n => /^linear-gradient\(to right,rgb\(/.test(ctxv.cssRamp(PAL[n]))));

console.log('\n--- GMT .cpt loading ---');
{
  const cpt = `# a GMT palette\n# COLOR_MODEL = RGB\n` +
              `-2000\t0\t0\t120\t-1000\t0\t120\t255\n` +
              `-1000 0 120 255 0 255 255 255\n` +
              `0 255 255 255 2000 200 0 0\n` +
              `B 0 0 0\nF 255 255 255\nN 128 128 128\n`;
  const a = ctxv.parseCPT(cpt, 'test.cpt');
  check('parses "z r g b z r g b" rows and skips B/F/N',
        a.length === 6 && a[0][0] === 0 && a[a.length-1][0] === 1,
        `${a.length} stops`);
  // z spans -2000..2000, so -1000 lands at 0.25 and 0 at 0.5
  check('z is normalised to 0..1 with the colours kept verbatim',
        a[0].slice(1).join() === '0,0,120' && a[a.length-1].slice(1).join() === '200,0,0' &&
        Math.abs(a[2][0] - 0.25) < 1e-12 && Math.abs(a[3][0] - 0.5) < 1e-12,
        a.map(s => s[0].toFixed(2)).join(' '));
  const b = ctxv.parseCPT('0 0/0/255 1 255/0/0\n', 'slash.cpt');
  check('parses the r/g/b slash form', b.length === 2 &&
        b[0].slice(1).join() === '0,0,255' && b[1].slice(1).join() === '255,0,0');
  const h = ctxv.parseCPT('# COLOR_MODEL = +HSV\n0 0 1 1 1 240 1 1\n', 'hsv.cpt');
  check('converts an HSV palette to RGB',
        Math.abs(h[0][1] - 255) < 1e-9 && Math.abs(h[1][3] - 255) < 1e-9,
        `${h[0].slice(1).map(Math.round)} → ${h[1].slice(1).map(Math.round)}`);
  for (const [txt, why] of [['# only a comment\n', /no usable colour rows/],
                            ['5 0 0 0 5 255 255 255\n', /does not increase/]]){
    let ok = false;
    try { ctxv.parseCPT(txt, 'bad.cpt'); } catch (e){ ok = why.test(e.message); }
    check(`rejects: ${why.source.slice(0,22)}`, ok);
  }
  // a real palette, loaded and drawn
  PAL['custom_test'] = a;
  badCoord = null; $('cmap').value = 'custom_test'; ctxv.drawQC();
  check('a loaded .cpt renders the QC figure', badCoord === null, badCoord || '');
}

console.log('\n--- every palette drawn for real ---');
{
  let bad = [];
  for (const n of ORDER){
    badCoord = null; $('cmap').value = n; ctxv.drawQC();
    if (badCoord) bad.push(n + ':' + badCoord);
    badCoord = null; $('cmap').value = n + '_r'; ctxv.drawQC();
    if (badCoord) bad.push(n + '_r:' + badCoord);
  }
  check(`all ${ORDER.length} palettes and their reverses render the full figure`,
        bad.length === 0, bad.slice(0,3).join(' '));
}
for (const cm of ['jet','viridis','gray','jet_r']){
  badCoord = null; $('cmap').value = cm; ctxv.drawQC();
  check(`redraw with ${cm}`, badCoord === null, badCoord || '');
}
for (const o of [0, 2, -3.5]){
  badCoord = null; $('poff').value = String(o); ctxv.drawQC();
  check(`profile offset ${o}° renders`, badCoord === null, badCoord || '');
}
badCoord = null; $('vmin').value = ''; $('vmax').value = ''; ctxv.drawQC();
check('empty colour-scale boxes fall back instead of drawing nothing', badCoord === null);
ctxv.autoScale();

badCoord = null; ctxv.setL(0.05); ctxv.runMerge();
check('very small L still renders', badCoord === null && S.merged.lon.length === 34881);
badCoord = null; ctxv.setL(60); ctxv.runMerge();
check('very large L still renders', badCoord === null);
check('over-smooth L turns the verdict red', $('verdict').className === 'verdict bad',
      $('verdict').className);
ctxv.recommend();

console.log('\n--- north-south split through the app ---');
S.A = So; S.B = No; S.merged = null;
ctxv.initialiseFromData();
check('north-south detected', /north–south/.test($('axisHint').textContent),
      $('axisHint').textContent);
check('centres share a longitude', $('ax').value === $('bx').value);
badCoord = null; ctxv.runMerge();
check('north-south merge renders', badCoord === null && S.od.idx.length > 0,
      `${S.od.idx.length} blended`);

console.log('\n--- degenerate inputs ---');
const far = ctxv.parseXYZ(Array.from({length: 60}, (_, i) =>
  `${-90 + (i % 6)*0.5} ${40 + Math.floor(i/6)*0.5} 4.2`).join('\n'), 'far');
S.A = A; S.B = far; S.merged = null;
ctxv.initialiseFromData();
badCoord = null; ctxv.runMerge();
check('disjoint maps concatenate without crashing',
      S.merged.lon.length === 21291 + 60 && S.od.idx.length === 0 && badCoord === null);
check('verdict explains there is nothing to blend', /concatenat/.test($('verdict').textContent),
      $('verdict').textContent.slice(0, 60));

S.A = base; S.B = off; S.merged = null;
ctxv.setMatchTol(1e-4); $('tol').value = '0.0001';
ctxv.initialiseFromData();
check('half-cell-offset grids are refused, not silently paired',
      +$('tol').value === 1e-4 && $('loadNote').className === 'verdict bad' &&
      /different grids/.test($('loadNote').textContent),
      $('loadNote').textContent.slice(0, 72));
badCoord = null; ctxv.runMerge();
check('and the merge still completes as a plain concatenation',
      S.od.idx.length === 0 && badCoord === null);

S.A = base; S.B = jit; S.merged = null;
ctxv.setMatchTol(1e-4); $('tol').value = '0.0001';
ctxv.initialiseFromData();
check('the same grid printed differently IS rescued',
      +$('tol').value > 1e-4 && +$('tol').value < 1e-3 &&
      $('loadNote').className === 'verdict warn',
      `tol now ${$('tol').value}`);
badCoord = null; ctxv.runMerge();
check('and then blends every shared node',
      S.od.idx.length === base.lon.length && badCoord === null,
      `${S.od.idx.length} blended`);

console.log('\n' + (fails.length ? 'FAILURES: ' + fails.join(' | ') : 'ALL CHECKS PASSED'));
process.exit(fails.length ? 1 : 0);
