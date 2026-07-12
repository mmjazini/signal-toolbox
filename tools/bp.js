/* ==========================================================================
   bp.js — Oscillometric BP reasoner (rigorous edition).

   Drop-in, self-contained. Lives under section #bp on toolbox.html and
   replaces the simpler single-canvas envelope demo with a proper
   physics-grounded two-panel simulation.

   Panel 1 (#bp-tool-cuff): the raw cuff pressure trace during deflation.
     A slow ramp from ~180 mmHg to ~30 mmHg over ~30 s with the small
     (~1-4 mmHg pk-pk) oscillometric pulses superimposed. This is what
     the actual transducer sees before any processing.

   Panel 2 (#bp-tool-envelope): the extracted oscillogram envelope —
     pulse pk-pk amplitude vs cuff pressure — with markers for true
     SBP/MAP/DBP and the chosen estimator's predictions.

   The envelope is derived from a sigmoidal arterial volume curve
       V(Ptm) = a / (1 + exp(-(Ptm - b)/c))
   where Ptm = arterial_pressure(t) - cuff_pressure. Pulse amplitude at
   cuff pressure Pc is V(Ptm_systolic) - V(Ptm_diastolic), which produces
   the asymmetric bell shape Mahdi's 2025 Sci Reports paper exploits.

   Estimators implemented: fixed-ratio, max-amplitude (MAP-only),
   variable-ratio (Dhamotharan, Jazini et al. 2025), and a derivative
   method using slope-extrema on each side of the envelope peak.

   No external libraries. ~900 lines.
   ========================================================================== */

(function bpTool() {

  // Theme-aware canvas colors (reads data-theme attr cached on root — no DOM query overhead)
  const _cc = (() => {
    const r = document.documentElement;
    return {
      isLight:  () => r.dataset.theme === 'light',
      grid:     () => r.dataset.theme === 'light' ? 'rgba(0,0,0,0.06)'  : 'rgba(255,255,255,0.06)',
      axis:     () => r.dataset.theme === 'light' ? 'rgba(0,0,0,0.10)'  : 'rgba(255,255,255,0.10)',
      frame:    () => r.dataset.theme === 'light' ? 'rgba(0,0,0,0.12)'  : 'rgba(255,255,255,0.12)',
      label:    () => r.dataset.theme === 'light' ? 'rgba(0,0,0,0.55)'  : 'rgba(255,255,255,0.55)',
      labelDim: () => r.dataset.theme === 'light' ? 'rgba(0,0,0,0.35)'  : 'rgba(255,255,255,0.35)',
      ghost:    () => r.dataset.theme === 'light' ? 'rgba(0,0,0,0.15)'  : 'rgba(255,255,255,0.12)',
      dot:      () => r.dataset.theme === 'light' ? '#1a1a1f'           : '#e8e8ed',
    };
  })();
  'use strict';

  // ----- bail if this page doesn't host the tool ---------------------------
  const cuffEl = document.getElementById('bp-tool-cuff');
  if (!cuffEl) return;
  const envEl  = document.getElementById('bp-tool-envelope');
  if (!envEl)  return;

  // ----- DPR-aware canvas helper ------------------------------------------
  const DPR = Math.min(window.devicePixelRatio || 1, 2);

  function makeCanvas(canvasEl) {
    const ctx = canvasEl.getContext('2d');
    let W = 0, H = 0;
    function resize() {
      // Honour the CSS-set width and height of the canvas itself rather
      // than its parent, so tool-specific heights set in styles.css apply.
      const r = canvasEl.getBoundingClientRect();
      W = Math.max(120, Math.floor(r.width));
      H = Math.max(80,  Math.floor(r.height));
      canvasEl.width  = Math.round(W * DPR);
      canvasEl.height = Math.round(H * DPR);
      ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
    }
    resize();
    return {
      ctx,
      get W() { return W; },
      get H() { return H; },
      canvas: canvasEl,
      resize,
    };
  }

  const cuffCV = makeCanvas(cuffEl);
  const envCV  = makeCanvas(envEl);

  // ----- DOM controls ------------------------------------------------------
  const $ = (id) => document.getElementById(id);

  const inputs = {
    sbp:    $('bp-sbp'),
    dbp:    $('bp-dbp'),
    hr:     $('bp-hr'),
    rate:   $('bp-rate'),
    cstiff: $('bp-c'),
    shape:  $('bp-shape'),
    noise:  $('bp-noise'),
    method: $('bp-method'),
    rs:     $('bp-rs'),
    rd:     $('bp-rd'),
    play:   $('bp-play'),
  };

  const labels = {
    sbp:    $('bp-sbp-val'),
    dbp:    $('bp-dbp-val'),
    hr:     $('bp-hr-val'),
    rate:   $('bp-rate-val'),
    cstiff: $('bp-c-val'),
    shape:  $('bp-shape-val'),
    noise:  $('bp-noise-val'),
    rs:     $('bp-rs-val'),
    rd:     $('bp-rd-val'),
  };

  const readouts = {
    trueSbp:  $('bp-true-sbp'),
    trueMap:  $('bp-true-map'),
    trueDbp:  $('bp-true-dbp'),
    estSbp:   $('bp-est-sbp'),
    estMap:   $('bp-est-map'),
    estDbp:   $('bp-est-dbp'),
    errSbp:   $('bp-err-sbp'),
    errMap:   $('bp-err-map'),
    errDbp:   $('bp-err-dbp'),
    peakAmp:  $('bp-peak-amp'),
    peakLoc:  $('bp-peak-loc'),
    conf:     $('bp-conf'),
    method:   $('bp-method-note'),
  };

  // The estimator-specific control row(s) get shown/hidden based on the
  // <select> value. We tag them with data-method on the HTML side, but
  // since we can't modify HTML, we just toggle the inline display of the
  // ratio inputs' parent .tool-row by id.
  const ratioRow = $('bp-ratio-row');

  // ----- physics & simulation parameters -----------------------------------
  const FS_CUFF        = 100;       // simulated sample rate of the cuff trace, Hz
  const P_START        = 180;       // cuff pressure at start of deflation, mmHg
  const P_END          = 30;        // cuff pressure at end of deflation,   mmHg
  const ENV_NN         = 401;       // envelope sample count
  const ENV_PMIN       = 30;        // envelope x-axis min
  const ENV_PMAX       = 200;       // envelope x-axis max

  // ----- Box-Muller Gaussian -----------------------------------------------
  let _spareGauss = null;
  function gauss() {
    if (_spareGauss !== null) {
      const v = _spareGauss; _spareGauss = null; return v;
    }
    let u1 = 0, u2 = 0;
    while (u1 === 0) u1 = Math.random();
    while (u2 === 0) u2 = Math.random();
    const mag = Math.sqrt(-2.0 * Math.log(u1));
    const z0 = mag * Math.cos(2.0 * Math.PI * u2);
    const z1 = mag * Math.sin(2.0 * Math.PI * u2);
    _spareGauss = z1;
    return z0;
  }

  // ----- shape function for the arterial pressure pulse --------------------
  // Returns a normalized 0..1 waveform over phase in [0,1] with a tunable
  // rise-to-decay ratio. shape = 0 -> symmetric raised-cosine pulse;
  // shape = 1 -> sharp upstroke, slow falloff (more like aortic pressure).
  function pulseShape(phase, shape) {
    if (phase < 0 || phase >= 1) return 0;
    // rise fraction goes from 0.5 at shape=0 down to 0.18 at shape=1
    const rise = 0.5 - 0.32 * shape;
    let y;
    if (phase < rise) {
      // cosine rise from 0 to 1
      const u = phase / rise;
      y = 0.5 - 0.5 * Math.cos(Math.PI * u);
    } else {
      // cosine decay from 1 back to 0
      const u = (phase - rise) / (1 - rise);
      y = 0.5 + 0.5 * Math.cos(Math.PI * u);
    }
    // Tail flattens slightly to mimic dicrotic shoulder when shape > 0.6
    if (shape > 0.6 && phase > rise + 0.25) {
      const dq = (phase - rise - 0.25) * 4.0;
      const dn = Math.exp(-dq * dq) * 0.06 * (shape - 0.6) / 0.4;
      y += dn * (1 - phase);
    }
    return Math.max(0, Math.min(1, y));
  }

  // ----- arterial pressure waveform ---------------------------------------
  // Returns mmHg given seconds and (sbp, dbp, hr_bpm, shape).
  function arterialP(t, sbp, dbp, hr, shape) {
    const T = 60 / hr;
    const phase = (t % T) / T;
    return dbp + (sbp - dbp) * pulseShape(phase, shape);
  }

  // ----- arterial volume curve (sigmoidal compliance) ---------------------
  // V(Ptm) = a / (1 + exp(-(Ptm - b)/c)). Defaults a=1, b=0, c per slider.
  function volume(Ptm, c) {
    return 1.0 / (1.0 + Math.exp(-Ptm / c));
  }

  // ----- analytic envelope point ------------------------------------------
  // For a single beat at cuff pressure Pc, the pulse amplitude under the
  // cuff is V(P_sys - Pc) - V(P_dia - Pc). This is the *clean* envelope.
  function envelopePoint(Pc, sbp, dbp, c) {
    return volume(sbp - Pc, c) - volume(dbp - Pc, c);
  }

  // ----- compute the analytic envelope across cuff pressures --------------
  function computeEnvelope(sbp, dbp, c) {
    const env = new Float64Array(ENV_NN);
    const xs  = new Float64Array(ENV_NN);
    let peak = 0, peakIdx = 0;
    for (let i = 0; i < ENV_NN; i++) {
      const Pc = ENV_PMIN + (ENV_PMAX - ENV_PMIN) * (i / (ENV_NN - 1));
      xs[i]  = Pc;
      env[i] = envelopePoint(Pc, sbp, dbp, c);
      if (env[i] > peak) { peak = env[i]; peakIdx = i; }
    }
    // Reference peak amplitude for a normotensive 120/80 envelope at the SAME
    // arterial stiffness c. Used by the variable-ratio estimator so its
    // amplitude term is dimensionless (peak / refPeak) rather than comparing an
    // arbitrary-unit sigmoid peak (≪1) against a hard-coded 1.0.
    const refPeak = maxEnvelopeAmp(120, 80, c);
    return { xs, env, peak, peakIdx, peakP: xs[peakIdx], refPeak };
  }

  // Peak (over cuff pressure) of the analytic envelope for a given sbp/dbp/c.
  function maxEnvelopeAmp(sbp, dbp, c) {
    let mx = 0;
    for (let i = 0; i < ENV_NN; i++) {
      const Pc = ENV_PMIN + (ENV_PMAX - ENV_PMIN) * (i / (ENV_NN - 1));
      const v = envelopePoint(Pc, sbp, dbp, c);
      if (v > mx) mx = v;
    }
    return mx;
  }

  // ----- simulate the full cuff deflation trace ---------------------------
  // Returns { t, Pc, raw } sampled at FS_CUFF.
  // raw = Pc + (oscillometric AC component) + Gaussian noise.
  // The oscillometric AC component for a beat at cuff pressure Pc is:
  //   delta(t) = K * [ V(P_art(t) - Pc, c) - V(<V>, c) ]
  // i.e. instantaneous deviation from the mean. K converts a.u. volume
  // into mmHg pressure deviation under the cuff. We pick K so the peak
  // pulse is ~3 mmHg pk-pk near MAP — matches a real arm cuff.
  function simulateDeflation(p) {
    const dur = (P_START - P_END) / p.rate; // seconds
    const N = Math.max(8, Math.floor(dur * FS_CUFF));
    const t = new Float64Array(N);
    const Pc = new Float64Array(N);
    const raw = new Float64Array(N);
    const ac = new Float64Array(N);
    const vol = new Float64Array(N);   // arterial volume oscillation (a.u.), = V - mean
    const dt = 1 / FS_CUFF;
    // Calibrate K so the maximum pk-pk pulse is ~3 mmHg.
    // The oscillometric envelope peaks at the cuff pressure where mean
    // transmural pressure = 0 and compliance is maximal — i.e. MAP. In this
    // symmetric sigmoid model that peak is exactly Pc = (SBP+DBP)/2 (holds for
    // any c, since V' is symmetric about Ptm=0). -> K = 3 / peak.
    const map = (p.sbp + p.dbp) / 2;
    const peakAmp = envelopePoint(map, p.sbp, p.dbp, p.cstiff);
    const K = 3.0 / Math.max(0.01, peakAmp);

    for (let i = 0; i < N; i++) {
      const ti = i * dt;
      t[i] = ti;
      Pc[i] = P_START - p.rate * ti;
      const Part = arterialP(ti, p.sbp, p.dbp, p.hr, p.shape);
      // Instantaneous volume relative to the mean of V over the cardiac
      // cycle. We approximate the mean as the volume at MAP for speed.
      const v   = volume(Part - Pc[i], p.cstiff);
      const vbar = volume(map - Pc[i], p.cstiff);
      vol[i] = v - vbar;                 // raw volume oscillation (a.u.)
      ac[i] = K * vol[i];
      raw[i] = Pc[i] + ac[i] + p.noise * gauss();
    }
    return { t, Pc, raw, ac, vol, dt, N, dur };
  }

  // ----- estimators --------------------------------------------------------
  // All estimators consume the analytic envelope. They return
  //   { sbp, map, dbp } in mmHg, where any field may be NaN if the method
  // doesn't produce that estimate.
  function estimateFixedRatio(envObj, rs, rd) {
    const { xs, env, peak, peakIdx } = envObj;
    const NN = env.length;
    let estSbp = NaN, estDbp = NaN;
    // SBP: scan from peak toward higher cuff pressures, find first
    // crossing of rs * peak.
    for (let i = peakIdx; i < NN - 1; i++) {
      if (env[i] >= rs * peak && env[i + 1] < rs * peak) {
        // linear interpolate
        const frac = (rs * peak - env[i + 1]) / (env[i] - env[i + 1] + 1e-12);
        estSbp = xs[i + 1] + frac * (xs[i] - xs[i + 1]);
        break;
      }
    }
    // DBP: scan from peak toward lower cuff pressures.
    for (let i = peakIdx; i > 0; i--) {
      if (env[i] >= rd * peak && env[i - 1] < rd * peak) {
        const frac = (rd * peak - env[i - 1]) / (env[i] - env[i - 1] + 1e-12);
        estDbp = xs[i - 1] + frac * (xs[i] - xs[i - 1]);
        break;
      }
    }
    return { sbp: estSbp, map: xs[peakIdx], dbp: estDbp };
  }

  function estimateMaxAmplitude(envObj) {
    return { sbp: NaN, map: envObj.peakP, dbp: NaN };
  }

  // Variable-ratio estimator (echoing Dhamotharan, Jazini et al. Sci
  // Reports 2025): SBP ratio scales with peak amplitude relative to a
  // population reference; DBP ratio moves the opposite way. The amplitude
  // term must be dimensionless — we normalize this patient's envelope peak
  // by the peak of a normotensive 120/80 envelope at the same arterial
  // stiffness (envObj.refPeak). So A_ref = 1.0 means "same as reference",
  // and a wider pulse pressure (larger relative peak) shifts the ratios.
  function estimateVariableRatio(envObj) {
    const A_ref = 1.0;
    const A = envObj.refPeak > 1e-9 ? (envObj.peak / envObj.refPeak) : 1.0;
    let r_s = 0.55 + 0.05 * (A - A_ref) / A_ref;
    let r_d = 0.75 - 0.05 * (A - A_ref) / A_ref;
    if (r_s < 0.45) r_s = 0.45;
    if (r_s > 0.65) r_s = 0.65;
    if (r_d < 0.65) r_d = 0.65;
    if (r_d > 0.85) r_d = 0.85;
    const r = estimateFixedRatio(envObj, r_s, r_d);
    r.rs = r_s; r.rd = r_d;
    return r;
  }

  // Derivative method: SBP at the cuff pressure where dE/dPc has its
  // *minimum* (steepest fall) on the systolic side of the peak;
  // DBP where dE/dPc has its maximum (steepest rise) on the diastolic
  // side. These slope-extrema bracket the inflection points of the
  // sigmoidal envelope, which empirically track SBP and DBP well.
  function estimateDerivative(envObj) {
    const { xs, env, peakIdx } = envObj;
    const NN = env.length;
    // 1st derivative
    const d = new Float64Array(NN);
    const dx = xs[1] - xs[0];
    for (let i = 1; i < NN - 1; i++) {
      d[i] = (env[i + 1] - env[i - 1]) / (2 * dx);
    }
    d[0] = d[1]; d[NN - 1] = d[NN - 2];
    // Systolic side (i > peakIdx): d is negative; we want its minimum.
    let minD = Infinity, minI = peakIdx;
    for (let i = peakIdx + 1; i < NN; i++) {
      if (d[i] < minD) { minD = d[i]; minI = i; }
    }
    // Diastolic side (i < peakIdx): d is positive; we want its maximum.
    let maxD = -Infinity, maxI = peakIdx;
    for (let i = peakIdx - 1; i >= 0; i--) {
      if (d[i] > maxD) { maxD = d[i]; maxI = i; }
    }
    return { sbp: xs[minI], map: xs[peakIdx], dbp: xs[maxI] };
  }

  // ----- confidence rating ------------------------------------------------
  // Heuristic: confidence rises with envelope sharpness (higher peak)
  // and symmetry (similar widths above and below the peak at half-max).
  function confidenceRating(envObj) {
    const { xs, env, peak, peakIdx } = envObj;
    const NN = env.length;
    if (peak <= 0) return { tag: 'low', score: 0 };
    const half = 0.5 * peak;
    let upI = NN - 1, downI = 0;
    for (let i = peakIdx; i < NN; i++) if (env[i] < half) { upI = i; break; }
    for (let i = peakIdx; i >= 0; i--) if (env[i] < half) { downI = i; break; }
    const widthUp   = xs[upI]   - xs[peakIdx];
    const widthDown = xs[peakIdx] - xs[downI];
    const total = widthUp + widthDown;
    if (total <= 0) return { tag: 'low', score: 0 };
    const sym = 1 - Math.abs(widthUp - widthDown) / total; // 0..1
    // sharpness: peak / median envelope (avoids div-by-zero)
    let sumEnv = 0;
    for (let i = 0; i < NN; i++) sumEnv += env[i];
    const mean = sumEnv / NN;
    const sharp = Math.min(1, (peak - mean) / (peak + 1e-9) * 1.5);
    const score = 0.6 * sym + 0.4 * sharp;
    let tag = 'low';
    if (score > 0.55) tag = 'medium';
    if (score > 0.75) tag = 'high';
    return { tag, score };
  }

  // ----- read & validate slider state -------------------------------------
  function getParams() {
    const sbp   = +inputs.sbp.value;
    let   dbp   = +inputs.dbp.value;
    if (dbp > sbp - 10) {
      dbp = sbp - 10;
      inputs.dbp.value = dbp;
    }
    return {
      sbp,
      dbp,
      hr:     +inputs.hr.value,
      rate:   +inputs.rate.value,
      cstiff: +inputs.cstiff.value,
      shape:  +inputs.shape.value,
      noise:  +inputs.noise.value,
      method: inputs.method.value,
      rs:     +inputs.rs.value,
      rd:     +inputs.rd.value,
    };
  }

  function syncLabels(p) {
    if (labels.sbp)    labels.sbp.textContent    = p.sbp + ' mmHg';
    if (labels.dbp)    labels.dbp.textContent    = p.dbp + ' mmHg';
    if (labels.hr)     labels.hr.textContent     = p.hr + ' bpm';
    if (labels.rate)   labels.rate.textContent   = p.rate.toFixed(1) + ' mmHg/s';
    if (labels.cstiff) labels.cstiff.textContent = p.cstiff.toFixed(0);
    if (labels.shape)  labels.shape.textContent  = p.shape.toFixed(2);
    if (labels.noise)  labels.noise.textContent  = p.noise.toFixed(2) + ' mmHg';
    if (labels.rs)     labels.rs.textContent     = p.rs.toFixed(2);
    if (labels.rd)     labels.rd.textContent     = p.rd.toFixed(2);
  }

  function setMethodNote(method) {
    if (!readouts.method) return;
    let txt = '';
    switch (method) {
      case 'fixed':
        txt = 'Fixed-ratio: SBP = cuff pressure where envelope crosses rs · peak; DBP at rd · peak. The ratios are population averages — they only fit one patient by accident.';
        break;
      case 'maxamp':
        txt = 'Maximum-amplitude: MAP only. The peak of the envelope is the cleanest landmark there is, but you get no SBP/DBP without further assumptions.';
        break;
      case 'variable':
        txt = 'Variable-ratio (Dhamotharan, Jazini et al., Sci Reports 2025): rs and rd shift with peak amplitude. Reduces error at extreme BPs vs fixed ratios.';
        break;
      case 'deriv':
        txt = 'Derivative method: SBP at the steepest-fall slope on the systolic side of the envelope; DBP at the steepest-rise slope on the diastolic side. No ratios.';
        break;
    }
    readouts.method.textContent = txt;
  }

  function showRatioRow(method) {
    if (!ratioRow) return;
    ratioRow.style.display = (method === 'fixed') ? '' : 'none';
  }

  // ----- core update -------------------------------------------------------
  let cachedSim = null;
  let cachedEnv = null;
  let cachedEst = null;
  let cachedConf = null;
  let cachedParams = null;
  let cursorIdx = 0;       // current sample index along simulated trace
  let playing = false;
  let lastPlayT = 0;

  function recompute() {
    const p = getParams();
    syncLabels(p);
    setMethodNote(p.method);
    showRatioRow(p.method);

    cachedParams = p;
    cachedSim = simulateDeflation(p);
    cachedEnv = computeEnvelope(p.sbp, p.dbp, p.cstiff);

    let est;
    switch (p.method) {
      case 'maxamp':   est = estimateMaxAmplitude(cachedEnv); break;
      case 'variable': est = estimateVariableRatio(cachedEnv); break;
      case 'deriv':    est = estimateDerivative(cachedEnv); break;
      case 'fixed':
      default:         est = estimateFixedRatio(cachedEnv, p.rs, p.rd); break;
    }
    cachedEst = est;
    cachedConf = confidenceRating(cachedEnv);

    // readouts. True MAP is defined here as the oscillometric landmark: the
    // cuff pressure of maximum pulse amplitude = zero mean transmural pressure
    // = (SBP+DBP)/2 in this compliance model (NOT the clinical (SBP+2DBP)/3).
    const trueMap = (p.sbp + p.dbp) / 2;
    if (readouts.trueSbp) readouts.trueSbp.textContent = p.sbp;
    if (readouts.trueMap) readouts.trueMap.textContent = trueMap.toFixed(0);
    if (readouts.trueDbp) readouts.trueDbp.textContent = p.dbp;
    if (readouts.estSbp)  readouts.estSbp.textContent  = isFinite(est.sbp) ? est.sbp.toFixed(0) : '—';
    if (readouts.estMap)  readouts.estMap.textContent  = isFinite(est.map) ? est.map.toFixed(0) : '—';
    if (readouts.estDbp)  readouts.estDbp.textContent  = isFinite(est.dbp) ? est.dbp.toFixed(0) : '—';
    if (readouts.errSbp)  readouts.errSbp.textContent  = isFinite(est.sbp) ? (est.sbp - p.sbp).toFixed(1) : '—';
    if (readouts.errMap)  readouts.errMap.textContent  = isFinite(est.map) ? (est.map - trueMap).toFixed(1) : '—';
    if (readouts.errDbp)  readouts.errDbp.textContent  = isFinite(est.dbp) ? (est.dbp - p.dbp).toFixed(1) : '—';
    if (readouts.peakAmp) readouts.peakAmp.textContent = cachedEnv.peak.toFixed(3);
    if (readouts.peakLoc) readouts.peakLoc.textContent = cachedEnv.peakP.toFixed(0) + ' mmHg';
    if (readouts.conf) {
      readouts.conf.textContent = cachedConf.tag;
      readouts.conf.dataset.level = cachedConf.tag;
    }

    // reset cursor if it's beyond the new trace
    if (cursorIdx >= cachedSim.N) cursorIdx = 0;
    redraw();
  }

  // ----- drawing -----------------------------------------------------------
  function drawGrid(ctx, W, H, stepX, stepY) {
    ctx.save();
    ctx.strokeStyle = _cc.grid();
    ctx.lineWidth = 1;
    for (let x = 0; x < W; x += stepX) {
      ctx.beginPath();
      ctx.moveTo(Math.round(x) + 0.5, 0);
      ctx.lineTo(Math.round(x) + 0.5, H);
      ctx.stroke();
    }
    for (let y = 0; y < H; y += stepY) {
      ctx.beginPath();
      ctx.moveTo(0, Math.round(y) + 0.5);
      ctx.lineTo(W, Math.round(y) + 0.5);
      ctx.stroke();
    }
    ctx.restore();
  }

  function drawCuffTrace() {
    if (!cachedSim) return;
    const { ctx, W, H } = cuffCV;
    const sim = cachedSim;
    const padL = 38, padR = 12, padT = 12, padB = 26;
    const plotW = W - padL - padR;
    const plotH = H - padT - padB;

    ctx.clearRect(0, 0, W, H);
    drawGrid(ctx, W, H, Math.max(40, W / 14), Math.max(28, H / 5));

    // y-axis: cuff pressure 0..200 mmHg
    const yMin = 0, yMax = 200;
    const yToPx = (y) => padT + plotH * (1 - (y - yMin) / (yMax - yMin));

    // axis labels
    ctx.fillStyle = _cc.label();
    ctx.font = '11px JetBrains Mono, Consolas, monospace';
    [0, 50, 100, 150, 200].forEach(p => {
      const yy = yToPx(p);
      ctx.fillText(p.toString().padStart(3, ' '), 4, yy + 4);
      ctx.strokeStyle = _cc.grid();
      ctx.beginPath(); ctx.moveTo(padL, yy + 0.5); ctx.lineTo(W - padR, yy + 0.5); ctx.stroke();
    });
    ctx.fillStyle = _cc.labelDim();
    ctx.fillText('mmHg', 4, padT - 2);

    // x-axis: time 0..dur seconds
    const tMax = sim.dur;
    const xToPx = (t) => padL + plotW * (t / tMax);
    [0, tMax * 0.25, tMax * 0.5, tMax * 0.75, tMax].forEach(tt => {
      const xx = xToPx(tt);
      ctx.fillText(tt.toFixed(1) + 's', xx + 2, H - 6);
    });

    // baseline cuff ramp (faint)
    ctx.strokeStyle = _cc.ghost();
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(xToPx(0),    yToPx(sim.Pc[0]));
    ctx.lineTo(xToPx(tMax), yToPx(sim.Pc[sim.N - 1]));
    ctx.stroke();
    ctx.setLineDash([]);

    // raw trace (cuff + AC + noise)
    ctx.strokeStyle = 'rgba(255,255,255,0.7)';
    ctx.lineWidth = 1.0;
    ctx.beginPath();
    // Subsample for speed: at most ~2*W points across the canvas.
    const stride = Math.max(1, Math.floor(sim.N / (W * 2)));
    for (let i = 0; i < sim.N; i += stride) {
      const x = xToPx(sim.t[i]);
      const y = yToPx(sim.raw[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // ---- Arterial VOLUME oscillations (the pulse that builds the envelope) --
    // Draw the raw volume pulsation V(Part−Pc) as its own trace, riding along
    // the falling cuff ramp and amplified so it's actually visible on the
    // 0–200 mmHg axis. It swells to a maximum near MAP (zero transmural
    // pressure / peak compliance) and shrinks toward SBP and DBP — this is the
    // oscillogram that panel 2 turns into SBP/MAP/DBP.
    let volMax = 1e-9;
    for (let i = 0; i < sim.N; i++) { const a = Math.abs(sim.vol[i]); if (a > volMax) volMax = a; }
    // Amplify so the largest oscillation spans ~30 mmHg of the axis (visible
    // but not overwhelming), centred on the cuff ramp.
    const volGain = 30 / volMax;
    ctx.strokeStyle = 'rgba(10,132,255,0.85)';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let i = 0; i < sim.N; i += stride) {
      const x = xToPx(sim.t[i]);
      const y = yToPx(sim.Pc[i] + sim.vol[i] * volGain);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // moving cursor at current cuff pressure
    if (cursorIdx >= 0 && cursorIdx < sim.N) {
      const cx = xToPx(sim.t[cursorIdx]);
      const cy = yToPx(sim.Pc[cursorIdx]);
      ctx.strokeStyle = 'rgba(255,159,10,0.9)';
      ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(cx, padT); ctx.lineTo(cx, H - padB); ctx.stroke();
      ctx.fillStyle = '#ff9f0a';
      ctx.beginPath(); ctx.arc(cx, cy, 3.5, 0, Math.PI * 2); ctx.fill();
      // readout box
      const txt = sim.Pc[cursorIdx].toFixed(0) + ' mmHg';
      ctx.font = '11px JetBrains Mono, Consolas, monospace';
      const tw = ctx.measureText(txt).width;
      const bx = Math.min(W - padR - tw - 8, cx + 6);
      const by = padT + 4;
      ctx.fillStyle = 'rgba(0,0,0,0.55)';
      ctx.fillRect(bx, by, tw + 8, 16);
      ctx.fillStyle = '#ff9f0a';
      ctx.fillText(txt, bx + 4, by + 12);
    }

    // title strip
    ctx.fillStyle = _cc.label();
    ctx.font = '11px JetBrains Mono, Consolas, monospace';
    ctx.fillText('cuff pressure during deflation (raw transducer signal)', padL, padT - 2);
    // legend for the volume-oscillation overlay
    ctx.fillStyle = 'rgba(10,132,255,0.95)';
    ctx.fillText('— arterial volume oscillation (a.u., amplified)', W - padR - 250, padT - 2);
  }

  function drawEnvelope() {
    if (!cachedEnv || !cachedParams) return;
    const { ctx, W, H } = envCV;
    const env = cachedEnv;
    const p   = cachedParams;
    const est = cachedEst;

    // padT gives clear headroom below the panel title so the envelope PEAK
    // (near the top of the plot) is never hidden behind the title text.
    const padL = 38, padR = 12, padT = 26, padB = 28;
    const plotW = W - padL - padR;
    const plotH = H - padT - padB;

    ctx.clearRect(0, 0, W, H);
    drawGrid(ctx, W, H, Math.max(40, W / 14), Math.max(28, H / 5));

    // y axis: amplitude 0..max. 1.18 headroom keeps the peak marker + its
    // implied top well inside the plot area on the shorter panel.
    const yMax = Math.max(0.01, env.peak * 1.18);
    const yToPx = (y) => padT + plotH * (1 - y / yMax);

    // x axis: cuff pressure ENV_PMIN..ENV_PMAX
    const xToPx = (Pc) => padL + plotW * (Pc - ENV_PMIN) / (ENV_PMAX - ENV_PMIN);

    // y-axis ticks
    ctx.fillStyle = _cc.label();
    ctx.font = '11px JetBrains Mono, Consolas, monospace';
    for (let k = 0; k <= 4; k++) {
      const v = (k / 4) * yMax;
      const yy = yToPx(v);
      ctx.fillText(v.toFixed(2), 4, yy + 4);
      ctx.strokeStyle = _cc.grid();
      ctx.beginPath(); ctx.moveTo(padL, yy + 0.5); ctx.lineTo(W - padR, yy + 0.5); ctx.stroke();
    }
    ctx.fillText('A (a.u.)', 4, padT - 4);

    // x-axis ticks
    [50, 80, 100, 120, 150, 180].forEach(Pc => {
      if (Pc >= ENV_PMIN && Pc <= ENV_PMAX) {
        const xx = xToPx(Pc);
        ctx.fillStyle = _cc.label();
        ctx.fillText(Pc + '', xx - 8, H - 8);
      }
    });
    ctx.fillStyle = _cc.labelDim();
    ctx.fillText('cuff pressure (mmHg)', W - padR - 130, H - 8);

    // threshold horizontal lines for fixed-ratio
    if (p.method === 'fixed') {
      const thrS = p.rs * env.peak;
      const thrD = p.rd * env.peak;
      ctx.strokeStyle = 'rgba(255,55,95,0.3)';
      ctx.setLineDash([3, 4]);
      ctx.beginPath(); ctx.moveTo(padL, yToPx(thrS)); ctx.lineTo(W - padR, yToPx(thrS)); ctx.stroke();
      ctx.fillStyle = 'rgba(255,55,95,0.6)';
      ctx.fillText('rs·peak = ' + thrS.toFixed(2), W - padR - 110, yToPx(thrS) - 3);

      ctx.strokeStyle = 'rgba(52,199,89,0.3)';
      ctx.beginPath(); ctx.moveTo(padL, yToPx(thrD)); ctx.lineTo(W - padR, yToPx(thrD)); ctx.stroke();
      ctx.fillStyle = 'rgba(52,199,89,0.6)';
      ctx.fillText('rd·peak = ' + thrD.toFixed(2), W - padR - 110, yToPx(thrD) - 3);
      ctx.setLineDash([]);
    }
    if (p.method === 'variable' && est && isFinite(est.rs)) {
      const thrS = est.rs * env.peak;
      const thrD = est.rd * env.peak;
      ctx.strokeStyle = 'rgba(255,55,95,0.3)';
      ctx.setLineDash([3, 4]);
      ctx.beginPath(); ctx.moveTo(padL, yToPx(thrS)); ctx.lineTo(W - padR, yToPx(thrS)); ctx.stroke();
      ctx.fillStyle = 'rgba(255,55,95,0.6)';
      ctx.fillText('rs(adapt) = ' + est.rs.toFixed(2), W - padR - 130, yToPx(thrS) - 3);
      ctx.strokeStyle = 'rgba(52,199,89,0.3)';
      ctx.beginPath(); ctx.moveTo(padL, yToPx(thrD)); ctx.lineTo(W - padR, yToPx(thrD)); ctx.stroke();
      ctx.fillStyle = 'rgba(52,199,89,0.6)';
      ctx.fillText('rd(adapt) = ' + est.rd.toFixed(2), W - padR - 130, yToPx(thrD) - 3);
      ctx.setLineDash([]);
    }

    // envelope curve
    ctx.strokeStyle = '#0a84ff';
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    for (let i = 0; i < env.env.length; i++) {
      const x = xToPx(env.xs[i]);
      const y = yToPx(env.env[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    // soft fill under envelope
    ctx.lineTo(xToPx(env.xs[env.xs.length - 1]), yToPx(0));
    ctx.lineTo(xToPx(env.xs[0]), yToPx(0));
    ctx.closePath();
    ctx.fillStyle = 'rgba(10,132,255,0.10)';
    ctx.fill();

    // peak marker
    {
      const px = xToPx(env.peakP);
      const py = yToPx(env.peak);
      ctx.fillStyle = '#0a84ff';
      ctx.beginPath(); ctx.arc(px, py, 3.5, 0, Math.PI * 2); ctx.fill();
    }

    // ----- vertical markers -----------------------------------------------
    function vline(Pc, color, label, dashed, labelTop) {
      if (!isFinite(Pc) || Pc < ENV_PMIN || Pc > ENV_PMAX) return;
      const xx = xToPx(Pc);
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.4;
      if (dashed) ctx.setLineDash([5, 4]); else ctx.setLineDash([]);
      ctx.beginPath();
      ctx.moveTo(xx, padT);
      ctx.lineTo(xx, H - padB);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = color;
      ctx.font = '11px JetBrains Mono, Consolas, monospace';
      const ty = labelTop ? padT + 11 : H - padB - 6;
      ctx.fillText(label, xx + 3, ty);
    }

    // True MAP = envelope peak = zero transmural pressure = (SBP+DBP)/2 here.
    const trueMap = (p.sbp + p.dbp) / 2;
    // truth lines (solid, vivid)
    vline(p.sbp,   '#ff375f', 'SBP ' + p.sbp,           false, true);
    vline(trueMap, '#0a84ff', 'MAP ' + trueMap.toFixed(0), false, true);
    vline(p.dbp,   '#34c759', 'DBP ' + p.dbp,           false, true);

    // estimates (dashed, faded)
    if (est) {
      if (isFinite(est.sbp)) vline(est.sbp, 'rgba(255,55,95,0.55)',  '~SBP ' + est.sbp.toFixed(0), true, false);
      if (isFinite(est.map)) vline(est.map, 'rgba(10,132,255,0.55)', '~MAP ' + est.map.toFixed(0), true, false);
      if (isFinite(est.dbp)) vline(est.dbp, 'rgba(52,199,89,0.55)',  '~DBP ' + est.dbp.toFixed(0), true, false);
    }

    // panel title
    ctx.fillStyle = _cc.label();
    ctx.font = '11px JetBrains Mono, Consolas, monospace';
    ctx.fillText('extracted oscillogram envelope (sigmoid model, c = ' + p.cstiff + ')', padL, padT - 4);
  }

  function redraw() {
    drawCuffTrace();
    drawEnvelope();
  }

  // ----- play / scrub through deflation -----------------------------------
  let playRAF = 0;
  function setPlayLabel() {
    if (!inputs.play) return;
    inputs.play.textContent = playing ? 'Pause deflation' : 'Play deflation';
    inputs.play.setAttribute('aria-pressed', playing ? 'true' : 'false');
  }

  function playStep(now) {
    if (!playing) { playRAF = 0; return; }
    if (!cachedSim) { playRAF = requestAnimationFrame(playStep); return; }
    if (lastPlayT === 0) lastPlayT = now;
    const dt = (now - lastPlayT) / 1000;
    lastPlayT = now;
    // advance cursor by real-time seconds * FS_CUFF, scaled so the
    // playback duration matches the simulated deflation duration.
    const advance = Math.max(1, Math.round(dt * FS_CUFF));
    cursorIdx += advance;
    if (cursorIdx >= cachedSim.N) {
      cursorIdx = cachedSim.N - 1;
      playing = false;
      setPlayLabel();
    }
    drawCuffTrace();
    if (playing) playRAF = requestAnimationFrame(playStep);
  }

  function togglePlay() {
    if (!cachedSim) return;
    playing = !playing;
    if (playing && cursorIdx >= cachedSim.N - 1) cursorIdx = 0;
    lastPlayT = 0;
    setPlayLabel();
    if (playing && !playRAF) playRAF = requestAnimationFrame(playStep);
  }

  // ----- click on cuff canvas to scrub the cursor manually ----------------
  cuffEl.addEventListener('click', (ev) => {
    if (!cachedSim) return;
    const r = cuffEl.getBoundingClientRect();
    const x = ev.clientX - r.left;
    const padL = 38, padR = 12;
    const plotW = (cuffCV.W) - padL - padR;
    if (x < padL || x > cuffCV.W - padR) return;
    const frac = (x - padL) / plotW;
    cursorIdx = Math.max(0, Math.min(cachedSim.N - 1, Math.floor(frac * cachedSim.N)));
    if (playing) { playing = false; setPlayLabel(); }
    drawCuffTrace();
  });

  // ----- debounced recompute on slider change ------------------------------
  let recomputeTimer = 0;
  function scheduleRecompute() {
    if (recomputeTimer) clearTimeout(recomputeTimer);
    recomputeTimer = setTimeout(() => {
      recomputeTimer = 0;
      recompute();
    }, 50);
  }

  function bindInputs() {
    Object.values(inputs).forEach(el => {
      if (!el) return;
      if (el === inputs.play) {
        el.addEventListener('click', togglePlay);
      } else if (el.tagName === 'SELECT') {
        el.addEventListener('change', scheduleRecompute);
      } else {
        el.addEventListener('input', scheduleRecompute);
      }
    });
  }
  bindInputs();

  // ----- IntersectionObserver: skip resize/redraw work when offscreen -----
  let visible = true;
  if ('IntersectionObserver' in window) {
    const io = new IntersectionObserver((entries) => {
      entries.forEach(e => {
        if (e.target === cuffEl || e.target === envEl) {
          visible = e.isIntersecting;
          if (visible) redraw();
        }
      });
    }, { threshold: 0.01 });
    io.observe(cuffEl);
    io.observe(envEl);
  }

  // ----- handle responsive resize ------------------------------------------
  let resizeT = 0;
  window.addEventListener('resize', () => {
    if (resizeT) clearTimeout(resizeT);
    resizeT = setTimeout(() => {
      resizeT = 0;
      cuffCV.resize();
      envCV.resize();
      if (visible) redraw();
    }, 80);
  });

  // Pause autoplay when tab is hidden; resume when visible again
  document.addEventListener('visibilitychange', () => {
    if (!document.hidden && visible) redraw();
  });

  // Repaint on theme toggle so the canvas grid/axis/label colors follow the
  // theme (they're read live per draw via _cc, but a static frame won't refresh
  // on its own). Matches journeyCanvas / ecgBackground / bpWidget in main.js.
  document.addEventListener('themechange', () => { if (visible) redraw(); });

  // ----- first compute -----------------------------------------------------
  setPlayLabel();
  recompute();

})();
