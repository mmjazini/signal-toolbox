/* ==========================================================================
   R-peak detector — Pan–Tompkins on a synthetic ECG with realistic artifacts.

   Pipeline (Pan & Tompkins, "A Real-Time QRS Detection Algorithm",
   IEEE Trans. Biomed. Eng., 1985):
     1. Bandpass 5–15 Hz   (cascade of two single-pole IIRs: HP @5, LP @15)
     2. 5-point derivative y[n] = (1/8)(2x[n] + x[n-1] − x[n-3] − 2x[n-4])
     3. Square (point-wise)
     4. Moving-window integration over 150 ms (~38 samples @250 Hz)
     5. Adaptive double-threshold (SPKI/NPKI), 200 ms refractory,
        T-wave discriminator (slope check inside 360 ms of last R),
        and search-back when no peak appears within 1.66·RR_avg.

   Signal model (synthetic ECG): McSharry et al. 2003 — each beat is the
   sum of 5 Gaussians (P, Q, R, S, T) parameterised by HR, with optional
   PVC injections, baseline wander, EMG, powerline, and motion bursts.

   No external libraries. ~Single IIFE. Bails silently if its host
   <section id="rpeak"> is not on the page.
   ========================================================================== */
(function rpeakTool() {

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
  const root = document.getElementById('ecg-tool-canvas');
  if (!root) return;

  /* --------------------------------------------------------------------- */
  /* DOM                                                                   */
  /* --------------------------------------------------------------------- */
  const $ = (id) => document.getElementById(id);

  const ui = {
    hr:        $('rpeak-hr'),
    hrLbl:     $('rpeak-hr-val'),
    rmssd:     $('rpeak-rmssd'),
    rmssdLbl:  $('rpeak-rmssd-val'),
    drift:     $('rpeak-drift'),
    driftLbl:  $('rpeak-drift-val'),
    emg:       $('rpeak-emg'),
    emgLbl:    $('rpeak-emg-val'),
    plMode:    $('rpeak-plmode'),     // <select> or radios — read .value
    plAmp:     $('rpeak-plamp'),
    plAmpLbl:  $('rpeak-plamp-val'),
    motion:    $('rpeak-motion'),
    motionLbl: $('rpeak-motion-val'),
    pvc:       $('rpeak-pvc'),         // checkbox
    reset:     $('rpeak-reset'),

    // readouts
    trueHr:    $('rpeak-out-truehr'),
    detHr:     $('rpeak-out-dethr'),
    sens:      $('rpeak-out-sens'),
    ppv:       $('rpeak-out-ppv'),
    meanRR:    $('rpeak-out-meanrr'),
    rmssdOut:  $('rpeak-out-rmssd'),
    pnn50:     $('rpeak-out-pnn50'),
    latency:   $('rpeak-out-latency'),
    snr:       $('rpeak-out-snr'),
  };

  function bindSlider(el, lbl, unit, fmt = (v) => v) {
    if (!el) return () => 0;
    const update = () => {
      if (lbl) lbl.textContent = `${fmt(parseFloat(el.value))} ${unit}`.trim();
    };
    el.addEventListener('input', update);
    update();
    return () => parseFloat(el.value);
  }

  const getHR     = bindSlider(ui.hr,     ui.hrLbl,     'bpm');
  const getRMSSD  = bindSlider(ui.rmssd,  ui.rmssdLbl,  'ms');
  const getDrift  = bindSlider(ui.drift,  ui.driftLbl,  '',  v => v.toFixed(2));
  const getEMG    = bindSlider(ui.emg,    ui.emgLbl,    '',  v => v.toFixed(2));
  const getPLAmp  = bindSlider(ui.plAmp,  ui.plAmpLbl,  '',  v => v.toFixed(2));
  const getMotion = bindSlider(ui.motion, ui.motionLbl, '',  v => v.toFixed(2));
  const getPLMode = () => (ui.plMode ? ui.plMode.value : 'off');     // 'off' | '50' | '60'
  const getPVC    = () => (ui.pvc ? !!ui.pvc.checked : false);

  /* --------------------------------------------------------------------- */
  /* Canvas helpers — DPR-aware                                            */
  /* --------------------------------------------------------------------- */
  function makeCanvas(id) {
    const c = document.getElementById(id);
    if (!c) return null;
    const ctx = c.getContext('2d');
    let W = 0, H = 0, DPR = Math.min(window.devicePixelRatio || 1, 2);
    function resize() {
      // Honour the canvas' CSS-set size so tool-specific heights apply.
      const r = c.getBoundingClientRect();
      W = Math.max(120, Math.floor(r.width));
      H = Math.max(80,  Math.floor(r.height));
      c.width  = Math.round(W * DPR);
      c.height = Math.round(H * DPR);
      ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
    }
    resize();
    window.addEventListener('resize', resize);
    return { ctx, canvas: c, get W(){ return W; }, get H(){ return H; }, resize };
  }

  const top = makeCanvas('ecg-tool-canvas');
  const bot = makeCanvas('ecg-tool-features');
  if (!top) return;

  /* --------------------------------------------------------------------- */
  /* Constants                                                             */
  /* --------------------------------------------------------------------- */
  const Fs = 250;                              // sampling rate (Hz)
  const WINDOW_SEC = 10;
  const WIN = WINDOW_SEC * Fs;                 // 2500 samples
  const REFRACT_S  = Math.round(0.20 * Fs);    // 200 ms refractory
  const TWAVE_S    = Math.round(0.36 * Fs);    // 360 ms T-wave window
  const MWIN       = Math.round(0.150 * Fs);   // 150 ms integrator (≈38)

  /* --------------------------------------------------------------------- */
  /* PRNG — deterministic when prefers-reduced-motion is set               */
  /* --------------------------------------------------------------------- */
  const reduceMotion = window.matchMedia &&
                       window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  let rngState = 0x12345678 >>> 0;
  function rand() {
    if (reduceMotion) {
      // mulberry32
      rngState = (rngState + 0x6D2B79F5) >>> 0;
      let t = rngState;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    }
    return Math.random();
  }
  function gauss() {
    // Box–Muller
    let u = 0, v = 0;
    while (u === 0) u = rand();
    while (v === 0) v = rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  /* --------------------------------------------------------------------- */
  /* Synthetic ECG generator — McSharry-style 5-Gaussian beat              */
  /* --------------------------------------------------------------------- */
  // Each (theta_deg, amplitude, width_rad) for [P, Q, R, S, T]
  const BEAT_NORMAL = [
    { th: -70, a:  0.12, b: 0.25 },   // P
    { th: -15, a: -0.20, b: 0.10 },   // Q
    { th:   0, a:  1.20, b: 0.10 },   // R
    { th:  15, a: -0.30, b: 0.10 },   // S
    { th:  90, a:  0.30, b: 0.40 },   // T
  ];
  // PVC: wide, inverted, no P, large monomorphic deflection.
  const BEAT_PVC = [
    { th:   0, a: -1.40, b: 0.22 },
    { th:  30, a:  0.50, b: 0.30 },
  ];

  // Generator state — schedule beats in advance so detection latency can be
  // measured against the *true* R-time of each emitted beat.
  let beats = [];        // array of { tR (s, absolute), pvc:bool }
  let nextBeatT = 0;     // next scheduled R-time (s)
  let beatCounter = 0;
  let absT = 0;          // absolute time of newest sample (s)

  function scheduleBeats(uptoT, hrBpm, rmssdMs, pvcOn) {
    const meanRR = 60 / hrBpm;
    const sigma = (rmssdMs * 1e-3) / Math.SQRT2;   // RMSSD ≈ √2·σ_RR
    while (nextBeatT <= uptoT + 2.0) {              // schedule 2 s ahead
      const isPVC = pvcOn && (beatCounter % 10 === 9);
      beats.push({ tR: nextBeatT, pvc: isPVC });
      // Next interval: PVCs are early (compensatory pause follows)
      let rr;
      if (isPVC) {
        rr = meanRR * 0.65 + sigma * gauss();
      } else if (beats.length >= 2 && beats[beats.length - 2].pvc) {
        rr = meanRR * 1.35 + sigma * gauss();      // compensatory pause
      } else {
        rr = meanRR + sigma * gauss();
      }
      rr = Math.max(0.25, Math.min(2.5, rr));      // sanity clamp
      nextBeatT += rr;
      beatCounter++;
    }
    // garbage-collect beats older than the visible window + slack
    const cutoff = uptoT - (WINDOW_SEC + 1.0);
    while (beats.length && beats[0].tR < cutoff) beats.shift();
  }

  // sum-of-Gaussians around each scheduled beat
  function ecgSampleAt(t, hrBpm) {
    const meanRR = 60 / hrBpm;
    let v = 0;
    // only beats within ±1 s contribute meaningfully
    for (let i = 0; i < beats.length; i++) {
      const b = beats[i];
      const dt = t - b.tR;
      if (dt < -0.5 || dt > 0.6) continue;
      // map dt → angle relative to this beat's RR
      const th = (dt / meanRR) * 2 * Math.PI;       // radians
      const model = b.pvc ? BEAT_PVC : BEAT_NORMAL;
      for (let k = 0; k < model.length; k++) {
        const m = model[k];
        const dth = th - (m.th * Math.PI / 180);
        // wrap to (-pi, pi)
        const wd = ((dth + Math.PI) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) - Math.PI;
        v += m.a * Math.exp(- (wd * wd) / (2 * m.b * m.b));
      }
    }
    return v;
  }

  /* --------------------------------------------------------------------- */
  /* Artifact generators                                                   */
  /* --------------------------------------------------------------------- */
  // Motion-burst manager: bursts last 0.4–1.2 s, shaped like a Hann window
  let motionBursts = [];   // { t0, t1, amp, sign }
  let motionLastCheckT = 0;

  function maybeStartMotion(t, prob) {
    // probability per second → probability per dt
    const dtCheck = t - motionLastCheckT;
    if (dtCheck < 1 / 60) return;
    motionLastCheckT = t;
    const pPerCheck = 1 - Math.pow(1 - prob, dtCheck);
    if (rand() < pPerCheck) {
      const dur = 0.4 + 0.8 * rand();
      const amp = 0.6 + 1.4 * rand();
      const sgn = rand() < 0.5 ? -1 : 1;
      motionBursts.push({ t0: t, t1: t + dur, amp, sign: sgn });
    }
    // prune
    while (motionBursts.length && motionBursts[0].t1 < t - 0.5) motionBursts.shift();
  }
  function motionAt(t, severity) {
    if (severity <= 0 || motionBursts.length === 0) return 0;
    let v = 0;
    for (let i = 0; i < motionBursts.length; i++) {
      const b = motionBursts[i];
      if (t < b.t0 || t > b.t1) continue;
      const u = (t - b.t0) / (b.t1 - b.t0);            // 0..1
      const env = 0.5 - 0.5 * Math.cos(2 * Math.PI * u); // Hann
      v += b.sign * b.amp * env * severity;
    }
    return v;
  }

  /* --------------------------------------------------------------------- */
  /* Rolling buffers — raw ECG, integrator output, threshold trace         */
  /* --------------------------------------------------------------------- */
  const ecgBuf  = new Float64Array(WIN);   // visible raw signal
  const intBuf  = new Float64Array(WIN);   // integrator output (visible)
  const t1Buf   = new Float64Array(WIN);   // THRESHOLD_I1 trace
  const t2Buf   = new Float64Array(WIN);   // THRESHOLD_I2 trace
  // small history needed for filters/derivative (length up to ~12)
  const HBP_LEN = 6, DRV_LEN = 6;
  let bp_x = new Float64Array(HBP_LEN);    // bandpass intermediate
  let bp_y = new Float64Array(HBP_LEN);
  let drv_x = new Float64Array(DRV_LEN);   // derivative input history
  let writeIdx = 0;                         // ring write pointer (not used for raw)
  let firstSampleT = -WINDOW_SEC;           // absolute time at index 0 of buffer

  // mark for ground-truth peak times that intersect the visible window
  function visibleBeats() {
    const out = [];
    const t0 = firstSampleT;
    const t1 = firstSampleT + WINDOW_SEC;
    for (let i = 0; i < beats.length; i++) {
      if (beats[i].tR >= t0 && beats[i].tR <= t1) out.push(beats[i]);
    }
    return out;
  }

  /* --------------------------------------------------------------------- */
  /* Pan–Tompkins detector — runs over the current 10 s window each frame  */
  /* Detector state is persisted between frames (SPKI, NPKI, RR_avg).      */
  /* --------------------------------------------------------------------- */
  // IIR coefficients — single-pole HP @5 Hz and LP @15 Hz, used in cascade
  function makePoleFilters() {
    // bilinear-mapped 1-pole approximations (good enough for bandpass shape)
    const fc_lp = 15, fc_hp = 5;
    const wlp = Math.tan(Math.PI * fc_lp / Fs);
    const alp = wlp / (1 + wlp);
    const whp = Math.tan(Math.PI * fc_hp / Fs);
    const ahp = 1 / (1 + whp);
    return { alp, ahp };
  }
  const POLE = makePoleFilters();

  // Detector scratch state. NOTE: the visible window is fully re-filtered from
  // zero state every frame (see applyBandpass "re-prime per window"), so the
  // detector must likewise run as a fresh pass over the whole window each call
  // — it is reset at the top of runPanTompkins(). Persisting + decaying these
  // across frames while re-walking the same window re-counted every beat dozens
  // of times, polluting rrHist/rrAvg. Kept as an object only to avoid churn.
  const det = {
    SPKI: 0, NPKI: 0,
    rrHist: [],                  // last 8 RR (samples), within this window
    rrAvg: Math.round(0.8 * Fs), // initial guess
    lastSlope: 0,                // for T-wave discriminator
  };

  function applyBandpass(out, x, n) {
    // HP first, then LP — single-pole each. Persist state across runs would
    // require global filter memory; instead we re-prime per window with a
    // short warm-up (the first ~50 samples are not used for thresholding).
    let lp = 0;     // LP state
    let hp = 0, hp_xprev = 0;   // HP (first-order high-pass state)
    const alp = POLE.alp, ahp = POLE.ahp;
    for (let i = 0; i < n; i++) {
      // High-pass: y = ahp*(y_prev + x − x_prev)
      hp = ahp * (hp + x[i] - hp_xprev);
      hp_xprev = x[i];
      // Low-pass: y += alp*(x − y)
      lp = lp + alp * (hp - lp);
      out[i] = lp;
    }
  }

  function deriv5(out, x, n) {
    // y[n] = (1/8)(2x[n] + x[n-1] − x[n-3] − 2x[n-4])
    for (let i = 0; i < n; i++) {
      const x0 = x[i];
      const x1 = i >= 1 ? x[i - 1] : x[0];
      const x3 = i >= 3 ? x[i - 3] : x[0];
      const x4 = i >= 4 ? x[i - 4] : x[0];
      out[i] = (2 * x0 + x1 - x3 - 2 * x4) / 8;
    }
  }

  function squareInPlace(x, n) {
    for (let i = 0; i < n; i++) x[i] = x[i] * x[i];
  }

  function movingIntegrate(out, x, n, w) {
    // boxcar moving average length w
    let sum = 0;
    const inv = 1 / w;
    for (let i = 0; i < n; i++) {
      sum += x[i];
      if (i >= w) sum -= x[i - w];
      out[i] = sum * inv;
    }
  }

  // Detected R indices for the current window (sample indices into ecgBuf)
  let detectedR = [];     // { idx, isPVC?, latencyMs }
  // For sensitivity/PPV we keep a global record over the visible window
  // matched against ground-truth times.

  function runPanTompkins(integ, n, windowStartT) {
    // integ: Float64Array length n (≈WIN)
    // Fresh pass every call — the whole window is re-filtered from zero each
    // frame, so the adaptive estimators are re-primed here rather than carried
    // (and decayed) across frames, which would re-count every beat repeatedly.
    detectedR.length = 0;
    det.rrHist.length = 0;
    det.rrAvg = Math.round(0.8 * Fs);
    det.lastSlope = 0;

    // 1) find local maxima (peaks) in integ — candidates for R
    //    use a simple "greater than both neighbours by ≥small" rule.
    const candIdx = [];
    for (let i = 1; i < n - 1; i++) {
      if (integ[i] > integ[i - 1] && integ[i] >= integ[i + 1] && integ[i] > 1e-9) {
        candIdx.push(i);
      }
    }

    // Initialise SPKI/NPKI from the first 2 s of this window.
    {
      let mx = 0;
      const init_n = Math.min(n, 2 * Fs);
      for (let i = 0; i < init_n; i++) {
        if (integ[i] > mx) mx = integ[i];
      }
      det.SPKI = mx * 0.25;
      det.NPKI = mx * 0.05;
    }

    let lastAcceptIdx = -Infinity;     // sample idx within window
    let lastAcceptVal = 0;             // integrator height
    let lastRR = det.rrAvg;
    const rrAccept = [];

    // Helper — accept candidate at i (sample idx) using value v
    function accept(i, v, isSearchback) {
      // T-wave discriminator: if too soon after last R but past refractory,
      // compare slope of integrator around this peak to the prior beat's slope.
      const since = i - lastAcceptIdx;
      if (since > REFRACT_S && since < TWAVE_S && lastAcceptIdx > -Infinity) {
        // current slope ≈ integ[i] − integ[i − ~5]
        const k = Math.max(0, i - 5);
        const slope = integ[i] - integ[k];
        if (slope < 0.5 * det.lastSlope) {
          // probable T-wave: still update NPKI, do not accept
          det.NPKI = 0.125 * v + 0.875 * det.NPKI;
          return false;
        }
      }
      // Accept
      detectedR.push({ idx: i });
      // RR bookkeeping (skip the very first accepted)
      if (lastAcceptIdx > -Infinity) {
        const rr = i - lastAcceptIdx;
        if (rr > REFRACT_S) {
          rrAccept.push(rr);
          det.rrHist.push(rr);
          if (det.rrHist.length > 8) det.rrHist.shift();
          // update RR_avg
          let s = 0;
          for (let k2 = 0; k2 < det.rrHist.length; k2++) s += det.rrHist[k2];
          det.rrAvg = Math.round(s / det.rrHist.length);
          lastRR = det.rrAvg;
        }
      }
      // Slope memory for next T-wave check
      const k = Math.max(0, i - 5);
      det.lastSlope = integ[i] - integ[k];
      // Update SPKI (search-back uses softer update)
      if (isSearchback) det.SPKI = 0.25 * v + 0.75 * det.SPKI;
      else              det.SPKI = 0.125 * v + 0.875 * det.SPKI;
      lastAcceptIdx = i;
      lastAcceptVal = v;
      return true;
    }

    // Walk candidates with refractory + adaptive threshold
    let candPtr = 0;
    let lastWalk = -REFRACT_S;
    while (candPtr < candIdx.length) {
      const i = candIdx[candPtr];
      const v = integ[i];
      const TH1 = det.NPKI + 0.25 * (det.SPKI - det.NPKI);

      if ((i - lastWalk) < REFRACT_S) { candPtr++; continue; }

      // Search-back: if it's been > 1.66·RR_avg since last accept, drop to TH2
      const sinceLast = (lastAcceptIdx > -Infinity)
        ? (i - lastAcceptIdx)
        : (i - 0);
      let TH = TH1;
      let sb = false;
      if (lastAcceptIdx > -Infinity && sinceLast > Math.round(1.66 * det.rrAvg)) {
        TH = 0.5 * TH1;
        sb = true;
      }

      if (v >= TH) {
        const ok = accept(i, v, sb);
        if (ok) lastWalk = i;
        else    lastWalk = i;        // still respect refractory
      } else {
        // Not a peak — update noise estimate
        det.NPKI = 0.125 * v + 0.875 * det.NPKI;
      }
      candPtr++;
    }
  }

  /* --------------------------------------------------------------------- */
  /* HRV metrics over detected peaks (visible window)                      */
  /* --------------------------------------------------------------------- */
  function hrvMetrics(rIdxArr) {
    if (rIdxArr.length < 2) return { meanRR: 0, rmssd: 0, pnn50: 0, hr: 0 };
    const rr = [];
    for (let i = 1; i < rIdxArr.length; i++) {
      rr.push((rIdxArr[i] - rIdxArr[i - 1]) / Fs * 1000);  // ms
    }
    let s = 0; for (let i = 0; i < rr.length; i++) s += rr[i];
    const mean = s / rr.length;
    let sqd = 0, n50 = 0;
    for (let i = 1; i < rr.length; i++) {
      const d = rr[i] - rr[i - 1];
      sqd += d * d;
      if (Math.abs(d) > 50) n50++;
    }
    const rmssd = rr.length > 1 ? Math.sqrt(sqd / (rr.length - 1)) : 0;
    const pnn50 = rr.length > 1 ? (100 * n50 / (rr.length - 1)) : 0;
    // median RR for HR (robust against outliers)
    const sorted = rr.slice().sort((a, b) => a - b);
    const med = sorted[Math.floor(sorted.length / 2)];
    const hr = med > 0 ? 60000 / med : 0;
    return { meanRR: mean, rmssd, pnn50, hr };
  }

  /* --------------------------------------------------------------------- */
  /* Match detected → ground-truth (sensitivity / PPV / latency)           */
  /* Tolerance: ±100 ms.                                                   */
  /* --------------------------------------------------------------------- */
  function matchPeaks(detectedTimes, truthTimes, tolMs = 100) {
    const tol = tolMs / 1000;
    const matchedTruth = new Set();
    const matchedDet   = new Set();
    const latencies = [];
    // greedy nearest-match
    for (let i = 0; i < detectedTimes.length; i++) {
      let bestJ = -1, bestD = Infinity;
      for (let j = 0; j < truthTimes.length; j++) {
        if (matchedTruth.has(j)) continue;
        const d = detectedTimes[i] - truthTimes[j];
        if (Math.abs(d) < bestD) { bestD = Math.abs(d); bestJ = j; }
      }
      if (bestJ >= 0 && bestD <= tol) {
        matchedTruth.add(bestJ);
        matchedDet.add(i);
        // Positive = detector fires after truth R; negative = detector fires before.
    // A small negative value (≈ −10–20 ms) is expected here: the group-delay
    // correction d.idx − MWIN/2 maps to ≈ the true R index, but the synthetic
    // generator schedules tR at the Gaussian peak center while the detector's
    // threshold crossing can fire a few samples early on a fast upstroke.
    latencies.push((detectedTimes[i] - truthTimes[bestJ]) * 1000);
      }
    }
    const TP = matchedDet.size;
    const FN = truthTimes.length - matchedTruth.size;
    const FP = detectedTimes.length - matchedDet.size;
    const sens = (TP + FN) > 0 ? TP / (TP + FN) : 0;
    const ppv  = (TP + FP) > 0 ? TP / (TP + FP) : 0;
    let medLat = 0;
    if (latencies.length) {
      latencies.sort((a, b) => a - b);
      medLat = latencies[Math.floor(latencies.length / 2)];
    }
    return { sens, ppv, TP, FP, FN, medLat };
  }

  /* --------------------------------------------------------------------- */
  /* SNR estimate — RMS(beats around R ±60 ms) vs RMS(rest)                */
  /* --------------------------------------------------------------------- */
  function snrEstimate(ecg, n, peakIdxs) {
    if (peakIdxs.length < 2) return { mean: 0, sd: 0 };
    const inBeat = new Uint8Array(n);
    const halfW = Math.round(0.06 * Fs);
    for (let k = 0; k < peakIdxs.length; k++) {
      const i = peakIdxs[k];
      for (let j = Math.max(0, i - halfW); j < Math.min(n, i + halfW); j++) inBeat[j] = 1;
    }
    // sliding 1-s windows for variability
    const winN = Fs;
    const snrs = [];
    for (let s = 0; s + winN < n; s += Fs) {
      let sigSq = 0, sigCnt = 0, noiseSq = 0, noiseCnt = 0;
      for (let i = s; i < s + winN; i++) {
        const v = ecg[i];
        if (inBeat[i]) { sigSq += v * v; sigCnt++; }
        else           { noiseSq += v * v; noiseCnt++; }
      }
      if (sigCnt > 0 && noiseCnt > 0) {
        const sig = Math.sqrt(sigSq / sigCnt);
        const noi = Math.sqrt(noiseSq / noiseCnt);
        if (noi > 1e-9) snrs.push(20 * Math.log10(sig / noi));
      }
    }
    if (snrs.length === 0) return { mean: 0, sd: 0 };
    let m = 0; for (let i = 0; i < snrs.length; i++) m += snrs[i]; m /= snrs.length;
    let v = 0; for (let i = 0; i < snrs.length; i++) v += (snrs[i] - m) * (snrs[i] - m);
    return { mean: m, sd: Math.sqrt(v / snrs.length) };
  }

  /* --------------------------------------------------------------------- */
  /* Drawing                                                               */
  /* --------------------------------------------------------------------- */
  function drawGrid(ctx, W, H, secondsVisible) {
    ctx.save();
    // major: 1-s ticks
    ctx.strokeStyle = 'rgba(255,255,255,0.06)';
    ctx.lineWidth = 1;
    for (let s = 0; s <= secondsVisible; s++) {
      const x = (s / secondsVisible) * W;
      ctx.beginPath(); ctx.moveTo(x + 0.5, 0); ctx.lineTo(x + 0.5, H); ctx.stroke();
    }
    // minor horizontal
    ctx.strokeStyle = _cc.grid();
    for (let k = 1; k < 6; k++) {
      const y = (k / 6) * H;
      ctx.beginPath(); ctx.moveTo(0, y + 0.5); ctx.lineTo(W, y + 0.5); ctx.stroke();
    }
    // 1-s tick labels along bottom
    ctx.fillStyle = _cc.labelDim();
    ctx.font = '10px ui-sans-serif, system-ui, sans-serif';
    ctx.textBaseline = 'bottom';
    for (let s = 0; s <= secondsVisible; s++) {
      const x = (s / secondsVisible) * W;
      ctx.fillText(`${s}s`, x + 3, H - 2);
    }
    ctx.restore();
  }

  function drawECGPanel(ecg, n, truthBeats, detIdxs, motionBands) {
    const ctx = top.ctx, W = top.W, H = top.H;
    ctx.clearRect(0, 0, W, H);
    drawGrid(ctx, W, H, WINDOW_SEC);

    // motion-burst shading (each band: t0,t1 in absolute seconds)
    if (motionBands.length) {
      ctx.fillStyle = 'rgba(255,55,95,0.10)';
      const t0win = firstSampleT;
      for (let k = 0; k < motionBands.length; k++) {
        const m = motionBands[k];
        const x0 = ((m.t0 - t0win) / WINDOW_SEC) * W;
        const x1 = ((m.t1 - t0win) / WINDOW_SEC) * W;
        const X0 = Math.max(0, x0), X1 = Math.min(W, x1);
        if (X1 > X0) ctx.fillRect(X0, 0, X1 - X0, H);
      }
    }

    // Auto-scale ECG to visible amplitude
    let mn = Infinity, mx = -Infinity;
    for (let i = 0; i < n; i++) {
      if (ecg[i] < mn) mn = ecg[i];
      if (ecg[i] > mx) mx = ecg[i];
    }
    const amp = Math.max(0.6, mx - mn);
    const yMid = H * 0.55;
    const yScale = (H * 0.78) / amp;

    // Trace
    ctx.strokeStyle = '#34c759';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const x = (i / (n - 1)) * W;
      const y = yMid - (ecg[i] - (mn + mx) / 2) * yScale;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // Ground-truth peaks: open circles
    ctx.strokeStyle = _cc.label();
    ctx.lineWidth = 1;
    for (let k = 0; k < truthBeats.length; k++) {
      const t = truthBeats[k].tR;
      const xf = (t - firstSampleT) / WINDOW_SEC;
      if (xf < 0 || xf > 1) continue;
      const x = xf * W;
      // y at scheduled R: just put a dot at the top arc
      const y = yMid - (1.0) * yScale * 0.7;
      ctx.beginPath(); ctx.arc(x, y, 5, 0, Math.PI * 2); ctx.stroke();
      if (truthBeats[k].pvc) {
        ctx.fillStyle = 'rgba(255,255,255,0.7)';
        ctx.font = '10px ui-sans-serif, system-ui';
        ctx.fillText('PVC', x + 6, y - 6);
      }
    }

    // Detected peaks: filled red dots snapped to the local R-peak maximum.
    // The integrator fires MWIN/2 samples late; we correct the index then search
    // a ±40 ms window for the signal maximum so the dot sits on the R-peak tip.
    const snapHalf = Math.round(0.040 * Fs);  // ±40 ms search window
    ctx.fillStyle = '#ff375f';
    for (let k = 0; k < detIdxs.length; k++) {
      const i = detIdxs[k];
      const iCorr = Math.max(0, Math.min(n - 1, i - Math.floor(MWIN / 2)));
      const lo = Math.max(0, iCorr - snapHalf);
      const hi = Math.min(n - 1, iCorr + snapHalf);
      let peakIdx = iCorr;
      for (let j = lo; j <= hi; j++) {
        if (ecg[j] > ecg[peakIdx]) peakIdx = j;
      }
      const x = (peakIdx / (n - 1)) * W;
      const y = yMid - (ecg[peakIdx] - (mn + mx) / 2) * yScale;
      ctx.beginPath(); ctx.arc(x, y, 3.5, 0, Math.PI * 2); ctx.fill();
    }

    // Title strip
    ctx.fillStyle = _cc.label();
    ctx.font = '11px ui-sans-serif, system-ui';
    ctx.textBaseline = 'top';
    ctx.fillText('ECG (synthetic) — green: signal · open circles: truth · red: detected', 8, 6);
  }

  function drawFeaturePanel(integ, n, t1Trace, t2Trace, detIdxs) {
    if (!bot) return;
    const ctx = bot.ctx, W = bot.W, H = bot.H;
    ctx.clearRect(0, 0, W, H);
    drawGrid(ctx, W, H, WINDOW_SEC);

    // Auto-scale
    let mx = 0;
    for (let i = 0; i < n; i++) if (integ[i] > mx) mx = integ[i];
    for (let i = 0; i < n; i++) if (t1Trace[i] > mx) mx = t1Trace[i];
    if (mx < 1e-6) mx = 1;
    const yScale = (H * 0.8) / mx;
    const yBase = H - 8;

    // Integrator
    ctx.strokeStyle = 'rgba(10,132,255,0.7)';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const x = (i / (n - 1)) * W;
      const y = yBase - integ[i] * yScale;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // THRESHOLD_I1 (dashed)
    ctx.strokeStyle = '#ffd60a';
    ctx.setLineDash([6, 4]);
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const x = (i / (n - 1)) * W;
      const y = yBase - t1Trace[i] * yScale;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.setLineDash([]);

    // THRESHOLD_I2 (dotted, fainter)
    ctx.strokeStyle = 'rgba(255,214,10,0.45)';
    ctx.setLineDash([2, 4]);
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const x = (i / (n - 1)) * W;
      const y = yBase - t2Trace[i] * yScale;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.setLineDash([]);

    // Detected peak ticks at the bottom
    ctx.fillStyle = '#ff375f';
    for (let k = 0; k < detIdxs.length; k++) {
      const x = (detIdxs[k] / (n - 1)) * W;
      ctx.fillRect(x - 1, H - 6, 2, 6);
    }

    // Title
    ctx.fillStyle = _cc.label();
    ctx.font = '11px ui-sans-serif, system-ui';
    ctx.textBaseline = 'top';
    ctx.fillText('Pan–Tompkins integrator — dashed: TH1 (NPKI + 0.25·(SPKI−NPKI))   ·   dotted: TH2 (0.5·TH1)', 8, 6);
  }

  /* --------------------------------------------------------------------- */
  /* Main loop                                                             */
  /* --------------------------------------------------------------------- */
  // scratch arrays sized to WIN
  const scratchBP    = new Float64Array(WIN);
  const scratchDRV   = new Float64Array(WIN);
  const scratchSQ    = new Float64Array(WIN);
  const scratchINT   = new Float64Array(WIN);

  let running = true;
  let rafId = 0;

  // IntersectionObserver — pause when section is offscreen
  const section = document.getElementById('rpeak') || top.canvas;
  if (section && 'IntersectionObserver' in window) {
    const io = new IntersectionObserver((entries) => {
      for (const e of entries) running = e.isIntersecting;
    }, { threshold: 0.05 });
    io.observe(section);
  }

  function reset() {
    beats = [];
    nextBeatT = 0;
    beatCounter = 0;
    absT = 0;
    firstSampleT = -WINDOW_SEC;
    motionBursts = [];
    motionLastCheckT = 0;
    ecgBuf.fill(0);
    intBuf.fill(0);
    t1Buf.fill(0);
    t2Buf.fill(0);
    det.SPKI = 0; det.NPKI = 0; det.rrHist = []; det.rrAvg = Math.round(0.8 * Fs);
    det.lastSlope = 0;
  }
  if (ui.reset) ui.reset.addEventListener('click', reset);

  // Tick: append `samplesPerFrame` new samples
  const samplesPerFrame = Math.max(1, Math.round(Fs / 60));   // ≈4 @60 Hz
  let lastWallT = performance.now();

  function step() {
    rafId = requestAnimationFrame(step);
    if (!running) { lastWallT = performance.now(); return; }
    const now = performance.now();
    let dtMs = now - lastWallT; lastWallT = now;
    if (dtMs > 80) dtMs = 80;            // cap on long stalls

    // We push exactly samplesPerFrame samples to keep the X-axis smooth.
    const nNew = samplesPerFrame;

    // Read UI
    const hrBpm   = getHR();
    const rmssdMs = getRMSSD();
    const drift   = getDrift();
    const emg     = getEMG();
    const pl      = getPLAmp();
    const plMode  = getPLMode();
    const motionP = getMotion();
    const pvcOn   = getPVC();

    // Schedule beats forward of the new samples
    scheduleBeats(absT + nNew / Fs + 1.0, hrBpm, rmssdMs, pvcOn);

    // Generate samples
    for (let s = 0; s < nNew; s++) {
      const t = absT + s / Fs;
      // possibly start a motion burst
      maybeStartMotion(t, motionP);

      let v = ecgSampleAt(t, hrBpm);

      // Baseline wander
      if (drift > 0) {
        v += drift * (0.6 * Math.sin(2 * Math.PI * 0.15 * t) +
                      0.25 * Math.sin(2 * Math.PI * 0.4 * t));
      }
      // EMG (white Gaussian)
      if (emg > 0) v += emg * 0.25 * gauss();
      // Powerline
      if (plMode !== 'off' && pl > 0) {
        const f = (plMode === '50') ? 50 : 60;
        v += pl * 0.25 * Math.sin(2 * Math.PI * f * t);
      }
      // Motion
      v += motionAt(t, motionP);

      // Slide buffer
      ecgBuf.copyWithin(0, 1);
      ecgBuf[WIN - 1] = v;
    }

    absT += nNew / Fs;
    firstSampleT = absT - WINDOW_SEC;

    // ---- Detector pipeline (over current 10 s window) ------------------
    applyBandpass(scratchBP, ecgBuf, WIN);
    deriv5(scratchDRV, scratchBP, WIN);
    // square into scratchSQ
    for (let i = 0; i < WIN; i++) scratchSQ[i] = scratchDRV[i] * scratchDRV[i];
    movingIntegrate(scratchINT, scratchSQ, WIN, MWIN);

    // Run Pan–Tompkins decision logic
    runPanTompkins(scratchINT, WIN, firstSampleT);

    // Build threshold traces (TH1, TH2 are scalar at any moment but we draw
    // the *current* values across the visible axis as horizontal lines —
    // the clean, honest visualisation).
    const TH1 = det.NPKI + 0.25 * (det.SPKI - det.NPKI);
    const TH2 = 0.5 * TH1;
    for (let i = 0; i < WIN; i++) { t1Buf[i] = TH1; t2Buf[i] = TH2; intBuf[i] = scratchINT[i]; }

    // ---- Metrics --------------------------------------------------------
    const truth = visibleBeats();
    const truthTimes = truth.map(b => b.tR);
    // detected times (absolute) — undo integrator group delay (≈MWIN/2) to
    // line them up with R-times in the generator's clock.
    const detTimes = detectedR.map(d => firstSampleT + (d.idx - MWIN / 2) / Fs);

    const m = matchPeaks(detTimes, truthTimes, 100);
    const detIdxs = detectedR.map(d => d.idx);
    const hrv = hrvMetrics(detIdxs);
    const snr = snrEstimate(ecgBuf, WIN, detIdxs);

    // ---- Draw -----------------------------------------------------------
    const motionBands = motionBursts.filter(b => b.t1 >= firstSampleT && b.t0 <= absT);
    drawECGPanel(ecgBuf, WIN, truth, detIdxs, motionBands);
    drawFeaturePanel(intBuf, WIN, t1Buf, t2Buf, detIdxs);

    // ---- Readouts -------------------------------------------------------
    const fmt = (x, n = 1) => (isFinite(x) ? x.toFixed(n) : '—');
    if (ui.trueHr)   ui.trueHr.textContent   = `${fmt(hrBpm, 0)} bpm`;
    if (ui.detHr)    ui.detHr.textContent    = `${fmt(hrv.hr, 0)} bpm`;
    if (ui.sens)     ui.sens.textContent     = `${fmt(m.sens * 100, 1)}%`;
    if (ui.ppv)      ui.ppv.textContent      = `${fmt(m.ppv * 100, 1)}%`;
    if (ui.meanRR)   ui.meanRR.textContent   = `${fmt(hrv.meanRR, 0)} ms`;
    if (ui.rmssdOut) ui.rmssdOut.textContent = `${fmt(hrv.rmssd, 1)} ms`;
    if (ui.pnn50)    ui.pnn50.textContent    = `${fmt(hrv.pnn50, 1)}%`;
    if (ui.latency)  ui.latency.textContent  = `${fmt(m.medLat, 1)} ms`;
    if (ui.snr)      ui.snr.textContent      = `${fmt(snr.mean, 1)} ± ${fmt(snr.sd, 1)} dB`;
  }

  // Kick off
  reset();
  rafId = requestAnimationFrame(step);

  // -- Visibility hide → pause RAF entirely (saves cycles)
  document.addEventListener('visibilitychange', () => {
    if (document.hidden) {
      if (rafId) cancelAnimationFrame(rafId);
      rafId = 0;
    } else if (!rafId) {
      lastWallT = performance.now();
      rafId = requestAnimationFrame(step);
    }
  });
})();
