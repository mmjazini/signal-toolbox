/* ============================================================
 * fft.js — Interactive FFT / spectrum explorer
 *
 * Required DOM (all under <section id="fft">):
 *   Canvases:    #fft-time, #fft-freq
 *   Controls:    #fft-preset, #fft-fs (+ #fft-fs-val),
 *                #fft-nsamples, #fft-window,
 *                #fft-noise (+ #fft-noise-val)
 *   3sin:        #fft-f1..3, #fft-a1..3 (+ -val labels)
 *   Optional:    #fft-hr, #fft-saw-decay, #fft-chirp-f0/-f1,
 *                #fft-am-fc, #fft-am-fm (+ -val labels)
 *   Readouts:    #fft-readout-peak1..3, #fft-readout-thd
 *   Wrapper:     #fft-preset-wrap (data-preset attr)
 *
 * Bails silently if #fft-time is absent.
 * ============================================================ */
(function fftTool() {
    'use strict';

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

    const timeCanvas = document.getElementById('fft-time');
    if (!timeCanvas) return;
    const freqCanvas = document.getElementById('fft-freq');
    if (!freqCanvas) return;

    const timeCtx = timeCanvas.getContext('2d');
    const freqCtx = freqCanvas.getContext('2d');

    // --- DOM helpers ---------------------------------------------------
    const $ = (id) => document.getElementById(id);
    const setText = (id, v) => { const el = $(id); if (el) el.textContent = v; };
    const numVal = (id, fallback) => {
        const el = $(id);
        if (!el) return fallback;
        const n = parseFloat(el.value);
        return Number.isFinite(n) ? n : fallback;
    };
    const strVal = (id, fallback) => {
        const el = $(id);
        return el ? el.value : fallback;
    };

    // --- Twiddle cache (per N) ----------------------------------------
    const twiddleCache = new Map();
    function getTwiddles(N) {
        let t = twiddleCache.get(N);
        if (t) return t;
        // Pre-compute cos(-2πk/N), sin(-2πk/N) for k = 0..N/2-1.
        // Cooley-Tukey uses these as the "rotation factors" multiplying
        // the odd-indexed butterfly half.
        const half = N >> 1;
        const cos = new Float64Array(half);
        const sin = new Float64Array(half);
        for (let k = 0; k < half; k++) {
            const a = -2 * Math.PI * k / N;
            cos[k] = Math.cos(a);
            sin[k] = Math.sin(a);
        }
        t = { cos, sin };
        twiddleCache.set(N, t);
        return t;
    }

    // --- Bit reversal --------------------------------------------------
    function bitReverseCopy(real, re, im) {
        const N = real.length;
        const bits = Math.log2(N) | 0;
        for (let i = 0; i < N; i++) {
            // Reverse the bits of i in `bits` slots, then place x[i]
            // into the bit-reversed position. This is the standard
            // permutation that makes Cooley-Tukey in-place work.
            let j = 0, x = i;
            for (let b = 0; b < bits; b++) {
                j = (j << 1) | (x & 1);
                x >>= 1;
            }
            re[j] = real[i];
            im[j] = 0;
        }
    }

    // --- Iterative radix-2 in-place Cooley-Tukey FFT -------------------
    // Real input -> complex {re, im}. N must be a power of two.
    // Decomposes the DFT into log2(N) stages of butterfly pairs:
    //   X_even + W_N^k · X_odd   /   X_even - W_N^k · X_odd.
    function fft(real) {
        const N = real.length;
        const re = new Float64Array(N);
        const im = new Float64Array(N);
        bitReverseCopy(real, re, im);

        const { cos, sin } = getTwiddles(N);

        // size = current butterfly span (2,4,8,...,N).
        for (let size = 2; size <= N; size <<= 1) {
            const half = size >> 1;
            const step = N / size; // index stride into the global twiddle table
            for (let i = 0; i < N; i += size) {
                for (let k = 0; k < half; k++) {
                    const ti = k * step;
                    const c = cos[ti];
                    const s = sin[ti];
                    const aRe = re[i + k];
                    const aIm = im[i + k];
                    const bRe = re[i + k + half];
                    const bIm = im[i + k + half];
                    // t = W * b
                    const tRe = c * bRe - s * bIm;
                    const tIm = c * bIm + s * bRe;
                    re[i + k] = aRe + tRe;
                    im[i + k] = aIm + tIm;
                    re[i + k + half] = aRe - tRe;
                    im[i + k + half] = aIm - tIm;
                }
            }
        }
        return { re, im };
    }

    // --- Box-Muller Gaussian ------------------------------------------
    // Returns one N(0,1) sample per call. Two uniforms get one normal;
    // we discard the cosine partner for simplicity.
    function boxMuller() {
        let u1 = 0, u2 = 0;
        while (u1 === 0) u1 = Math.random();
        while (u2 === 0) u2 = Math.random();
        return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    }

    // --- Window functions ---------------------------------------------
    // Each of these is the standard textbook form on n = 0..N-1.
    // Flat-top is the SRS/Heinzel 5-term variant (low amplitude error).
    // Tukey is a cosine-tapered rectangle with α controlling the taper.
    function applyWindow(x, kind) {
        const N = x.length;
        const y = new Float64Array(N);
        if (kind === 'rect') {
            y.set(x);
            return y;
        }
        for (let n = 0; n < N; n++) {
            let w = 1;
            const r = n / (N - 1);
            if (kind === 'hann') {
                w = 0.5 * (1 - Math.cos(2 * Math.PI * r));
            } else if (kind === 'hamming') {
                w = 0.54 - 0.46 * Math.cos(2 * Math.PI * r);
            } else if (kind === 'blackman') {
                w = 0.42 - 0.5 * Math.cos(2 * Math.PI * r)
                        + 0.08 * Math.cos(4 * Math.PI * r);
            } else if (kind === 'flattop') {
                const a0 = 0.21557895, a1 = 0.41663158, a2 = 0.277263158;
                const a3 = 0.083578947, a4 = 0.006947368;
                w = a0
                    - a1 * Math.cos(2 * Math.PI * r)
                    + a2 * Math.cos(4 * Math.PI * r)
                    - a3 * Math.cos(6 * Math.PI * r)
                    + a4 * Math.cos(8 * Math.PI * r);
            } else if (kind === 'tukey') {
                const alpha = 0.5;
                if (r < alpha / 2) {
                    w = 0.5 * (1 + Math.cos(Math.PI * (2 * r / alpha - 1)));
                } else if (r > 1 - alpha / 2) {
                    w = 0.5 * (1 + Math.cos(Math.PI * (2 * r / alpha - 2 / alpha + 1)));
                } else {
                    w = 1;
                }
            }
            y[n] = x[n] * w;
        }
        return y;
    }

    // --- Parabolic interpolation around bin β -------------------------
    // Given dB magnitudes at k-1, k, k+1, fit a parabola through them
    // and return the sub-bin offset δ of its vertex from k. Standard
    // result: δ = 0.5·(α-γ)/(α-2β+γ).
    function parabolic(alpha, beta, gamma) {
        const denom = (alpha - 2 * beta + gamma);
        if (denom === 0) return 0;
        const d = 0.5 * (alpha - gamma) / denom;
        if (!Number.isFinite(d)) return 0;
        return Math.max(-1, Math.min(1, d));
    }

    // --- Peak finding --------------------------------------------------
    function findPeaks(specDb, freqs, k) {
        const peaks = [];
        const N = specDb.length;
        const floor = -80;
        // Local-maxima sweep over the positive-frequency half. We skip DC
        // (i=0) but scan down to bin 1 so a low-frequency fundamental (e.g. a
        // ~1 Hz ECG/PPG tone that lands in bin 1 at large N) is still found.
        // Neighbor lookups are clamped so the 5-point test works at the edges.
        const at = (j) => specDb[Math.max(0, Math.min(N - 1, j))];
        for (let i = 1; i < N - 1; i++) {
            const v = specDb[i];
            if (v < floor + 5) continue;
            if (v > at(i - 1) && v > at(i + 1)
                && v >= at(i - 2) && v >= at(i + 2)) {
                const off = parabolic(at(i - 1), v, at(i + 1));
                const interpFreq = freqs[i] + off * (freqs[1] - freqs[0]);
                peaks.push({ bin: i, freq: interpFreq, dB: v });
            }
        }
        peaks.sort((a, b) => b.dB - a.dB);
        return peaks.slice(0, k);
    }

    // --- THD ----------------------------------------------------------
    // RMS of fundamental + first 4 harmonics divided by RMS of all
    // remaining bins, expressed as a percentage. Only meaningful when
    // a single tonal component dominates the spectrum.
    function computeTHD(reIm, freqs, peaks) {
        if (peaks.length < 4) return null;
        // Confidence check: top peak must be 6 dB above the 4th-strongest peak.
        if (peaks[0].dB - peaks[3].dB < 6) return null;
        const f0 = peaks[0].freq;
        if (f0 <= 0) return null;

        const { re, im } = reIm;
        const N = re.length;
        const half = N >> 1;
        const df = freqs[1] - freqs[0];

        // Collect the bins belonging to the fundamental + 2..5 harmonics,
        // each over a small band to capture leakage spillover. (Energy is
        // summed below into fundEnergy / upperEnergy.)
        const harmonicBins = new Set();
        for (let h = 1; h <= 5; h++) {
            const fh = f0 * h;
            if (fh >= freqs[half - 1]) break;
            const center = Math.round(fh / df);
            const radius = 2;
            for (let b = center - radius; b <= center + radius; b++) {
                if (b < 1 || b >= half) continue;
                harmonicBins.add(b);
            }
        }

        const fundBins = new Set();
        const center0 = Math.round(f0 / df);
        for (let b = center0 - 2; b <= center0 + 2; b++) {
            if (b >= 1 && b < half) fundBins.add(b);
        }
        let fundEnergy = 0;
        for (const b of fundBins) {
            fundEnergy += re[b] * re[b] + im[b] * im[b];
        }
        if (fundEnergy <= 0) return null;
        let upperEnergy = 0;
        for (const b of harmonicBins) {
            if (!fundBins.has(b)) {
                upperEnergy += re[b] * re[b] + im[b] * im[b];
            }
        }
        // Classic definition: THD = sqrt(sum of harmonic powers) / fundamental amplitude.
        const thd = Math.sqrt(upperEnergy / fundEnergy) * 100;
        if (!Number.isFinite(thd)) return null;
        return thd;
    }

    // --- Signal builder -----------------------------------------------
    function buildSignal(params) {
        const { N, Fs, preset, noise } = params;
        const x = new Float64Array(N);
        const dt = 1 / Fs;
        const T = N * dt;

        if (preset === '3sin') {
            const f1 = params.f1, f2 = params.f2, f3 = params.f3;
            const a1 = params.a1, a2 = params.a2, a3 = params.a3;
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                x[n] = a1 * Math.sin(2 * Math.PI * f1 * t)
                     + a2 * Math.sin(2 * Math.PI * f2 * t)
                     + a3 * Math.sin(2 * Math.PI * f3 * t);
            }
        } else if (preset === 'ecg') {
            // 5-Gaussian PQRST. The (μ, σ, A) values are tuned for a
            // recognizable beat at HR ~70 bpm; widths in seconds.
            const hr = params.hr;
            const period = 60 / hr;
            const beats = [
                { mu: -0.20, sigma: 0.025, a: 0.10 }, // P
                { mu: -0.04, sigma: 0.012, a: -0.15 }, // Q
                { mu:  0.00, sigma: 0.010, a: 1.20 }, // R
                { mu:  0.04, sigma: 0.012, a: -0.25 }, // S
                { mu:  0.20, sigma: 0.040, a: 0.30 }  // T
            ];
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                // phase relative to nearest R-peak in (-T/2, T/2]
                let phase = (t % period);
                if (phase > period / 2) phase -= period;
                let v = 0;
                for (const g of beats) {
                    const z = (phase - g.mu) / g.sigma;
                    v += g.a * Math.exp(-0.5 * z * z);
                }
                x[n] = v;
            }
        } else if (preset === 'ppg') {
            // Two-harmonic Fourier model — fundamental + 2nd at 0.4·A.
            // Real PPGs have a slow upstroke and dicrotic notch; this
            // is a coarse but spectrally-honest stand-in.
            const hr = params.hr;
            const f0 = hr / 60;
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                const phi = 2 * Math.PI * f0 * t;
                x[n] = -Math.cos(phi) - 0.4 * Math.cos(2 * phi - 0.5);
            }
        } else if (preset === 'saw') {
            const decay = params.sawDecay;
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                x[n] = Math.exp(-decay * t) * Math.sin(2 * Math.PI * 40 * t);
            }
        } else if (preset === 'chirp') {
            // Linear chirp: instantaneous freq sweeps f0 -> f1 across
            // T = N/Fs. Phase is integral of 2π·f(t).
            const f0 = params.chirpF0;
            const f1 = params.chirpF1;
            const k = (f1 - f0) / T;
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                const phase = 2 * Math.PI * (f0 * t + 0.5 * k * t * t);
                x[n] = Math.sin(phase);
            }
        } else if (preset === 'square') {
            // 5 odd harmonics of 10 Hz — band-limited Fourier square.
            const f0 = 10;
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                let v = 0;
                for (let kh = 1; kh <= 9; kh += 2) {
                    v += (1 / kh) * Math.sin(2 * Math.PI * kh * f0 * t);
                }
                x[n] = v * (4 / Math.PI);
            }
        } else if (preset === 'am') {
            const fc = params.amFc;
            const fm = params.amFm;
            for (let n = 0; n < N; n++) {
                const t = n * dt;
                x[n] = Math.cos(2 * Math.PI * fc * t)
                     * (1 + 0.5 * Math.cos(2 * Math.PI * fm * t));
            }
        }

        if (noise > 0) {
            for (let n = 0; n < N; n++) {
                x[n] += noise * boxMuller();
            }
        }

        return x;
    }

    // --- DPR-aware canvas sizing --------------------------------------
    function resizeCanvas(canvas) {
        const dpr = Math.min(2, window.devicePixelRatio || 1);
        const rect = canvas.getBoundingClientRect();
        const w = Math.max(2, Math.floor(rect.width));
        const h = Math.max(2, Math.floor(rect.height));
        if (canvas.width !== w * dpr || canvas.height !== h * dpr) {
            canvas.width = w * dpr;
            canvas.height = h * dpr;
        }
        return { dpr, w, h };
    }

    // --- Time-domain plot ---------------------------------------------
    function drawTime(rawX, windowedX, Fs, N) {
        const { dpr, w, h } = resizeCanvas(timeCanvas);
        const ctx = timeCtx;
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, w, h);

        // Compute symmetric y range from the larger of the two traces.
        let amax = 1e-9;
        for (let i = 0; i < N; i++) {
            const a = Math.abs(rawX[i]);
            if (a > amax) amax = a;
            const b = Math.abs(windowedX[i]);
            if (b > amax) amax = b;
        }
        amax *= 1.1;

        // Grid
        ctx.strokeStyle = _cc.grid();
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let i = 0; i <= 8; i++) {
            const x = (i / 8) * w;
            ctx.moveTo(x, 0); ctx.lineTo(x, h);
        }
        for (let i = 0; i <= 4; i++) {
            const y = (i / 4) * h;
            ctx.moveTo(0, y); ctx.lineTo(w, y);
        }
        ctx.stroke();

        // Zero line
        ctx.strokeStyle = _cc.axis();
        ctx.beginPath();
        ctx.moveTo(0, h / 2);
        ctx.lineTo(w, h / 2);
        ctx.stroke();

        const xToPx = (i) => (i / (N - 1)) * w;
        const yToPx = (v) => h / 2 - (v / amax) * (h / 2 - 4);

        // Raw signal (ghost)
        ctx.strokeStyle = _cc.frame();
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let i = 0; i < N; i++) {
            const px = xToPx(i);
            const py = yToPx(rawX[i]);
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Windowed signal
        ctx.strokeStyle = '#0a84ff';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        for (let i = 0; i < N; i++) {
            const px = xToPx(i);
            const py = yToPx(windowedX[i]);
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        ctx.fillStyle = _cc.label();
        ctx.font = '11px -apple-system, system-ui, sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        ctx.fillText(`Fs = ${Fs} Hz · N = ${N}`, 8, 6);
    }

    // --- Frequency-domain plot ----------------------------------------
    function drawFreq(specDb, freqs, peaks, Fs) {
        const { dpr, w, h } = resizeCanvas(freqCanvas);
        const ctx = freqCtx;
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, w, h);

        const N2 = specDb.length; // = N/2
        const fmax = Fs / 2;
        const dBmin = -80;
        const dBmax = 5;

        // Pick a "sane" x-axis label increment.
        let labelStep = 100;
        if (fmax <= 100) labelStep = 10;
        else if (fmax <= 250) labelStep = 25;
        else if (fmax <= 600) labelStep = 50;
        else labelStep = 100;

        const padL = 36, padR = 8, padT = 8, padB = 22;
        const plotW = w - padL - padR;
        const plotH = h - padT - padB;

        const fToPx = (f) => padL + (f / fmax) * plotW;
        const dBToPx = (d) => padT + (1 - (d - dBmin) / (dBmax - dBmin)) * plotH;

        // Grid
        ctx.strokeStyle = _cc.grid();
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let f = 0; f <= fmax + 0.001; f += labelStep) {
            const x = fToPx(f);
            ctx.moveTo(x, padT); ctx.lineTo(x, padT + plotH);
        }
        for (let d = dBmin; d <= dBmax; d += 20) {
            const y = dBToPx(d);
            ctx.moveTo(padL, y); ctx.lineTo(padL + plotW, y);
        }
        ctx.stroke();

        // Bars
        ctx.fillStyle = 'rgba(10,132,255,0.85)';
        const baseY = dBToPx(dBmin);
        for (let i = 1; i < N2; i++) {
            const x = fToPx(freqs[i]);
            const xNext = fToPx(freqs[Math.min(i + 1, N2 - 1)]);
            const bw = Math.max(1, xNext - x);
            const y = dBToPx(specDb[i]);
            ctx.fillRect(x - bw / 2, y, bw, Math.max(0, baseY - y));
        }

        // Axis labels
        ctx.fillStyle = _cc.label();
        ctx.font = '10px -apple-system, system-ui, sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        for (let f = 0; f <= fmax + 0.001; f += labelStep) {
            ctx.fillText(`${f}`, fToPx(f), padT + plotH + 4);
        }
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        for (let d = dBmin; d <= dBmax; d += 20) {
            ctx.fillText(`${d}`, padL - 4, dBToPx(d));
        }
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        ctx.fillStyle = _cc.labelDim();
        ctx.fillText('Hz', padL + plotW - 14, padT + plotH + 4);
        ctx.save();
        ctx.translate(10, padT + plotH / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.textAlign = 'center';
        ctx.fillText('dB', 0, 0);
        ctx.restore();

        // Peak markers
        ctx.strokeStyle = '#ff375f';
        ctx.fillStyle = '#ff375f';
        ctx.lineWidth = 1;
        ctx.font = '11px -apple-system, system-ui, sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        for (const p of peaks) {
            const x = fToPx(p.freq);
            const y = dBToPx(p.dB);
            ctx.beginPath();
            ctx.arc(x, y, 3, 0, Math.PI * 2);
            ctx.fill();
            ctx.beginPath();
            ctx.moveTo(x, y);
            ctx.lineTo(x, padT + plotH);
            ctx.stroke();
            const label = `${p.freq.toFixed(1)} Hz`;
            const ty = Math.max(padT + 10, y - 6);
            ctx.fillText(label, x, ty);
        }
    }

    // --- Magnitude spectrum in dB -------------------------------------
    function magnitudeDb(reIm) {
        const { re, im } = reIm;
        const N = re.length;
        const half = N >> 1;
        const mag = new Float64Array(half);
        let peak = 1e-12;
        for (let i = 0; i < half; i++) {
            const m = Math.sqrt(re[i] * re[i] + im[i] * im[i]);
            mag[i] = m;
            if (m > peak) peak = m;
        }
        const dB = new Float64Array(half);
        for (let i = 0; i < half; i++) {
            // Floor of -80 dB keeps the log from blowing up at silent bins.
            const d = 20 * Math.log10(mag[i] / peak + 1e-12);
            dB[i] = d < -80 ? -80 : d;
        }
        return dB;
    }

    // --- Read all params from the DOM ---------------------------------
    function readParams() {
        const preset = strVal('fft-preset', '3sin');
        const Fs = numVal('fft-fs', 1000);
        const N = parseInt(strVal('fft-nsamples', '1024'), 10);
        const win = strVal('fft-window', 'hann');
        const noise = numVal('fft-noise', 0.05);
        return {
            preset, Fs, N, window: win, noise,
            f1: numVal('fft-f1', 50),
            f2: numVal('fft-f2', 120),
            f3: numVal('fft-f3', 200),
            a1: numVal('fft-a1', 1),
            a2: numVal('fft-a2', 0.5),
            a3: numVal('fft-a3', 0.25),
            hr: numVal('fft-hr', 70),
            sawDecay: numVal('fft-saw-decay', 5),
            chirpF0: numVal('fft-chirp-f0', 5),
            chirpF1: numVal('fft-chirp-f1', 100),
            amFc: numVal('fft-am-fc', 60),
            amFm: numVal('fft-am-fm', 3)
        };
    }

    // --- Update slider value labels -----------------------------------
    function syncLabels(p) {
        setText('fft-fs-val', `${p.Fs}`);
        setText('fft-noise-val', p.noise.toFixed(2));
        setText('fft-f1-val', `${p.f1.toFixed(0)}`);
        setText('fft-f2-val', `${p.f2.toFixed(0)}`);
        setText('fft-f3-val', `${p.f3.toFixed(0)}`);
        setText('fft-a1-val', p.a1.toFixed(2));
        setText('fft-a2-val', p.a2.toFixed(2));
        setText('fft-a3-val', p.a3.toFixed(2));
        setText('fft-hr-val', `${p.hr.toFixed(0)}`);
        setText('fft-saw-decay-val', p.sawDecay.toFixed(1));
        setText('fft-chirp-f0-val', p.chirpF0.toFixed(0));
        setText('fft-chirp-f1-val', p.chirpF1.toFixed(0));
        setText('fft-am-fc-val', p.amFc.toFixed(0));
        setText('fft-am-fm-val', p.amFm.toFixed(1));
    }

    // --- Main compute + draw pass -------------------------------------
    let visible = true;
    function compute() {
        if (!visible) return;
        const p = readParams();
        syncLabels(p);

        const wrap = $('fft-preset-wrap');
        if (wrap) wrap.setAttribute('data-preset', p.preset);

        const N = p.N;
        const Fs = p.Fs;

        const raw = buildSignal(p);
        const wx = applyWindow(raw, p.window);

        const X = fft(wx);
        const dB = magnitudeDb(X);

        const half = N >> 1;
        const freqs = new Float64Array(half);
        for (let i = 0; i < half; i++) freqs[i] = i * Fs / N;

        const peaks = findPeaks(dB, freqs, 3);
        const thd = computeTHD(X, freqs, peaks);

        drawTime(raw, wx, Fs, N);
        drawFreq(dB, freqs, peaks, Fs);

        const fmt = (pk) => pk ? `${pk.freq.toFixed(2)} Hz` : '—';
        setText('fft-readout-peak1', fmt(peaks[0]));
        setText('fft-readout-peak2', fmt(peaks[1]));
        setText('fft-readout-peak3', fmt(peaks[2]));
        setText('fft-readout-thd', thd === null ? '—' : `${thd.toFixed(2)} %`);
    }

    // --- Wiring -------------------------------------------------------
    let debounceId = null;
    function scheduleCompute() {
        if (debounceId) clearTimeout(debounceId);
        debounceId = setTimeout(() => {
            debounceId = null;
            compute();
        }, 30);
    }

    function wireUI() {
        const ids = [
            'fft-preset', 'fft-fs', 'fft-nsamples', 'fft-window', 'fft-noise',
            'fft-f1', 'fft-f2', 'fft-f3', 'fft-a1', 'fft-a2', 'fft-a3',
            'fft-hr', 'fft-saw-decay',
            'fft-chirp-f0', 'fft-chirp-f1',
            'fft-am-fc', 'fft-am-fm'
        ];
        for (const id of ids) {
            const el = $(id);
            if (!el) continue;
            el.addEventListener('input', scheduleCompute);
            el.addEventListener('change', scheduleCompute);
        }

        window.addEventListener('resize', scheduleCompute);

        // Pause when offscreen — saves CPU on long pages.
        const section = document.getElementById('fft') || timeCanvas;
        if ('IntersectionObserver' in window) {
            const io = new IntersectionObserver((entries) => {
                for (const e of entries) {
                    visible = e.isIntersecting;
                    if (visible) scheduleCompute();
                }
            }, { threshold: 0.05 });
            io.observe(section);
        }

        // prefers-reduced-motion: only-on-input is already our policy,
        // so this is effectively a no-op — but listen anyway in case
        // we ever add an idle animation.
        if (window.matchMedia) {
            const mq = window.matchMedia('(prefers-reduced-motion: reduce)');
            if (mq.addEventListener) mq.addEventListener('change', scheduleCompute);
        }
    }

    wireUI();
    // First paint after layout settles.
    if (document.readyState === 'complete' || document.readyState === 'interactive') {
        scheduleCompute();
    } else {
        window.addEventListener('DOMContentLoaded', scheduleCompute, { once: true });
    }
})();
