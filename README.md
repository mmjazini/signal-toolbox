# Signal-Processing Toolbox

Three interactive, from-scratch signal-processing tools that run entirely in the
browser — no libraries, no frameworks, no build step. Each is a self-contained
`<canvas>` widget that recomputes live as you drag a slider.

**▶ Live demo: [mahdijazini.com/toolbox](https://www.mahdijazini.com/toolbox.html)**

Built by [Mahdi Jazini](https://www.mahdijazini.com) — PhD in bioengineering,
non-invasive cardiovascular monitoring. These are teaching versions of the same
math I use in the lab for cuff-based hemodynamics.

---

## The tools

### 1. Oscillometric BP reasoner — [`tools/bp.js`](tools/bp.js)
How an arm cuff turns pressure oscillations into SBP / MAP / DBP.

- Arterial volume modeled with a sigmoidal compliance curve
  `V(Ptm) = 1 / (1 + e^(−Ptm/c))`; the oscillogram envelope at cuff pressure `Pc`
  is `V(SBP−Pc) − V(DBP−Pc)`.
- **MAP is the envelope peak = zero mean transmural pressure = maximum
  compliance**, which in this symmetric model is exactly `(SBP+DBP)/2` (holds for
  any `c`).
- Four estimators you can switch between: fixed-ratio, maximum-amplitude,
  variable-ratio (amplitude-dependent, per Dhamotharan, Jazini et al., *Sci
  Reports* 2025), and derivative (slope-extrema of `dE/dPc`).
- Draws the arterial **volume oscillations** swelling toward MAP and shrinking at
  the pressure extremes — the time-domain picture behind the envelope.

### 2. FFT explorer — [`tools/fft.js`](tools/fft.js)
A radix-2 Cooley–Tukey FFT written from scratch.

- Iterative in-place transform with cached twiddle factors and bit-reversal.
- Window functions (rectangular, Hann, Hamming, Blackman, flat-top, Tukey),
  Box–Muller Gaussian noise, and **parabolic sub-bin peak interpolation**.
- dB-normalized spectrum, top-3 peak markers, and a THD readout
  (`√(Σ harmonic² / fundamental²) × 100%`).
- Presets: three sinusoids, synthetic ECG/PPG, damped SAW resonator, linear
  chirp, square wave, AM tone.

### 3. Pan–Tompkins R-peak detector — [`tools/rpeak.js`](tools/rpeak.js)
The classic QRS detector, running on a synthetic ECG with known ground truth.

- McSharry-style 5-Gaussian PQRST beat generator, with optional PVCs, HR jitter,
  baseline wander, EMG, and 50/60 Hz powerline noise.
- Full detection pipeline: 5–15 Hz band-pass → 5-point derivative → squaring →
  150 ms moving-window integration → adaptive dual thresholds (SPKI/NPKI) with
  refractory period, T-wave discrimination, and search-back.
- Sensitivity / PPV computed against the **generator's ground-truth** R-times
  (not the detector's own output), plus HRV metrics (mean RR, RMSSD, pNN50).

---

## Design notes

- **Vanilla JS, zero dependencies.** Each tool is a single IIFE that renders to
  its own canvases and bails silently if they aren't present, so any tool can be
  dropped onto any page in isolation.
- **Everything is computed, not faked.** The FFT is a real FFT; the detector is a
  real adaptive-threshold pipeline; the BP envelope comes from the compliance
  model, not a drawn curve.
- The tools are theme-aware (read `data-theme` at draw time) and pause their
  animation loops when scrolled off-screen.

## Run locally
No build step. Serve the folder and open an HTML page that hosts the canvases —
or just see them wired up on the [live demo](https://www.mahdijazini.com/toolbox.html).
```
python -m http.server 8000
```

## License
MIT — see [LICENSE](LICENSE).
