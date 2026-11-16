// Interactive 1D Schrodinger / Gross-Pitaevskii equation (hbar = m = 1).
// Crank-Nicolson stepping, (1 + cH) psi' = (1 - cH) psi, with c = i dt/2 for
// real time and c = dtau/2 for imaginary time (plus renormalisation, which
// relaxes any state to the ground state). The tridiagonal system is factorised
// once per potential/mode. The grid extends beyond the visible window so the
// complex absorbing potential at the edges is never on screen.
//
// Every run follows the same recipe as an experiment:
//   prepare (a wavepacket, or the ground state of some potential at some g0)
//   -> optionally imprint (momentum kick, dark soliton, density bump)
//   -> evolve with the current potential and interaction g.
(function () {
    const root = document.getElementById("tdse-widget");
    if (!root) return;

    // ------------------------------------------------------------------
    // grid
    // ------------------------------------------------------------------
    let VMIN = -40, VMAX = 40;             // visible window (a ring shows all of it)
    const XMIN = -56, XMAX = 56;           // computational domain (psi = 0 at the ends)
    const DEFAULT_DX = 0.1, DEFAULT_DT = 0.01;
    const WALL = 200;                      // "infinite" wall height for the box
    const CAP_START = 44, CAP_STRENGTH = 1.5;
    const K_MAX = 4;

    // Everything that lives on the grid. Changing dx (in the numerics settings)
    // reallocates all of it; nothing else needs to know.
    let dx = DEFAULT_DX, dt = DEFAULT_DT, N = 0, RING_L = 0;
    let x;                                  // grid points
    let pR, pI;                             // psi
    let V, Ve, W;                           // potential, what the solver sees (V + g|psi|^2), absorber (H = ... - iW)
    let mR, mI, cR, cI;                     // 1/pivots and modified super-diagonal of the factorisation
    let rR, rI, p0R, p0I;                   // work arrays; psi at the start of a nonlinear step
    let e0R, e0I, zR, zI;                   // an eigenstate to orthogonalise against; ring correction vector
    let j0 = 0, j1 = 0;                     // visible index range

    function setView(min, max) {
        VMIN = min;
        VMAX = max;
        j0 = 0;
        j1 = N - 1;
        while (x[j0] < VMIN) j0++;
        while (x[j1] > VMAX) j1--;
    }

    function allocateGrid(spacing) {
        dx = spacing;
        N = Math.round((XMAX - XMIN) / dx) - 1;
        RING_L = N * dx;
        const arr = () => new Float64Array(N);
        x = arr();
        for (let j = 0; j < N; j++) x[j] = XMIN + (j + 1) * dx;
        [pR, pI, V, Ve, W, mR, mI, cR, cI, rR, rI, p0R, p0I, e0R, e0I, zR, zI] = Array.from({ length: 17 }, arr);
        setView(VMIN, VMAX);
    }
    allocateGrid(DEFAULT_DX);

    let imaginary = false;
    let g = 0;                                                  // mean-field interaction strength
    let sR = 0, sI = 0;                                         // off-diagonal of LHS

    // ------------------------------------------------------------------
    // Crank-Nicolson (+ predictor-corrector for the Gross-Pitaevskii term)
    // ------------------------------------------------------------------
    let relaxScale = 1;   // larger imaginary-time steps while preparing eigenstates

    function coefficient() {
        return imaginary ? [0.5 * dt * relaxScale, 0] : [0, 0.5 * dt];
    }

    // Diagonal of (1 + cH) at grid point j.
    function diagonal(j, kr, ki) {
        const hR = 1 / (dx * dx) + Ve[j], hI = imaginary ? 0 : -W[j];
        return [1 + kr * hR - ki * hI, kr * hI + ki * hR];
    }

    // On a ring the first and last points are neighbours, which adds two corner
    // elements s to the tridiagonal matrix A. Sherman-Morrison handles them:
    // A = A' + u v^T with u = (gamma, 0, ..., 0, s), v = (1, 0, ..., 0, s/gamma),
    // gamma = -b_0, and A' tridiagonal with b_0 -> 2 b_0, b_N-1 -> b_N-1 + s^2/b_0.
    // Then A^-1 r = y - (v.y / (1 + v.z)) z, with y = A'^-1 r and z = A'^-1 u.
    let periodic = false;
    const ring = { vR: 0, vI: 0, fR: 0, fI: 0 };   // v_N-1 = s/gamma, f = 1/(1 + v.z)

    const cdiv = (aR, aI, bR, bI) => {
        const n = bR * bR + bI * bI;
        return [(aR * bR + aI * bI) / n, (aI * bR - aR * bI) / n];
    };

    function factorize() {
        const [kr, ki] = coefficient();
        const off = -0.5 / (dx * dx);
        sR = kr * off;
        sI = ki * off;
        const [b0R, b0I] = diagonal(0, kr, ki);
        let pcR = 0, pcI = 0;
        for (let j = 0; j < N; j++) {
            let [bR, bI] = diagonal(j, kr, ki);
            if (periodic && j === 0) {
                bR *= 2;
                bI *= 2;
            } else if (periodic && j === N - 1) {
                const [qR, qI] = cdiv(sR * sR - sI * sI, 2 * sR * sI, b0R, b0I);
                bR += qR;
                bI += qI;
            }
            // pivot = b_j - s c'_{j-1}
            const dR = bR - (sR * pcR - sI * pcI);
            const dI = bI - (sR * pcI + sI * pcR);
            const n = dR * dR + dI * dI;
            mR[j] = dR / n;
            mI[j] = -dI / n;
            cR[j] = sR * mR[j] - sI * mI[j];
            cI[j] = sR * mI[j] + sI * mR[j];
            pcR = cR[j];
            pcI = cI[j];
        }
        if (!periodic) return;
        // z = A'^-1 u, and the constants of the correction
        zR.fill(0);
        zI.fill(0);
        zR[0] = -b0R;
        zI[0] = -b0I;
        zR[N - 1] = sR;
        zI[N - 1] = sI;
        thomas(zR, zI, zR, zI);
        [ring.vR, ring.vI] = cdiv(sR, sI, -b0R, -b0I);
        const dR = 1 + zR[0] + (ring.vR * zR[N - 1] - ring.vI * zI[N - 1]);
        const dI = zI[0] + (ring.vR * zI[N - 1] + ring.vI * zR[N - 1]);
        [ring.fR, ring.fI] = cdiv(1, 0, dR, dI);
    }

    // the linear problem only needs factorising when the potential or mode changes
    function refresh() {
        if (g === 0) Ve.set(V);
        factorize();
    }

    // Solve A' out = in with the current factorisation (in is used as scratch).
    function thomas(inR, inI, outR, outI) {
        // forward sweep: q_j = (r_j - s q_{j-1}) m_j
        let qR = 0, qI = 0;
        for (let j = 0; j < N; j++) {
            const aR = inR[j] - (sR * qR - sI * qI), aI = inI[j] - (sR * qI + sI * qR);
            qR = aR * mR[j] - aI * mI[j];
            qI = aR * mI[j] + aI * mR[j];
            inR[j] = qR;
            inI[j] = qI;
        }
        // back substitution: x_j = q_j - c'_j x_{j+1}
        outR[N - 1] = inR[N - 1];
        outI[N - 1] = inI[N - 1];
        for (let j = N - 2; j >= 0; j--) {
            outR[j] = inR[j] - (cR[j] * outR[j + 1] - cI[j] * outI[j + 1]);
            outI[j] = inI[j] - (cR[j] * outI[j + 1] + cI[j] * outR[j + 1]);
        }
    }

    // neighbours of j (zero beyond the ends, unless the grid is a ring)
    const left = (a, j) => (j > 0 ? a[j - 1] : periodic ? a[N - 1] : 0);
    const right = (a, j) => (j < N - 1 ? a[j + 1] : periodic ? a[0] : 0);

    let absorbedLeft = 0, absorbedRight = 0;

    // one Crank-Nicolson solve with the current factorisation: psi = (1 + cH)^-1 (1 - cH) src
    function solve(srcR, srcI, countLoss) {
        const [kr, ki] = coefficient();
        for (let j = 0; j < N; j++) {
            const hR = 1 / (dx * dx) + Ve[j], hI = imaginary ? 0 : -W[j];
            const gR = 1 - (kr * hR - ki * hI), gI = -(kr * hI + ki * hR);
            const nR = left(srcR, j) + right(srcR, j);
            const nI = left(srcI, j) + right(srcI, j);
            rR[j] = gR * srcR[j] - gI * srcI[j] - (sR * nR - sI * nI);
            rI[j] = gR * srcI[j] + gI * srcR[j] - (sR * nI + sI * nR);
            if (countLoss && !imaginary && W[j] > 0) {
                const loss = 2 * W[j] * (srcR[j] * srcR[j] + srcI[j] * srcI[j]) * dt * dx;
                if (x[j] < 0) absorbedLeft += loss; else absorbedRight += loss;
            }
        }
        thomas(rR, rI, pR, pI);
        if (!periodic) return;
        // psi -= (v.y) f z
        const yR = pR[0] + (ring.vR * pR[N - 1] - ring.vI * pI[N - 1]);
        const yI = pI[0] + (ring.vR * pI[N - 1] + ring.vI * pR[N - 1]);
        const aR = yR * ring.fR - yI * ring.fI, aI = yR * ring.fI + yI * ring.fR;
        for (let j = 0; j < N; j++) {
            pR[j] -= aR * zR[j] - aI * zI[j];
            pI[j] -= aR * zI[j] + aI * zR[j];
        }
    }

    // Time-dependent potentials (a few presets only, e.g. a moving obstacle): the
    // potential is evaluated at the middle of each step, which keeps
    // Crank-Nicolson second-order accurate, and the system is re-factorised.
    // A preset may also modulate the interaction (gMod) or carry its own state that
    // moves on each step (advance), e.g. an obstacle following the mouse.
    let tNow = 0;
    let gScale = 1;   // g(t) = g * gScale
    const stir = { pos: -25, target: -25, vel: 0 };
    function applyTime(t) {
        const pr = presets[presetKey];
        if (!pr || !pr.timeDep) return;
        for (let j = 0; j < N; j++) V[j] = pr.V(x[j], curA, curB, t);
        capPotential();
        gScale = pr.gMod ? pr.gMod(t, curA, curB) : 1;
        if (g === 0) {
            Ve.set(V);
            factorize();
        }
    }

    function step() {
        const pr = presets[presetKey];
        const timeDep = !imaginary && pr && pr.timeDep;
        if (timeDep && pr.advance) pr.advance(dt);
        if (timeDep) applyTime(tNow + 0.5 * dt);
        const gs = g * (timeDep ? gScale : 1);
        if (g === 0) {
            solve(pR, pI, true);
        } else {
            // predictor: step with the density at the start of the step
            p0R.set(pR);
            p0I.set(pI);
            for (let j = 0; j < N; j++) Ve[j] = V[j] + gs * (p0R[j] * p0R[j] + p0I[j] * p0I[j]);
            factorize();
            solve(p0R, p0I, false);
            if (imaginary) normalize();
            // corrector: redo it with the density averaged over the step
            for (let j = 0; j < N; j++) {
                const r0 = p0R[j] * p0R[j] + p0I[j] * p0I[j], r1 = pR[j] * pR[j] + pI[j] * pI[j];
                Ve[j] = V[j] + 0.5 * gs * (r0 + r1);
            }
            factorize();
            solve(p0R, p0I, true);
        }
        if (imaginary) normalize();
        if (timeDep) tNow += dt;
    }

    function normalize() {
        let norm = 0;
        for (let j = 0; j < N; j++) norm += pR[j] * pR[j] + pI[j] * pI[j];
        const s = 1 / Math.sqrt(norm * dx);
        for (let j = 0; j < N; j++) {
            pR[j] *= s;
            pI[j] *= s;
        }
    }

    // Relax psi to the lowest state of the current H in imaginary time (kept
    // orthogonal to `orth` if given), behind the scenes. Imaginary-time Crank-
    // Nicolson converges to the exact same eigenstates for any step size, so a
    // larger step makes this fast.
    function relax(maxSteps, orth = null, pin = null) {
        const wasImaginary = imaginary;
        imaginary = true;
        relaxScale = 5;
        refresh();
        let last = Infinity;
        const n = Math.round((maxSteps * DEFAULT_DT) / dt);
        for (let i = 0; i < n; i++) {
            step();
            if (orth) {
                let cR = 0, cI = 0;   // <orth|psi>
                for (let j = 0; j < N; j++) {
                    cR += orth[0][j] * pR[j] + orth[1][j] * pI[j];
                    cI += orth[0][j] * pI[j] - orth[1][j] * pR[j];
                }
                cR *= dx;
                cI *= dx;
                for (let j = 0; j < N; j++) {
                    pR[j] -= cR * orth[0][j] - cI * orth[1][j];
                    pI[j] -= cR * orth[1][j] + cI * orth[0][j];
                }
                normalize();
            }
            if (pin !== null) {
                // keep a real state with its node (sign change) fixed at x = pin
                for (let j = 0; j < N; j++) {
                    pR[j] = Math.abs(pR[j]) * (x[j] < pin ? -1 : 1);
                    pI[j] = 0;
                }
                normalize();
            }
            if (i % 50 === 49) {
                const e = measure().E;
                if (Math.abs(e - last) < 1e-10) break;
                last = e;
            }
        }
        relaxScale = 1;
        imaginary = wasImaginary;
        refresh();
    }

    // energy per particle (incl. mean-field term g/2 |psi|^4), kinetic energy and norm
    function measure() {
        let kin = 0, pot = 0, int = 0, norm = 0, xm = 0, vis = 0, flow = 0;
        for (let j = 0; j < N; j++) {
            const lR = left(pR, j) + right(pR, j) - 2 * pR[j];
            const lI = left(pI, j) + right(pI, j) - 2 * pI[j];
            const r = pR[j] * pR[j] + pI[j] * pI[j];
            kin += -0.5 * (pR[j] * lR + pI[j] * lI) / (dx * dx);
            pot += V[j] * r;
            int += 0.5 * g * r * r;
            norm += r;
            xm += x[j] * r;
            if (j >= j0 && j <= j1) vis += r;
            // local wavenumber = phase advance per grid cell; its grid speed error,
            // weighted by density, is how fast the flowing parts drift. Standing
            // (real) states have none, however sharp they are.
            if (j < N - 1 && r > 1e-12) {
                const re = pR[j] * pR[j + 1] + pI[j] * pI[j + 1], im = pR[j] * pI[j + 1] - pI[j] * pR[j + 1];
                flow += r * gridSpeedError(Math.abs(Math.atan2(im, re)) / dx);
            }
        }
        return { E: (kin + pot + int) / norm, T: kin / norm, norm: norm * dx, visible: vis * dx, x: xm / norm, drift: flow / norm };
    }

    // <p> per particle (central differences, ring-aware)
    function meanMomentum() {
        let pm = 0, nn = 0;
        for (let j = 0; j < N; j++) {
            pm += pR[j] * (right(pI, j) - left(pI, j)) - pI[j] * (right(pR, j) - left(pR, j));
            nn += pR[j] * pR[j] + pI[j] * pI[j];
        }
        return pm / (2 * dx) / nn;
    }

    // how many times the phase winds around the ring (an integer by construction)
    function winding() {
        let total = 0;
        for (let j = 0; j < N; j++) {
            const k = (j + 1) % N;
            total += Math.atan2(pR[j] * pI[k] - pI[j] * pR[k], pR[j] * pR[k] + pI[j] * pI[k]);
        }
        return Math.round(total / (2 * Math.PI));
    }

    // |phi(k)|^2 at the wavenumbers ks, phi(k) = sum_j psi_j e^{-i k x_j} dx, with the
    // phase factor advanced by a rotation from point to point (no trig in the loop)
    function momentumDistribution(ks, out) {
        for (let m = 0; m < ks.length; m++) {
            const k = ks[m], cs = Math.cos(k * dx), sn = -Math.sin(k * dx);
            let eR = Math.cos(k * x[0]), eI = -Math.sin(k * x[0]), aR = 0, aI = 0;
            for (let j = 0; j < N; j++) {
                aR += pR[j] * eR - pI[j] * eI;
                aI += pR[j] * eI + pI[j] * eR;
                const t = eR * cs - eI * sn;
                eI = eR * sn + eI * cs;
                eR = t;
            }
            out[m] = (aR * aR + aI * aI) * dx * dx / (2 * Math.PI);
        }
    }

    function potentialAt(xq) {
        const f = (xq - XMIN) / dx - 1;
        const j = Math.max(0, Math.min(N - 2, Math.floor(f)));
        const t = Math.max(0, Math.min(1, f - j));
        return V[j] * (1 - t) + V[j + 1] * t;
    }

    // ring geometry: circumference, density of a uniform condensate, and the
    // (signed) shortest distance between two points on it
    const ringMu = () => Math.max(g, 1) / RING_L;
    const ringGap = (d) => d - RING_L * Math.round(d / RING_L);
    const ringReadout = () => `Circulation: w = ${winding()}`;

    // rms width of the cloud
    function cloudWidth() {
        let n = 0, m = 0, m2 = 0;
        for (let j = 0; j < N; j++) {
            const r = pR[j] * pR[j] + pI[j] * pI[j];
            n += r;
            m += r * x[j];
            m2 += r * x[j] * x[j];
        }
        m /= n;
        return Math.sqrt(Math.max(0, m2 / n - m * m));
    }

    // Opening a trap from omega0 to omegaf in a time T. The shortcut follows a
    // scaling solution: the width grows as b(t), a polynomial with b' = b'' = 0 at
    // both ends, and the Ermakov equation b'' + omega^2 b = omega0^2 / b^3 gives the
    // frequency that produces it (omega^2 even dips below zero for fast ramps).
    // s mixes it with a plain linear ramp: s = 1 shortcut, s = 0 plain.
    const RAMP = { w0: 0.5, wf: 0.05 };
    function rampOmega2(t, T, s) {
        const { w0, wf } = RAMP;
        if (t >= T) return wf * wf;
        const gam = Math.sqrt(w0 / wf), u = t / T;
        const b = 1 + (gam - 1) * (10 * u ** 3 - 15 * u ** 4 + 6 * u ** 5);
        const bpp = ((gam - 1) * (60 * u - 180 * u * u + 120 * u ** 3)) / (T * T);
        const shortcut = (w0 * w0) / b ** 4 - bpp / b, plain = (w0 + (wf - w0) * u) ** 2;
        return (1 - s) * plain + s * shortcut;
    }

    // the three Bragg pulses (pi/2, pi, pi/2) of the atom interferometer
    const BRAGG = { tp: Math.PI, T: 12 };
    function braggStage(t) {
        const { tp, T } = BRAGG;
        if (t < tp) return "Pulse 1: beam splitter";
        if (t < tp + T) return "The halves drift apart";
        if (t < 3 * tp + T) return "Pulse 2: mirror";
        if (t < 3 * tp + 2 * T) return "The halves meet again";
        if (t < 4 * tp + 2 * T) return "Pulse 3: beam splitter";
        const moving = Math.max(0, Math.min(1, meanMomentum() / 2));
        return `Out: ${(100 * (1 - moving)).toFixed(0)}% at rest, ${(100 * moving).toFixed(0)}% moving`;
    }

    // ------------------------------------------------------------------
    // potentials (energy scale fixed per preset so the view never jumps).
    // `bound`: the potential has a proper ground state on screen.
    // ------------------------------------------------------------------
    const presets = {
        barrier: {
            label: "Single barrier", absorbing: true, flux: true, scale: 5,
            eq: String.raw`V(x) = \begin{cases} V_0 & |x| < a/2 \\ 0 & \text{otherwise} \end{cases}`,
            p1: { name: "Height", sym: "V_0", min: 0, max: 6, step: 0.1, value: 2.4 },
            p2: { name: "Width", sym: "a", min: 0.2, max: 6, step: 0.1, value: 0.8 },
            V: (xx, h, w) => (Math.abs(xx) < w / 2 ? h : 0),
            packet: { x0: -15, k: 1.8, sigma: 3 },
        },
        double: {
            label: "Double barrier", absorbing: true, flux: true, scale: 5,
            eq: String.raw`V(x) = \begin{cases} V_0 & b/2 < |x| < b/2 + 0.6 \\ 0 & \text{otherwise} \end{cases}`,
            p1: { name: "Height", sym: "V_0", min: 0, max: 6, step: 0.1, value: 3 },
            p2: { name: "Gap", sym: "b", min: 1, max: 12, step: 0.1, value: 4 },
            V: (xx, h, gap) => (Math.abs(Math.abs(xx) - gap / 2 - 0.3) < 0.3 ? h : 0),
            packet: { x0: -18, k: 1.6, sigma: 4 },
        },
        step: {
            label: "Step", absorbing: true, flux: true, scale: 5,
            eq: String.raw`V(x) = \begin{cases} V_0 & x > 0 \\ 0 & x < 0 \end{cases}`,
            p1: { name: "Height", sym: "V_0", min: 0, max: 6, step: 0.1, value: 1.6 },
            V: (xx, h) => (xx > 0 ? h : 0),
            packet: { x0: -15, k: 2, sigma: 2.5 },
        },
        well: {
            label: "Square well", absorbing: true, flux: true, scale: 5, holds: true,
            eq: String.raw`V(x) = \begin{cases} 0 & |x| < a/2 \\ U & \text{otherwise} \end{cases}`,
            p1: { name: "Depth", sym: "U", min: 0, max: 5, step: 0.1, value: 2 },
            p2: { name: "Width", sym: "a", min: 0.5, max: 10, step: 0.1, value: 3 },
            // a dip of depth U in an otherwise flat floor (the energy zero is irrelevant)
            V: (xx, U, a) => (Math.abs(xx) < a / 2 ? 0 : U),
            packet: { x0: -24, k: 0.85, sigma: 7 },
        },
        free: {
            label: "Free space", absorbing: true, flux: false, scale: 5,
            eq: String.raw`V(x) = 0`,
            V: () => 0,
            packet: { x0: -15, k: 1.5, sigma: 1.2 },
        },
        harmonic: {
            label: "Harmonic trap", absorbing: false, flux: false, scale: 5, bound: true,
            eq: String.raw`V(x) = \tfrac{1}{2}\,\omega^2 x^2`,
            p1: { name: "Frequency", sym: String.raw`\omega`, min: 0.05, max: 0.6, step: 0.01, value: 0.2 },
            V: (xx, w) => 0.5 * w * w * xx * xx,
            packet: { x0: -10, k: 0, sigma: 1.6 },
        },
        box: {
            label: "Box", absorbing: false, flux: false, scale: 5, bound: true,
            eq: String.raw`V(x) = \begin{cases} 0 & |x| < L/2 \\ \infty & \text{otherwise} \end{cases}`,
            p1: { name: "Width", sym: "L", min: 6, max: 76, step: 1, value: 30 },
            V: (xx, L) => (Math.abs(xx) < L / 2 ? 0 : WALL),
            guess: (xx, L) => (Math.abs(xx) < L / 2 ? 1 : 0),
            packet: { x0: -5, k: 2, sigma: 2 },
        },
        doublewell: {
            label: "Double well", absorbing: false, flux: false, scale: 2.5, bound: true,
            eq: String.raw`V(x) = h\,\big[(2x/s)^2 - 1\big]^2`,
            p1: { name: "Barrier", sym: "h", min: 0.1, max: 2, step: 0.05, value: 0.3 },
            p2: { name: "Separation", sym: "s", min: 6, max: 30, step: 0.5, value: 10 },
            V: (xx, h, s) => h * Math.pow((xx / (s / 2)) ** 2 - 1, 2),
            // roughly the left well's own ground state, so it tunnels cleanly
            packet: { x0: -5, k: 0, sigma: 1.2 },
        },
        lattice: {
            label: "Tilted lattice", absorbing: true, flux: false, scale: 7,
            eq: String.raw`V(x) = V_0 \sin^2(\pi x/2) + F x`,
            p1: { name: "Depth", sym: "V_0", min: 0, max: 4, step: 0.1, value: 1.5 },
            p2: { name: "Tilt", sym: "F", min: 0, max: 0.15, step: 0.005, value: 0.05 },
            // lattice spacing 2, plus a constant force F pushing towards -x
            V: (xx, V0, F) => V0 * Math.sin((Math.PI * xx) / 2) ** 2 + F * (xx + 40),
            // packets are first settled briefly in the untilted lattice (as atoms are
            // loaded in experiments), which removes their higher-band components
            untilted: (xx, V0) => V0 * Math.sin((Math.PI * xx) / 2) ** 2,
            settle: 1,
            kRange: 3,
            packet: { x0: 8, k: 0, sigma: 5 },
        },
        speckle: {
            label: "Random speckle", absorbing: true, flux: false, scale: 5,
            eq: String.raw`V(x) = S\, I_\ell(x)`,
            eqSub: String.raw`\(I_\ell\): a fixed random speckle pattern of mean 1 and grain size \(\ell\)`,
            p1: { name: "Strength", sym: "S", min: 0, max: 3, step: 0.1, value: 1 },
            p2: { name: "Grain", sym: String.raw`\ell`, min: 0.4, max: 3, step: 0.1, value: 1 },
            field: (s, grain) => speckleField(grain).map((v) => s * v),
            packet: { x0: -1.6, k: 0, sigma: 3 },
        },
        obstacle: {
            label: "Moving obstacle", absorbing: false, flux: false, scale: 0.8, timeDep: true, bound: true, needsCondensate: true,
            // a small bump dragged through a condensate in a box of width 70, from
            // x = -28 to +28. Speed in units of the sound speed c = sqrt(g n), height
            // in units of the chemical potential mu = g n, with n = 1/70.
            eq: String.raw`V(x, t) = U\, e^{-(x - x_0 - v t)^2/2w^2}`,
            eqSub: String.raw`in a box \(|x| < 35\), with \(w = 1.5\), sound speed \(c = \sqrt{g n}\) and \(\mu = g n\)`,
            p1: { name: "Speed", sym: "v/c", min: 0, max: 1.5, step: 0.05, value: 0.3 },
            p2: { name: "Height", sym: String.raw`U/\mu`, min: 0, max: 1, step: 0.05, value: 0.3 },
            V: (xx, vf, uf, t) => {
                if (Math.abs(xx) >= 35) return WALL;
                const mu = Math.max(g, 1) / 70, pos = Math.min(28, -28 + vf * Math.sqrt(mu) * t);
                return uf * mu * Math.exp(-((xx - pos) ** 2) / (2 * 1.5 * 1.5));
            },
            guess: (xx) => (Math.abs(xx) < 35 ? 1 : 0),
            packet: { x0: 0, k: 0, sigma: 3 },
        },
        // ---- driven: the user only gets the one or two knobs that keep them well behaved
        shaken: {
            label: "Shaken double well", absorbing: false, flux: false, scale: 2.5, bound: true, timeDep: true,
            eq: String.raw`V(x, t) = h\,\big[(2x/s)^2 - 1\big]^2 + A\, x \sin \Omega t`,
            eqSub: String.raw`with \(h = 0.3\), \(s = 8\) and \(\Omega = 0.08\): the wells are tilted back and forth`,
            p1: { name: "Drive", sym: "A", min: 0, max: 0.05, step: 0.001, value: 0.029 },
            V: (xx, A, _, t) => 0.3 * ((xx / 4) ** 2 - 1) ** 2 + A * xx * Math.sin(0.08 * t),
            packet: { x0: -4, k: 0, sigma: 1 },
        },
        sweep: {
            label: "Swept double well", absorbing: false, flux: false, scale: 2.5, bound: true, timeDep: true,
            eq: String.raw`V(x, t) = h\,\big[(2x/s)^2 - 1\big]^2 + F(t)\, x`,
            eqSub: String.raw`with \(h = 0.3\), \(s = 8\); the tilt \(F\) is swept from \(+0.015\) to \(-0.015\) over a time \(T\)`,
            p1: { name: "Sweep time", sym: "T", min: 20, max: 800, step: 10, value: 600 },
            V: (xx, T, _, t) => 0.3 * ((xx / 4) ** 2 - 1) ** 2 + (t < T ? 0.015 * (1 - (2 * t) / T) : -0.015) * xx,
            packet: { x0: -4, k: 0, sigma: 1 },
        },
        paul: {
            label: "Paul trap", absorbing: true, flux: false, scale: 8, timeDep: true,
            eq: String.raw`V(x, t) = V_0 + \big(A \cos \Omega t - \tfrac{\kappa}{2}\big)\, x^2`,
            eqSub: String.raw`for \(|x| < 20\), flat beyond, with \(V_0 = 4\), \(\kappa = 0.02\), \(\Omega = 2\): a hill whose curvature is shaken`,
            // beyond ~0.6 the micromotion gets violent and the result depends on dt
            p1: { name: "Drive", sym: "A", min: 0, max: 0.6, step: 0.02, value: 0.4 },
            V: (xx, A, _, t) => 4 + (A * Math.cos(2 * t) - 0.01) * Math.min(xx * xx, 400),
            packet: { x0: 3, k: 0, sigma: 1.2 },
        },
        bragg: {
            label: "Bragg pulses", absorbing: true, flux: false, scale: 2, timeDep: true,
            eq: String.raw`V(x, t) = V_0\, f(t) \cos(2x - 2t + \phi_n)`,
            eqSub: String.raw`three pulses of a moving lattice (\(V_0 = 0.5\)): a \(\pi/2\) pulse, a \(\pi\) pulse after \(T = 12\), and a \(\pi/2\) pulse with phase \(\phi\)`,
            p1: { name: "Phase", sym: String.raw`\phi`, min: 0, max: 6.28, step: 0.02, value: 0 },
            readout: () => braggStage(tNow),
            kRange: 4,
            V: (xx, phi, _, t) => {
                const { tp, T } = BRAGG;   // pi/2-pulse length pi/(2 V0), free time T
                const on = t < tp ? 0 : t >= tp + T && t < 3 * tp + T ? 0 : t >= 3 * tp + 2 * T && t < 4 * tp + 2 * T ? phi : null;
                // (+0.5 keeps the drawn potential above zero; a constant shift is only a global phase)
                return on === null ? 0 : 0.5 + 0.5 * Math.cos(2 * xx - 2 * t + on);
            },
            packet: { x0: -25, k: 0, sigma: 6 },
        },
        belt: {
            label: "Conveyor belt", absorbing: true, flux: false, scale: 18, timeDep: true,
            eq: String.raw`V(x, t) = V_0 \sin^2\!\big(\pi\,(x - x_L(t))/2\big)`,
            eqSub: String.raw`the lattice is moved smoothly by 30 over a time \(T\), then stops`,
            p1: { name: "Depth", sym: "V_0", min: 0, max: 16, step: 0.5, value: 14 },
            p2: { name: "Transport time", sym: "T", min: 20, max: 120, step: 5, value: 60 },
            V: (xx, V0, T, t) => {
                const xL = t < T ? 15 * (1 - Math.cos((Math.PI * t) / T)) : 30;
                return V0 * Math.sin((Math.PI * (xx - xL)) / 2) ** 2;
            },
            untilted: (xx, V0) => V0 * Math.sin((Math.PI * xx) / 2) ** 2,
            settle: 3,
            packet: { x0: -15, k: 0, sigma: 4 },
        },
        stir: {
            label: "Stirred condensate", absorbing: false, flux: false, scale: 0.8, timeDep: true, bound: true, needsCondensate: true, stir: true,
            eq: String.raw`V(x, t) = U\, e^{-(x - x_{\rm ob}(t))^2/2w^2}`,
            eqSub: String.raw`in a box \(|x| < 35\), \(w = 1.5\), \(\mu = g n\). Drag on the plot to move the obstacle`,
            p1: { name: "Height", sym: String.raw`U/\mu`, min: 0, max: 1, step: 0.05, value: 0.4 },
            V: (xx, uf) => {
                if (Math.abs(xx) >= 35) return WALL;
                const mu = Math.max(g, 1) / 70;
                return uf * mu * Math.exp(-((xx - stir.pos) ** 2) / (2 * 1.5 * 1.5));
            },
            guess: (xx) => (Math.abs(xx) < 35 ? 1 : 0),
            packet: { x0: 0, k: 0, sigma: 3 },
            // the obstacle follows the pointer with a short lag; its speed is smoothed
            advance: (step) => {
                const prev = stir.pos;
                stir.pos += (stir.target - stir.pos) * (1 - Math.exp(-step / 0.4));
                stir.vel += ((stir.pos - prev) / step - stir.vel) * (1 - Math.exp(-step / 1.5));
            },
            readout: () => `Obstacle speed: ${(Math.abs(stir.vel) / Math.sqrt(Math.max(g, 1) / 70)).toFixed(2)} c`,
        },
        pulse: {
            label: "Light pulse", absorbing: true, flux: false, scale: 8, timeDep: true, kRange: 3,
            eq: String.raw`V(x, t) = V_0 \sin^2 (x/2) \ \text{ for } t < \tau`,
            eqSub: String.raw`a standing wave of light flashed on for \(\tau = 0.1\); its strength is the pulse area \(V_0 \tau\)`,
            p1: { name: "Pulse area", sym: String.raw`V_0 \tau`, min: 0, max: 3, step: 0.1, value: 3 },
            V: (xx, area, _, t) => (t < 0.1 ? (area / 0.1) * Math.sin(xx / 2) ** 2 : 0),
            packet: { x0: 0, k: 0, sigma: 2.5 },
        },
        lens: {
            label: "Lens pulse", absorbing: true, flux: false, scale: 2, timeDep: true, kRange: 2,
            eq: String.raw`V(x, t) = \tfrac{1}{2}\kappa x^2 \ \text{ for } t_1 < t < t_1 + \tau`,
            eqSub: String.raw`a harmonic trap flashed on at \(t_1 = 12\) for \(\tau = 0.5\), with power \(f\): \(\kappa = f/(t_1 \tau)\)`,
            p1: { name: "Lens power", sym: "f", min: 0, max: 2, step: 0.05, value: 1 },
            V: (xx, f, _, t) => (t >= 12 && t < 12.5 ? (f * xx * xx) / 12 : 0),
            packet: { x0: 0, k: 0, sigma: 1 },
        },
        faraday: {
            label: "Pulsating interactions", absorbing: false, flux: false, scale: 0.8, timeDep: true, bound: true, needsCondensate: true,
            eq: String.raw`g(t) = g\,\big(1 + \epsilon \sin \Omega t\big)`,
            eqSub: String.raw`in a box \(|x| < 30\). Experiments do this by squeezing the cloud from the side, which changes the effective 1D interaction.`,
            p1: { name: "Modulation", sym: String.raw`\epsilon`, min: 0, max: 0.6, step: 0.05, value: 0.4 },
            p2: { name: "Frequency", sym: String.raw`\Omega`, min: 0.6, max: 3, step: 0.1, value: 2 },
            V: (xx) => (Math.abs(xx) < 30 ? 0 : WALL),
            gMod: (t, eps, Om) => 1 + eps * Math.sin(Om * t),
            kRange: 3,
            guess: (xx) => (Math.abs(xx) < 30 ? 1 : 0),
            packet: { x0: 0, k: 0, sigma: 3 },
        },
        trapbarrier: {
            label: "Trap with a barrier", absorbing: false, flux: false, scale: 3, bound: true,
            eq: String.raw`V(x) = \tfrac{1}{2}\omega^2 x^2 + U e^{-x^2/2w^2} + \epsilon\, \Theta(x)`,
            eqSub: String.raw`with \(\omega = 0.05\) and \(w = 0.5\); \(\epsilon\) raises the right half a little`,
            p1: { name: "Barrier", sym: "U", min: 0, max: 2, step: 0.01, value: 0.78 },
            p2: { name: "Step", sym: String.raw`\epsilon`, min: 0, max: 0.1, step: 0.001, value: 0 },
            V: (xx, U, eps) => 0.00125 * xx * xx + U * Math.exp(-(xx * xx) / 0.5) + (xx > 0 ? eps : 0),
            packet: { x0: -20, k: 0, sigma: 2 },
        },
        ring: {
            label: "Ring with a weak link", ring: true, absorbing: false, flux: false, scale: 0.4, bound: true, needsCondensate: true, kRange: 0.6,
            eq: String.raw`V(x) = U\, e^{-x^2/2w^2}`,
            eqSub: String.raw`on a ring of circumference \(L = 112\): the plot's two ends are the same point. \(w = 1.5\), \(\mu = g n\)`,
            p1: { name: "Weak link", sym: String.raw`U/\mu`, min: 0, max: 1.5, step: 0.05, value: 0.7 },
            V: (xx, uf) => uf * ringMu() * Math.exp(-(xx * xx) / 4.5),
            readout: ringReadout,
            guess: () => 1,
            packet: { x0: -20, k: 0.5, sigma: 3 },
        },
        ringstir: {
            label: "Stirred ring", ring: true, absorbing: false, flux: false, scale: 0.4, bound: true, needsCondensate: true, timeDep: true, kRange: 0.6,
            eq: String.raw`V(x, t) = U(t)\, e^{-(x - x_b(t))^2/2w^2}`,
            eqSub: String.raw`on a ring of circumference \(L = 112\). The weak link is swept around at speed \(v\) until \(t = 300\), then lowered to zero.`,
            p1: { name: "Speed", sym: "v/c", min: 0, max: 1, step: 0.05, value: 0.6 },
            p2: { name: "Height", sym: String.raw`U/\mu`, min: 0, max: 1.5, step: 0.05, value: 0.8 },
            V: (xx, vf, uf, t) => {
                const v = vf * Math.sqrt(ringMu()), ramp = 40, hold = 300, lower = 40;
                const tt = Math.min(t, hold);
                const pos = v * (tt < ramp ? (tt * tt) / (2 * ramp) : ramp / 2 + (tt - ramp));
                const amp = t < hold ? 1 : Math.max(0, 1 - (t - hold) / lower);
                return uf * ringMu() * amp * Math.exp(-(ringGap(xx - pos) ** 2) / 4.5);
            },
            readout: ringReadout,
            guess: () => 1,
            packet: { x0: -20, k: 0.5, sigma: 3 },
        },
        ramp: {
            label: "Opening trap", absorbing: false, flux: false, scale: 3, timeDep: true, bound: true,
            eq: String.raw`V(x, t) = \tfrac{1}{2}\,\omega(t)^2 x^2`,
            eqSub: String.raw`\(\omega\) goes from 0.5 to 0.05 in a time \(T\): along a designed shortcut (\(s = 1\)), a straight ramp (\(s = 0\)), or a mix`,
            p1: { name: "Ramp time", sym: "T", min: 5, max: 200, step: 5, value: 10 },
            p2: { name: "Shortcut", sym: "s", min: 0, max: 1, step: 0.25, value: 1 },
            V: (xx, T, sc, t) => 0.5 * rampOmega2(t, T, sc) * xx * xx,
            readout: () => `Cloud width: ${cloudWidth().toFixed(2)} (new ground state: ${(Math.sqrt(RAMP.w0 / RAMP.wf) / Math.sqrt(2 * RAMP.w0)).toFixed(2)})`,
            packet: { x0: 0, k: 0, sigma: 1 },
        },
        custom: {
            label: "Custom", absorbing: true, flux: true, scale: 5,
            eq: String.raw`V(x)`,
            eqSub: "drawn by hand on the plot",
        },
    };

    // A laser-speckle-like random potential: |E|^2 of a smoothed complex Gaussian
    // field, normalised to mean 1. Seeded and built on the default grid, then
    // interpolated, so the landscape is the same whatever dx is.
    const SPECKLE_N = Math.round((XMAX - XMIN) / DEFAULT_DX) - 1;
    const noise = (() => {
        let seed = 12345;
        const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
        const gaussian = () => {
            const u = rnd() || 1e-9, v = rnd();
            return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
        };
        const re = Float64Array.from({ length: SPECKLE_N }, gaussian);
        const im = Float64Array.from({ length: SPECKLE_N }, gaussian);
        return { re, im };
    })();
    const speckleCache = new Map();
    function speckleField(grain) {
        const key = `${grain.toFixed(2)}|${dx}`;
        if (speckleCache.has(key)) return speckleCache.get(key);
        const h = DEFAULT_DX, M = SPECKLE_N, m = Math.ceil((3 * grain) / h);
        const ker = Array.from({ length: 2 * m + 1 }, (_, i) => Math.exp(-(((i - m) * h) ** 2) / (2 * grain * grain)));
        const I = new Float64Array(M);
        for (let j = 0; j < M; j++) {
            let a = 0, b = 0;
            for (let i = -m; i <= m; i++) {
                const q = j + i;
                if (q < 0 || q >= M) continue;
                a += ker[i + m] * noise.re[q];
                b += ker[i + m] * noise.im[q];
            }
            I[j] = a * a + b * b;
        }
        const avg = I.reduce((sum, v) => sum + v) / M;
        const field = resample(I.map((v) => v / avg), (j) => XMIN + (j + 1) * h, M);
        speckleCache.set(key, field);
        // keep only the few most recent landscapes (~9 KB each)
        if (speckleCache.size > 4) speckleCache.delete(speckleCache.keys().next().value);
        return field;
    }

    // values f on points xOf(0..M-1) (evenly spaced), linearly interpolated onto x
    function resample(f, xOf, M) {
        const out = new Float64Array(N), x0 = xOf(0), h = xOf(1) - x0;
        for (let j = 0; j < N; j++) {
            const u = (x[j] - x0) / h, i = Math.max(0, Math.min(M - 2, Math.floor(u))), t = Math.max(0, Math.min(1, u - i));
            out[j] = f[i] * (1 - t) + f[i + 1] * t;
        }
        return out;
    }

    let presetKey = "barrier";
    let absorbing = true;
    let potVersion = 0;       // bumps whenever V changes, so cached eigenstates know to rebuild
    let curA = 0, curB = 0;   // the preset's slider values

    // Far above the plotted energy range a smooth potential acts as a wall anyway, but
    // letting it grow without bound (a quartic reaches ~10^4-10^6 at the grid's edge)
    // gives a packet placed there absurd energies the grid cannot resolve. Hard walls
    // (WALL) are left alone.
    function capPotential() {
        const c = Math.max(20, 4 * presets[presetKey].scale);
        for (let j = 0; j < N; j++) if (V[j] > c && V[j] < WALL) V[j] = c;
    }

    // a custom (painted) potential counts as bound when its edges are closed
    const isBound = () => (presetKey === "custom" ? !absorbing : !!presets[presetKey].bound);

    function buildPotential() {
        const p = presets[presetKey];
        if (presetKey !== "custom") {
            curA = p.p1 ? Number(ui.p1.value) : 0;
            curB = p.p2 ? Number(ui.p2.value) : 0;
            if (p.field) V.set(p.field(curA, curB));
            else for (let j = 0; j < N; j++) V[j] = p.V(x[j], curA, curB, tNow);
            capPotential();
        }
        potVersion++;
        buildAbsorber();
    }

    // back to the t = 0 potential (for a fresh launch), without invalidating caches
    function rewindTime() {
        tNow = 0;
        gScale = 1;
        stir.pos = stir.target = -25;
        stir.vel = 0;
        const p = presets[presetKey];
        if (!p || !p.timeDep) return;
        for (let j = 0; j < N; j++) V[j] = p.V(x[j], curA, curB, 0);
        capPotential();
        refresh();
    }

    function buildAbsorber() {
        const width = XMAX - CAP_START;
        for (let j = 0; j < N; j++) {
            const d = Math.abs(x[j]) - CAP_START;
            W[j] = absorbing && d > 0 ? CAP_STRENGTH * (d / width) ** 2 : 0;
        }
        refresh();
    }

    // Run fn with a temporary potential (filled by setV) and interaction g, then put
    // the real ones back. Used to prepare states in a different potential or at a
    // different g than the one they evolve in.
    function withHamiltonian(setV, gTemp, fn) {
        const saved = Float64Array.from(V), gSaved = g;
        if (setV) setV(V);
        g = gTemp;
        refresh();
        fn();
        V.set(saved);
        g = gSaved;
        refresh();
    }

    // A brief stretch of imaginary time (in the preset's untilted potential, if it
    // has one) damps the high-energy parts of a freshly placed packet.
    function settleIn(pr) {
        if (!pr || !pr.settle) return;
        normalize();
        const setV = pr.untilted ? (arr) => { for (let j = 0; j < N; j++) arr[j] = pr.untilted(x[j], curA); } : null;
        withHamiltonian(setV, g, () => {
            const wasImaginary = imaginary;
            imaginary = true;
            relaxScale = 5;
            refresh();
            for (let i = 0; i < Math.round(pr.settle / (5 * dt)); i++) step();
            relaxScale = 1;
            imaginary = wasImaginary;
        });
    }

    // ------------------------------------------------------------------
    // preparation: wavepacket / ground state (+ imprint) / left-well state
    // ------------------------------------------------------------------
    const sech = (u) => 1 / Math.cosh(u);
    // Thomas-Fermi profile sqrt(1 - u^2) with a slightly softened edge, so the
    // kink at the cloud's edge does not inject grid-scale wavenumbers.
    function tfAmp(u) {
        const n = 1 - u * u, d = 0.03;
        return Math.sqrt(0.5 * (n + Math.sqrt(n * n + d * d)));
    }

    let prep = {};
    const defaultPrep = {
        kind: "packet",          // packet | ground | leftwell
        shape: "gaussian",       // gaussian | pair | sech | sech2
        x0: -15, sigma: 3, k: 1.8, d: 16, phi: 0,
        prepIn: "same",          // same | any static preset | twotraps
        pa: 0.2, pb: 0,          // that potential's parameters
        sameG: true, g0: 10,
        imprint: "none",         // none | kick | soliton | bump
        kick: 1, wind: 2, v: 0, eps: 0.15, bw: 4,
        zoom: 1,                 // a demo can zoom the density axis
    };

    const g0Of = () => (prep.sameG ? g : prep.g0);
    // bright-soliton width for one soliton holding a fraction `share` of the atoms
    const brightWidth = (share) => (g < 0 ? 2 / (-g * share) : prep.sigma);

    function packetAmp(xx, p) {
        const ck = Math.cos(p.k * xx), sk = Math.sin(p.k * xx);
        if (p.shape === "gaussian" || p.shape === "sech") {
            const a = p.shape === "gaussian"
                ? Math.exp(-((xx - p.x0) ** 2) / (4 * p.sigma * p.sigma))
                : sech((xx - p.x0) / brightWidth(1));
            return [a * ck, a * sk];
        }
        // two packets heading for each other, relative phase phi
        const f = p.shape === "pair"
            ? (u) => Math.exp(-(u * u) / (4 * p.sigma * p.sigma))
            : (u) => sech(u / brightWidth(0.5));
        const a = f(xx - p.x0 + p.d / 2), b = f(xx - p.x0 - p.d / 2);
        const c2 = Math.cos(-p.k * xx + p.phi), s2 = Math.sin(-p.k * xx + p.phi);
        return [a * ck + b * c2, a * sk + b * s2];
    }

    // the potential the ground state is prepared in: any static preset with its own
    // parameters, or two separate harmonic traps
    const twoTraps = {
        p1: { name: "Frequency", sym: String.raw`\omega_0`, min: 0.05, max: 1, step: 0.01, value: 0.5 },
        p2: { name: "Separation", sym: "s", min: 4, max: 40, step: 1, value: 20 },
    };
    const prepSpec = () => (prep.prepIn === "twotraps" ? twoTraps : presets[prep.prepIn]);
    // a ground state only makes sense where the particle is held in one place (a
    // finite well always holds one bound state); in an open potential it would just
    // spread over the whole simulation box
    const holdsParticle = (key) => !!(presets[key].bound || presets[key].holds);
    const prepChoices = () => Object.keys(presets).filter((k) => k !== "custom" && !presets[k].timeDep && holdsParticle(k));

    function prepPotential(arr) {
        if (prep.prepIn === "twotraps") {
            for (let j = 0; j < N; j++) arr[j] = 0.5 * prep.pa * prep.pa * (Math.abs(x[j]) - prep.pb / 2) ** 2;
            return;
        }
        const pr = presets[prep.prepIn];
        if (pr.field) arr.set(pr.field(prep.pa, prep.pb));
        else for (let j = 0; j < N; j++) arr[j] = pr.V(x[j], prep.pa, prep.pb, 0);
        const c = Math.max(20, 4 * pr.scale);
        for (let j = 0; j < N; j++) if (arr[j] > c && arr[j] < WALL) arr[j] = c;
    }

    // lowest state of an open potential: spread over the whole grid
    const gridSine = (xx) => Math.sin((Math.PI * (xx - XMIN)) / (XMAX - XMIN));

    // a good first guess makes the relaxation converge quickly
    function groundGuess(xx, g0) {
        if (prep.prepIn === "same") {
            const pr = presets[presetKey];
            if (pr.guess) return pr.guess(xx, curA, curB);
            if (presetKey === "harmonic" && g0 > 0) return tfAmp(xx / tfRadius(curA, g0));
            return isBound() ? Math.exp(-(xx * xx) / 50) : gridSine(xx);
        }
        if (prep.prepIn === "twotraps") {
            const u = Math.abs(xx) - prep.pb / 2;
            return g0 > 0 ? tfAmp(u / tfRadius(prep.pa, g0)) : Math.exp(-0.5 * prep.pa * u * u);
        }
        const pr = presets[prep.prepIn];
        if (pr.guess) return pr.guess(xx, prep.pa, prep.pb);
        if (prep.prepIn === "harmonic") return g0 > 0 ? tfAmp(xx / tfRadius(prep.pa, g0)) : Math.exp(-0.5 * prep.pa * xx * xx);
        return pr.bound ? Math.exp(-(xx * xx) / 50) : gridSine(xx);
    }
    // radius of a Thomas-Fermi cloud (unit norm) in a harmonic trap of stiffness w
    function tfRadius(w, g0) {
        const mu = Math.pow((3 * g0 * w) / (4 * Math.SQRT2), 2 / 3);
        return Math.max(1, Math.sqrt(2 * mu) / w);
    }

    // Eigenstates are expensive (~50-100 ms of imaginary time), so keep the last
    // one per preparation and reuse it on every relaunch (~18 KB each).
    const gsCache = { key: null }, leftCache = { key: null };
    function cached(c, key) {
        if (c.key !== key) return false;
        pR.set(c.R);
        pI.set(c.I);
        return true;
    }
    function store(c, key) {
        c.key = key;
        c.R = Float64Array.from(pR);
        c.I = Float64Array.from(pI);
    }

    function preparedGround() {
        const g0 = g0Of();
        const key = [prep.prepIn, prep.prepIn === "same" ? potVersion : "", prep.pa, prep.pb, g0].join("|");
        if (cached(gsCache, key)) return;
        withHamiltonian(prep.prepIn === "same" ? null : prepPotential, g0, () => {
            for (let j = 0; j < N; j++) {
                pR[j] = groundGuess(x[j], g0);
                pI[j] = 0;
            }
            normalize();
            relax(4000);
        });
        store(gsCache, key);
    }

    // (psi_0 + psi_1)/sqrt(2) of the current potential without interactions, so that
    // g only changes how it evolves (Josephson oscillations vs self-trapping)
    function leftWellState() {
        if (cached(leftCache, String(potVersion))) return;
        withHamiltonian(null, 0, () => {
            for (let j = 0; j < N; j++) {
                pR[j] = Math.exp(-(x[j] ** 2) / 50);
                pI[j] = 0;
            }
            normalize();
            relax(4000);
            e0R.set(pR);
            e0I.set(pI);
            for (let j = 0; j < N; j++) {
                pR[j] = x[j] * Math.exp(-(x[j] ** 2) / 50);
                pI[j] = 0;
            }
            normalize();
            relax(4000, [e0R, e0I]);
            // pick the sign that puts the particle in the left well
            let left = 0, right = 0;
            for (let j = 0; j < N; j++) {
                const s2 = (e0R[j] + pR[j]) ** 2;
                if (x[j] < 0) left += s2; else right += s2;
            }
            const s = left >= right ? 1 : -1;
            for (let j = 0; j < N; j++) {
                pR[j] = e0R[j] + s * pR[j];
                pI[j] = e0I[j] + s * pI[j];
            }
        });
        store(leftCache, String(potVersion));
    }

    function imprint() {
        if (prep.imprint === "kick") {
            // on a ring the phase must close on itself: k = 2 pi w / L
            const k = periodic ? (2 * Math.PI * prep.wind) / RING_L : prep.kick;
            for (let j = 0; j < N; j++) {
                const c = Math.cos(k * x[j]), s = Math.sin(k * x[j]);
                const a = pR[j], b = pI[j];
                pR[j] = a * c - b * s;
                pI[j] = a * s + b * c;
            }
        } else if (prep.imprint === "bump") {
            for (let j = 0; j < N; j++) {
                const f = Math.sqrt(1 + prep.eps * Math.exp(-((x[j] - prep.x0) ** 2) / (2 * prep.bw * prep.bw)));
                pR[j] *= f;
                pI[j] *= f;
            }
        } else if (prep.imprint === "soliton") {
            // psi_0 * (gamma tanh(gamma (x - x0)/xi) + i v/c), xi the local healing length
            normalize();
            const jn = Math.max(0, Math.min(N - 1, Math.round((prep.x0 - XMIN) / dx - 1)));
            const n = pR[jn] * pR[jn] + pI[jn] * pI[jn];
            const xi = g > 0 && n > 0 ? Math.min(5, 1 / Math.sqrt(g * n)) : 1;
            const vc = Math.max(-0.95, Math.min(0.95, prep.v));
            const gam = Math.sqrt(1 - vc * vc);
            for (let j = 0; j < N; j++) {
                const t = gam * Math.tanh((gam * (x[j] - prep.x0)) / xi);
                const a = pR[j], b = pI[j];
                pR[j] = a * t - b * vc;
                pI[j] = a * vc + b * t;
            }
            if (vc === 0) {
                // a black soliton has a true node: let the cloud settle around it
                // (node pinned at x0) so the imprint does not set it breathing
                normalize();
                withHamiltonian(prep.prepIn === "same" ? null : prepPotential, g, () => relax(600, null, prep.x0));
            }
            // centre of mass at the trap centre: in a harmonic trap it would otherwise
            // slosh at omega regardless of the soliton (Kohn's theorem)
            if (presetKey === "harmonic") {
                const X = measure().x;
                rR.set(pR);
                rI.set(pI);
                for (let j = 0; j < N; j++) {
                    const f = (x[j] + X - x[0]) / dx, i = Math.floor(f), t = f - i;
                    const ok = i >= 0 && i < N - 1;
                    pR[j] = ok ? rR[i] * (1 - t) + rR[i + 1] * t : 0;
                    pI[j] = ok ? rI[i] * (1 - t) + rI[i + 1] * t : 0;
                }
            }
        }
    }

    function prepare() {
        if (prep.kind === "packet") {
            for (let j = 0; j < N; j++) [pR[j], pI[j]] = packetAmp(x[j], prep);
            settleIn(presets[presetKey]);
        } else if (prep.kind === "ground") {
            preparedGround();
            imprint();
        } else {
            leftWellState();
        }
        normalize();
    }

    // which state choices make sense right now
    const slowPrep = () => prep.kind !== "packet";              // built by relaxation
    const placeable = () => !presets[presetKey].stir && (prep.kind === "packet" || (prep.kind === "ground" && (prep.imprint === "soliton" || prep.imprint === "bump")));
    const classicalAllowed = () => prep.kind === "packet" && prep.shape === "gaussian" && g === 0 && !imaginary && !presets[presetKey].timeDep && !periodic;
    const fluxAllowed = () => !imaginary && presets[presetKey].flux && prep.kind === "packet" && (prep.shape === "gaussian" || prep.shape === "sech");

    // ------------------------------------------------------------------
    // wavepacket, classical particle, accuracy
    // ------------------------------------------------------------------
    let E = 0, rho0 = 1, time = 0;
    const classical = { x: 0, dir: 1, alive: false };
    let status = { level: "ok", text: "Accurate", tip: "" };
    let converged = false, lastE = Infinity;
    let driftDist = 0;   // accumulated position error from grid dispersion
    let finished = false; // everything has been absorbed at the edges
    let finishedAt = 0;   // for fading in the end message
    const FADE_MS = 600;

    function launch() {
        if (finished) {
            // the run had ended on its own: a new one starts it again
            finished = false;
            if (!playing) {
                playing = true;
                updatePlayButton();
                requestAnimationFrame(loop);
            }
        }
        rewindTime();
        prepare();
        rho0 = 0;
        for (let j = 0; j < N; j++) rho0 = Math.max(rho0, pR[j] * pR[j] + pI[j] * pI[j]);
        rho0 *= prep.zoom || 1;
        pPeak = 0;
        kBuf.P0 = null;
        resetHistory();
        absorbedLeft = absorbedRight = 0;
        time = 0;
        history.next = 0;
        snapshotHistory();
        E = measure().E;
        driftDist = 0;
        converged = false;
        lastE = Infinity;
        resetClassical(prep.x0, prep.k >= 0 ? 1 : -1);
        if (!classicalAllowed()) classical.alive = false;
        updateStatus();
        draw();
    }

    // A new potential is a new Hamiltonian: the state's energy changes with it.
    function potentialChanged() {
        converged = false;
        if (imaginary) return;
        E = measure().E;
        // the classical particle cannot sit where the new potential forbids it
        if (classical.alive && potentialAt(classical.x) >= E) classical.alive = false;
    }

    function resetClassical(x0, dir) {
        classical.x = x0;
        classical.dir = dir;
        classical.alive = potentialAt(x0) < E;
    }

    function stepClassical() {
        if (!classical.alive || !classicalAllowed()) return;
        // ride along the energy line: |p| = sqrt(2(E - V)), turning back where V > E
        const pMin = 0.08;
        const speed = Math.sqrt(Math.max(2 * (E - potentialAt(classical.x)), pMin * pMin));
        const xn = classical.x + classical.dir * speed * dt;
        if (absorbing && (xn < VMIN || xn > VMAX)) {
            classical.alive = false;
        } else if (potentialAt(xn) > E || xn <= XMIN + dx || xn >= XMAX - dx) {
            classical.dir *= -1;
        } else {
            classical.x = xn;
        }
    }

    // How far has the simulated wave drifted from the exact solution? On the
    // grid, a wave of wavenumber k travels at sin(kΔx)/Δx / (1 + (dt·ω/2)²)
    // instead of k (finite differences + Crank-Nicolson). That speed error, at
    // the state's rms wavenumber k = sqrt(2<T>), is accumulated over time, so it
    // catches long runs and quenches that create short wavelengths.
    function gridSpeedError(k) {
        if (k <= 0) return 0;
        const w = (1 - Math.cos(k * dx)) / (dx * dx);
        const v = Math.sin(k * dx) / dx / (1 + (0.5 * dt * w) ** 2);
        return Math.abs(k - v);
    }

    function updateStatus(m = measure(), elapsed = 0) {
        if (!Number.isFinite(m.norm) || !Number.isFinite(m.E)) {
            status = { level: "bad", text: "Numerical blow-up", tip: "The state is no longer finite; restart." };
            return false;
        }
        if (imaginary) {
            status = converged
                ? { level: "ok", text: "Ground state found", tip: "Energy has stopped changing in imaginary time." }
                : { level: "info", text: "Relaxing…", tip: "Imaginary time damps every excited state; the slowest to disappear is the first excited state." };
            return true;
        }
        // only judge what is on screen: once the packet has (mostly) left, hold the
        // last verdict instead of grading the tail that is being absorbed off-screen
        if (m.visible < 0.05) return true;
        const k = Math.sqrt(Math.max(0, 2 * m.T));
        driftDist += m.drift * elapsed;
        const kdx = k * dx;
        const pts = kdx > 0 ? Math.round((2 * Math.PI) / kdx) : "∞";
        const tip = `Grid spacing Δx = ${dx.toFixed(1)}, about ${pts} grid points per wavelength. ` +
            `On a grid, waves travel slightly too slowly; so far the moving parts are off by about ${driftDist.toFixed(2)} length units.`;
        if (kdx > 0.9 || driftDist > 4) {
            status = { level: "bad", text: "Unreliable", tip: tip + " Restart for an accurate run." };
        } else if (kdx > 0.5 || driftDist > 1) {
            status = { level: "warn", text: "Approximate", tip };
        } else {
            status = { level: "ok", text: "Accurate", tip };
        }
        return true;
    }

    // ------------------------------------------------------------------
    // UI
    // ------------------------------------------------------------------
    // slider specs for the initial-state controls
    const fields = {
        sigma: { name: "Width", sym: String.raw`\sigma`, min: 0.5, max: 6, step: 0.1 },
        k: { name: "Momentum", sym: "k", min: -K_MAX, max: K_MAX, step: 0.1 },
        d: { name: "Separation", sym: "d", min: 4, max: 40, step: 0.5 },
        phi: { name: "Phase", sym: String.raw`\phi`, min: 0, max: 3.14, step: 0.02 },
        g0: { name: "Prepared at", sym: "g_0", min: -6, max: 20, step: 0.5 },
        kick: { name: "Kick", sym: "k", min: -K_MAX, max: K_MAX, step: 0.1 },
        wind: { name: "Winding", sym: "w", min: -8, max: 8, step: 1 },
        v: { name: "Speed", sym: "v/c", min: -0.9, max: 0.9, step: 0.05 },
        eps: { name: "Height", sym: String.raw`\epsilon`, min: 0, max: 3, step: 0.05 },
        bw: { name: "Width", sym: String.raw`\sigma`, min: 0.5, max: 8, step: 0.1 },
    };
    const fieldRow = (key) => `<label class="tdse-field" data-ui="row_${key}"><span data-ui="name_${key}"></span><input type="range" data-ui="in_${key}"><output data-ui="out_${key}"></output></label>`;

    // the topic tabs and the guided demos live in tdse_demos.js (data only, meant
    // to be edited); captions may be wrapped over several lines
    const { topics, demos } = window.TDSE_DEMOS || { topics: [], demos: [] };
    demos.forEach((d) => (d.caption = d.caption.replace(/\s+/g, " ").trim()));

    root.classList.add("tdse");
    root.innerHTML = `
        <div class="tdse-demos">
            <div class="tdse-tabs" role="tablist" aria-label="Demo topics" data-ui="tabs"></div>
            <div class="tdse-chips" role="group" aria-label="Demos" data-ui="chips"></div>
            <p class="tdse-demo-caption" data-ui="demoCaption" aria-live="polite"></p>
        </div>
        <div class="tdse-stage">
            <canvas class="tdse-canvas" aria-label="Wavefunction plot. Click to place a wavepacket, drag to give it momentum."></canvas>
            <canvas class="tdse-strip tdse-strip-tall" data-ui="hcanvas" hidden aria-label="Space-time history of the density"></canvas>
            <canvas class="tdse-strip" data-ui="kcanvas" hidden aria-label="Momentum distribution"></canvas>
            <canvas class="tdse-strip" data-ui="phcanvas" hidden aria-label="Phase of the wavefunction"></canvas>
        </div>
        <div class="tdse-toolbar">
            <button type="button" class="tdse-btn" data-ui="play"></button>
            <button type="button" class="tdse-btn" data-ui="reset"><i class="bi bi-arrow-counterclockwise"></i><span>Restart</span></button>
            <div class="tdse-seg" role="group" aria-label="Time evolution">
                <button type="button" data-ui="modeReal" aria-pressed="true">Real time</button>
                <button type="button" data-ui="modeImag" aria-pressed="false" title="Evolve in imaginary time, renormalising after every step: relaxes to the ground state">Imaginary time</button>
            </div>
            <span class="tdse-status" data-ui="status" tabindex="0" role="status" aria-live="polite"><i class="bi" aria-hidden="true"></i><span></span></span>
            <label class="tdse-field tdse-speed"><span>Speed</span><input type="range" min="1" max="60" step="1" value="10" data-ui="speed"></label>
        </div>
        <details class="tdse-controls">
        <summary class="tdse-controls-toggle"><i class="bi bi-sliders" aria-hidden="true"></i><span>Tweak the simulation</span><small>Hamiltonian, initial state, display</small></summary>
        <div class="tdse-seg tdse-panel-tabs" role="group" aria-label="Settings">
            <button type="button" data-tab="ham" aria-pressed="true">Hamiltonian</button>
            <button type="button" data-tab="state" aria-pressed="false">Initial state</button>
            <button type="button" data-tab="display" aria-pressed="false">Display</button>
        </div>
        <div class="tdse-panel" data-tab="ham">
            <div class="tdse-group" data-group="ham" role="group" aria-labelledby="tdse-g-ham">
                <div class="tdse-group-title" id="tdse-g-ham">Hamiltonian</div>
                <p class="tdse-eq"><span class="math inline">\\(H = -\\tfrac{1}{2}\\,\\partial_x^2 + V(x) + g|\\psi|^2\\)</span></p>
                <p class="tdse-eq-sub">in units where <span class="math inline">\\(\\hbar = m = 1\\)</span></p>
                <p class="tdse-note tdse-imag-note" data-ui="imagNote" hidden>In imaginary time <span class="math inline">\\(t \\to -i\\tau\\)</span>, every eigenstate decays: <span class="math inline">\\(\\psi(\\tau) = \\sum_n c_n e^{-E_n \\tau} \\psi_n\\)</span>. The ground state decays slowest, so it is all that is left. The norm shrinks exponentially along the way, so <span class="math inline">\\(\\psi\\)</span> is renormalised after every step.</p>
                <label class="tdse-field"><span>Potential</span><select class="tdse-select" data-ui="preset"></select></label>
                <p class="tdse-eq tdse-eq-v" data-ui="potEq"></p>
                <p class="tdse-eq-sub" data-ui="potEqSub" hidden></p>
                <label class="tdse-field" data-ui="p1wrap"><span data-ui="p1name"></span><input type="range" data-ui="p1"><output data-ui="p1out"></output></label>
                <label class="tdse-field" data-ui="p2wrap"><span data-ui="p2name"></span><input type="range" data-ui="p2"><output data-ui="p2out"></output></label>
                <div class="tdse-field tdse-field-reset" title="Mean-field interaction between many bosons in the same state (Gross-Pitaevskii equation). g > 0 repels, g < 0 attracts."><label for="tdse-g-in">Interaction <span class="math inline">\\(g\\)</span></label><input type="range" id="tdse-g-in" min="-6" max="20" step="0.1" value="0" data-ui="g"><output data-ui="gout">0.0</output><button type="button" class="tdse-mini" data-ui="gZero" title="Back to g = 0: a single particle, no interactions" aria-label="Set g to 0"><i class="bi bi-x-lg" aria-hidden="true"></i></button></div>
                <div class="tdse-inline">
                    <button type="button" class="tdse-btn" data-ui="draw" aria-pressed="false" title="Paint the potential with the mouse; drag low to erase"><i class="bi bi-pencil"></i><span>Draw your own</span></button>
                    <label class="tdse-check"><input type="checkbox" data-ui="absorb">Absorbing edges</label>
                </div>
                <p class="tdse-note" data-ui="hamNote" hidden></p>
            </div>
            <div class="tdse-group" data-group="state" role="group" aria-labelledby="tdse-g-state">
                <div class="tdse-group-title" id="tdse-g-state">Initial state</div>
                <p class="tdse-eq" data-ui="stateEq"></p>
                <p class="tdse-eq-sub" data-ui="stateEqSub" hidden></p>
                <div class="tdse-step"><span class="tdse-step-n">1</span>Prepare</div>
                <div class="tdse-seg tdse-seg-full" role="group" aria-label="Start from" data-ui="kindSeg">
                    <button type="button" data-v="packet">Wavepacket</button>
                    <button type="button" data-v="ground">Ground state</button>
                    <button type="button" data-v="leftwell" title="Equal mix of the two lowest states: a particle in the left well">Left well</button>
                </div>
                <div class="tdse-seg tdse-seg-full" role="group" aria-label="Shape" data-ui="shapeSeg">
                    <button type="button" data-v="gaussian">Gaussian</button>
                    <button type="button" data-v="pair">Pair</button>
                    <button type="button" data-v="sech" title="Bright soliton (needs g < 0)">Bright</button>
                    <button type="button" data-v="sech2" title="Two bright solitons (need g < 0)">Bright pair</button>
                </div>
                ${["sigma", "k", "d", "phi"].map(fieldRow).join("")}
                <label class="tdse-field" data-ui="row_prepIn"><span>Of</span><select class="tdse-select" data-ui="prepIn"></select></label>
                <label class="tdse-field" data-ui="prepAwrap"><span data-ui="prepAname"></span><input type="range" data-ui="prepA"><output data-ui="prepAout"></output></label>
                <label class="tdse-field" data-ui="prepBwrap"><span data-ui="prepBname"></span><input type="range" data-ui="prepB"><output data-ui="prepBout"></output></label>
                <label class="tdse-check" data-ui="row_quench"><input type="checkbox" data-ui="quench">Prepare at a different <span class="math inline">\\(g\\)</span> (a quench)</label>
                ${fieldRow("g0")}
                <div class="tdse-step" data-ui="row_step2"><span class="tdse-step-n">2</span>Then</div>
                <div class="tdse-seg tdse-seg-full" role="group" aria-label="Then imprint" data-ui="imprintSeg">
                    <button type="button" data-v="none">Nothing</button>
                    <button type="button" data-v="kick" title="Give the whole cloud momentum k">Kick</button>
                    <button type="button" data-v="soliton" title="Imprint a dark soliton (needs g > 0)">Soliton</button>
                    <button type="button" data-v="bump" title="Add a density bump">Bump</button>
                </div>
                ${["kick", "wind", "v", "eps", "bw"].map(fieldRow).join("")}
                <p class="tdse-note tdse-warn" data-ui="stateWarn" hidden></p>
                <p class="tdse-note" data-ui="stateNote"></p>
            </div>
            <div class="tdse-group" data-group="display" role="group" aria-labelledby="tdse-g-display">
                <div class="tdse-group-title" id="tdse-g-display">Display</div>
                <label class="tdse-check"><input type="checkbox" checked data-ui="showRho"><i class="sw sw-rho"></i>Probability density</label>
                <label class="tdse-check"><input type="checkbox" data-ui="showRe"><i class="sw sw-re"></i>Real part of &psi;</label>
                <label class="tdse-check"><input type="checkbox" data-ui="showIm"><i class="sw sw-im"></i>Imaginary part of &psi;</label>
                <label class="tdse-check"><input type="checkbox" checked data-ui="showE"><i class="sw sw-e"></i>Average energy</label>
                <label class="tdse-check" title="The density over time, newest at the top: trajectories, oscillations and carpets in one picture"><input type="checkbox" data-ui="showHist"><i class="sw sw-hist"></i>Space–time history</label>
                <label class="tdse-check" title="How much of the wave moves with each momentum k, in a strip below the plot"><input type="checkbox" data-ui="showP"><i class="sw sw-p"></i>Momentum distribution</label>
                <label class="tdse-check" title="The phase of psi wherever there is density, in a strip below the plot"><input type="checkbox" data-ui="showPh"><i class="sw sw-ph"></i>Phase of &psi;</label>
                <label class="tdse-check" data-ui="clWrap"><input type="checkbox" checked data-ui="showCl"><i class="sw sw-cl"></i>Classical particle</label>
                <details class="tdse-adv">
                    <summary>Numerics</summary>
                    <label class="tdse-field"><span>Grid <span class="math inline">\\(\\Delta x\\)</span></span><select class="tdse-select" data-ui="numDx">
                        <option value="0.2">0.2 (coarse, fast)</option>
                        <option value="0.1" selected>0.1 (default)</option>
                        <option value="0.05">0.05 (fine)</option>
                    </select></label>
                    <label class="tdse-field"><span>Step <span class="math inline">\\(\\Delta t\\)</span></span><select class="tdse-select" data-ui="numDt">
                        <option value="0.02">0.02 (coarse, fast)</option>
                        <option value="0.01" selected>0.01 (default)</option>
                        <option value="0.005">0.005</option>
                        <option value="0.0025">0.0025 (fine)</option>
                    </select></label>
                    <p class="tdse-note" data-ui="numNote"></p>
                </details>
            </div>
        </div>
        </details>
        <div class="tdse-credit">
            <svg class="tdse-claude" viewBox="0 0 24 24" aria-hidden="true"><path d="m4.7144 15.9555 4.7174-2.6471.079-.2307-.079-.1275h-.2307l-.7893-.0486-2.6956-.0729-2.3375-.0971-2.2646-.1214-.5707-.1215-.5343-.7042.0546-.3522.4797-.3218.686.0608 1.5179.1032 2.2767.1578 1.6514.0972 2.4468.255h.3886l.0546-.1579-.1336-.0971-.1032-.0972L6.973 9.8356l-2.55-1.6879-1.3356-.9714-.7225-.4918-.3643-.4614-.1578-1.0078.6557-.7225.8803.0607.2246.0607.8925.686 1.9064 1.4754 2.4893 1.8336.3643.3035.1457-.1032.0182-.0728-.164-.2733-1.3539-2.4467-1.445-2.4893-.6435-1.032-.17-.6194c-.0607-.255-.1032-.4674-.1032-.7285L6.287.1335 6.6997 0l.9957.1336.419.3642.6192 1.4147 1.0018 2.2282 1.5543 3.0296.4553.8985.2429.8318.091.255h.1579v-.1457l.1275-1.706.2368-2.0947.2307-2.6957.0789-.7589.3764-.9107.7468-.4918.5828.2793.4797.686-.0668.4433-.2853 1.8517-.5586 2.9021-.3643 1.9429h.2125l.2429-.2429.9835-1.3053 1.6514-2.0643.7286-.8196.85-.9046.5464-.4311h1.0321l.759 1.1293-.34 1.1657-1.0625 1.3478-.8804 1.1414-1.2628 1.7-.7893 1.36.0729.1093.1882-.0183 2.8535-.607 1.5421-.2794 1.8396-.3157.8318.3886.091.3946-.3278.8075-1.967.4857-2.3072.4614-3.4364.8136-.0425.0304.0486.0607 1.5482.1457.6618.0364h1.621l3.0175.2247.7892.522.4736.6376-.079.4857-1.2142.6193-1.6393-.3886-3.825-.9107-1.3113-.3279h-.1822v.1093l1.0929 1.0686 2.0035 1.8092 2.5075 2.3314.1275.5768-.3218.4554-.34-.0486-2.2039-1.6575-.85-.7468-1.9246-1.621h-.1275v.17l.4432.6496 2.3436 3.5214.1214 1.0807-.17.3521-.6071.2125-.6679-.1214-1.3721-1.9246L14.38 17.959l-1.1414-1.9428-.1397.079-.674 7.2552-.3156.3703-.7286.2793-.6071-.4614-.3218-.7468.3218-1.4753.3886-1.9246.3157-1.53.2853-1.9004.17-.6314-.0121-.0425-.1397.0182-1.4328 1.9672-2.1796 2.9446-1.7243 1.8456-.4128.164-.7164-.3704.0667-.6618.4008-.5889 2.386-3.0357 1.4389-1.882.929-1.0868-.0062-.1579h-.0546l-6.3385 4.1164-1.1293.1457-.4857-.4554.0608-.7467.2307-.2429 1.9064-1.3114Z"/></svg>
            <span>Made with the assistance of Claude Opus 5.5</span>
        </div>`;

    const ui = {};
    root.querySelectorAll("[data-ui]").forEach((el) => (ui[el.dataset.ui] = el));

    // On narrow screens the three groups of settings are tabs (styles in theme/_widgets.scss)
    const panelTabs = [...root.querySelectorAll(".tdse-panel-tabs button")];
    panelTabs.forEach((tab) => tab.addEventListener("click", () => {
        root.querySelector(".tdse-panel").dataset.tab = tab.dataset.tab;
        panelTabs.forEach((t) => t.setAttribute("aria-pressed", String(t === tab)));
    }));
    const canvas = root.querySelector(".tdse-canvas");
    const ctx = canvas.getContext("2d");

    // Optional strips under the main plot, each behind a Display checkbox.
    const strips = [
        { toggle: ui.showHist, cv: ui.hcanvas, draw: (cx, sw, sh) => drawHistory(cx, sw, sh) },
        { toggle: ui.showP, cv: ui.kcanvas, draw: (cx, sw, sh) => drawMomentum(cx, sw, sh) },
        { toggle: ui.showPh, cv: ui.phcanvas, draw: (cx, sw, sh) => drawPhase(cx, sw, sh) },
    ];
    strips.forEach((st) => {
        st.cx = st.cv.getContext("2d");
        st.toggle.addEventListener("change", () => syncStrips());
    });
    function syncStrips() {
        strips.forEach((st) => (st.cv.hidden = !st.toggle.checked));
        resize();
        resetHistory();
    }

    // the potential menu, grouped so the long list stays easy to scan
    const presetGroups = [
        ["Scattering", ["barrier", "double", "step", "well", "free"]],
        ["Traps", ["harmonic", "box", "doublewell", "trapbarrier"]],
        ["Lattices and disorder", ["lattice", "speckle"]],
        ["Driven", ["pulse", "lens", "ramp", "bragg", "shaken", "sweep", "paul", "belt", "faraday", "obstacle", "stir"]],
        ["Ring", ["ring", "ringstir"]],
        ["Your own", ["custom"]],
    ];
    for (const [title, keys] of presetGroups) {
        const group = document.createElement("optgroup");
        group.label = title;
        for (const key of keys) group.append(new Option(presets[key].label, key));
        ui.preset.append(group);
    }
    // where a ground state can be prepared: here, any static potential, or two traps
    ui.prepIn.append(new Option("This potential", "same"));
    for (const key of prepChoices()) ui.prepIn.append(new Option(presets[key].label, key));
    ui.prepIn.append(new Option("Two harmonic traps", "twotraps"));

    let playing = !window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    let drawMode = false;

    const fmt = (v) => Number(v).toFixed(Math.abs(v) < 10 ? (Math.abs(v) < 1 && v !== 0 ? 2 : 1) : 0);
    // a slider's readout shows as many decimals as its step (0.001 -> 3, 0.05 -> 2, 10 -> 0)
    const fmtStep = (v, step) => Number(v).toFixed(Math.max(0, -Math.floor(Math.log10(Number(step)) + 1e-9)));

    function typeset(el) {
        if (window.MathJax && window.MathJax.typesetPromise) window.MathJax.typesetPromise([el]);
    }

    // initial-state slider labels (static, typeset once with the page)
    for (const [key, spec] of Object.entries(fields)) {
        ui[`name_${key}`].innerHTML = `${spec.name} <span class="math inline">\\(${spec.sym}\\)</span>`;
        Object.assign(ui[`in_${key}`], { min: spec.min, max: spec.max, step: spec.step });
    }

    function configureSlider(which, spec, value = spec && spec.value) {
        const wrap = ui[which + "wrap"];
        if (!spec) {
            wrap.hidden = true;
            return;
        }
        wrap.hidden = false;
        const label = spec.sym ? `${spec.name} <span class="math inline">\\(${spec.sym}\\)</span>` : spec.name;
        if (ui[which + "name"].dataset.label !== label) {
            ui[which + "name"].dataset.label = label;
            ui[which + "name"].innerHTML = label;
            typeset(ui[which + "name"]);
        }
        Object.assign(ui[which], { min: spec.min, max: spec.max, step: spec.step });
        ui[which].value = value;
        ui[which + "out"].textContent = fmtStep(value, spec.step);
    }

    function setInteraction(value) {
        g = value;
        ui.g.value = value;
        ui.gout.textContent = fmt(value);
        ui.gZero.disabled = value === 0;
        refresh();
        potentialChanged();
    }

    // ---- the initial-state panel follows `prep`, and `prep` is kept consistent ----
    function visibleFields() {
        if (prep.kind === "packet") {
            const two = prep.shape === "pair" || prep.shape === "sech2";
            return [prep.shape === "gaussian" || prep.shape === "pair" ? "sigma" : null, "k", two ? "d" : null, two ? "phi" : null];
        }
        if (prep.kind === "ground") {
            return [
                prep.sameG ? null : "g0",
                prep.imprint === "kick" ? (periodic ? "wind" : "kick") : null,
                prep.imprint === "soliton" ? "v" : null,
                prep.imprint === "bump" ? "eps" : null,
                prep.imprint === "bump" ? "bw" : null,
            ];
        }
        return [];
    }

    function stateEquation() {
        if (prep.kind === "leftwell") return String.raw`\psi(x, 0) = \big(\psi_0(x) + \psi_1(x)\big)/\sqrt{2}`;
        if (prep.kind === "ground") {
            return {
                none: String.raw`\psi(x, 0) = \psi_0(x)`,
                kick: periodic ? String.raw`\psi(x, 0) = \psi_0(x)\, e^{2\pi i w x / L}` : String.raw`\psi(x, 0) = \psi_0(x)\, e^{i k x}`,
                soliton: String.raw`\psi(x, 0) \propto \psi_0(x) \Big[ i\tfrac{v}{c} + \gamma \tanh\!\tfrac{\gamma (x - x_0)}{\xi} \Big]`,
                bump: String.raw`\psi(x, 0) \propto \psi_0(x)\, \sqrt{1 + \epsilon\, e^{-(x - x_0)^2/2\sigma^2}}`,
            }[prep.imprint];
        }
        const two = String.raw`\psi \propto f(x - x_1)\, e^{i k x} + e^{i\phi} f(x - x_2)\, e^{-i k x}`;
        return {
            gaussian: String.raw`\psi(x, 0) \propto e^{-(x - x_0)^2 / 4\sigma^2}\, e^{i k x}`,
            pair: two,
            sech: String.raw`\psi(x, 0) \propto \mathrm{sech}\big((x - x_0)/\xi\big)\, e^{i k x}`,
            sech2: two,
        }[prep.shape];
    }

    // the definitions an equation above needs, in small print
    function stateEquationSub() {
        if (prep.kind !== "packet") return "";
        return {
            gaussian: "",
            pair: String.raw`with \(f(u) = e^{-u^2/4\sigma^2}\), \(x_{1,2} = x_0 \mp d/2\)`,
            sech: String.raw`with width \(\xi = 2/|g|\)`,
            sech2: String.raw`with \(f(u) = \mathrm{sech}(u/\xi)\), \(\xi = 4/|g|\), \(x_{1,2} = x_0 \mp d/2\)`,
        }[prep.shape];
    }

    function stateNote() {
        if (prep.kind === "leftwell") return String.raw`The two lowest states \(\psi_0, \psi_1\) of the potential (without interactions), added so that they cancel in the right well. It tunnels back and forth with period \(2\pi/(E_1 - E_0)\).`;
        if (prep.kind === "ground") {
            const where = prep.prepIn === "same" ? "this potential"
                : prep.prepIn === "twotraps" ? String.raw`two harmonic traps \(\tfrac{1}{2}\omega_0^2 (|x| - s/2)^2\)`
                : `the ${presets[prep.prepIn].label.toLowerCase()} potential`;
            const gTxt = prep.sameG ? "the same interaction g" : String.raw`interaction \(g_0\)`;
            const extra = {
                none: "",
                kick: periodic
                    ? " The phase winds w times around the ring: a current whose circulation can only take whole-number values."
                    : " The kick sets the whole cloud moving.",
                soliton: String.raw` The imprint is a soliton with \(\gamma = \sqrt{1 - v^2/c^2}\), sound speed \(c = \sqrt{g\,n}\) and healing length \(\xi = 1/\sqrt{g\,n}\) at \(x_0\); at \(v = 0\) it is black and stands still in a flat condensate.`,
                bump: " A small, wide bump splits into two sound pulses; a large one breaks into shock waves.",
            }[prep.imprint];
            return String.raw`\(\psi_0\) is the ground state of ${where} with ${gTxt}, then evolved with the current potential and \(g\).` + extra;
        }
        const hint = "Click the plot to set where it starts, or drag to throw it.";
        return {
            gaussian: hint,
            pair: "Two packets heading for each other; where they overlap they interfere. " + hint,
            sech: "A bright soliton holds itself together through attractive interactions. " + hint,
            sech2: String.raw`Two solitons with relative phase \(\phi\): in phase they pass through each other, out of phase (\(\phi = \pi\)) they repel. ` + hint,
        }[prep.shape];
    }

    // Keep the choices consistent with the Hamiltonian; returns warnings to show.
    function validate() {
        const warn = [];
        const bound = isBound();
        // options that make no sense for an open potential
        ui.kindSeg.querySelector('[data-v="leftwell"]').disabled = !bound;
        if (prep.kind === "leftwell" && !bound) {
            prep.kind = "packet";
            warn.push("The left-well state needs a potential that traps the particle, so a wavepacket is used instead.");
        }
        const sameOK = bound || (presetKey !== "custom" && !!presets[presetKey].holds);
        ui.prepIn.querySelector('[value="same"]').disabled = !sameOK;
        if (prep.kind === "ground" && prep.prepIn === "same" && !sameOK) {
            prep.prepIn = "harmonic";
            prep.pa = presets.harmonic.p1.value;
            warn.push("This potential does not hold the particle anywhere, so the ground state is prepared in a harmonic trap instead.");
        }
        if (prep.kind === "packet" && (prep.shape === "sech" || prep.shape === "sech2") && g >= 0)
            warn.push("Bright solitons only hold together with attractive interactions (g < 0). With this g the packet will spread.");
        if (prep.kind === "ground" && prep.imprint === "soliton" && g <= 0)
            warn.push("Dark solitons need repulsive interactions (g > 0).");
        // classical particle: one Gaussian packet without interactions, in real time
        const cl = classicalAllowed();
        ui.showCl.disabled = !cl;
        ui.clWrap.title = cl ? "" : "Only shown for a single Gaussian packet without interactions, in real time.";
        if (!cl) classical.alive = false;
        // the Hamiltonian's own caveat
        const hamWarn = presets[presetKey].needsCondensate && g <= 0
            ? "The obstacle's speed and height are measured in units of the condensate's sound speed and chemical potential, so they need g > 0."
            : "";
        ui.hamNote.textContent = hamWarn;
        ui.hamNote.hidden = !hamWarn;
        ui.hamNote.classList.toggle("tdse-warn", !!hamWarn);
        return warn;
    }

    // bring every initial-state control in line with `prep` (and validate it)
    function syncPrepUI(extraWarn = []) {
        const warn = validate().concat(extraWarn);
        const press = (seg, v) => seg.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.v === v)));
        press(ui.kindSeg, prep.kind);
        press(ui.shapeSeg, prep.shape);
        press(ui.imprintSeg, prep.imprint);
        ui.prepIn.value = prep.prepIn;
        ui.quench.checked = !prep.sameG;
        const ground = prep.kind === "ground";
        ui.shapeSeg.hidden = prep.kind !== "packet";
        ui.row_prepIn.hidden = ui.row_quench.hidden = ui.row_step2.hidden = ui.imprintSeg.hidden = !ground;
        const spec = ground && prep.prepIn !== "same" ? prepSpec() : null;
        configureSlider("prepA", spec && spec.p1, prep.pa);
        configureSlider("prepB", spec && spec.p2, prep.pb);
        const show = new Set(visibleFields().filter(Boolean));
        for (const key of Object.keys(fields)) {
            ui[`row_${key}`].hidden = !show.has(key);
            ui[`in_${key}`].value = prep[key];
            ui[`out_${key}`].textContent = fmtStep(prep[key], fields[key].step);
        }
        ui.stateEq.innerHTML = `<span class="math inline">\\(${stateEquation()}\\)</span>`;
        const sub = stateEquationSub();
        ui.stateEqSub.innerHTML = sub;
        ui.stateEqSub.hidden = !sub;
        typeset(ui.stateEqSub);
        ui.stateNote.innerHTML = stateNote();
        ui.stateWarn.textContent = warn.join(" ");
        ui.stateWarn.hidden = !warn.length;
        typeset(ui.stateEq);
        typeset(ui.stateNote);
    }

    // Numerics settings. A new dx rebuilds the grid (a painted potential is carried
    // over by interpolation); cached eigenstates and landscapes are dropped.
    function setNumerics(newDx, newDt) {
        if (newDx !== dx) {
            const oldX = x, oldV = V, oldN = N;
            allocateGrid(newDx);
            speckleCache.clear();
            gsCache.key = leftCache.key = null;
            if (presetKey === "custom") V.set(resample(oldV, (j) => oldX[j], oldN));
        }
        dt = newDt;
        buildPotential();
        launch();
        updateNumericsNote();
    }

    function updateNumericsNote() {
        const cost = ((DEFAULT_DX / dx) * (DEFAULT_DT / dt)).toFixed(2).replace(/\.?0+$/, "");
        ui.numNote.innerHTML = `${N} grid points on \\([${XMIN}, ${XMAX}]\\), absorbing beyond \\(|x| = ${CAP_START}\\). ` +
            `The simulation covers the same time per frame whatever \\(\\Delta t\\) is, so each frame costs ${cost}× the default. ` +
            "The accuracy indicator takes both into account.";
        typeset(ui.numNote);
    }

    function selectPreset(key, relaunch = true) {
        presetKey = key;
        ui.preset.value = key;
        const p = presets[key];
        if (key !== "custom") {
            configureSlider("p1", p.p1);
            configureSlider("p2", p.p2);
        } else {
            ui.p1wrap.hidden = ui.p2wrap.hidden = true;
        }
        ui.potEq.innerHTML = `<span class="math display">\\[${p.eq}\\]</span>`;
        ui.potEqSub.innerHTML = p.eqSub || "";
        ui.potEqSub.hidden = !p.eqSub;
        typeset(ui.potEq);
        typeset(ui.potEqSub);
        absorbing = p.absorbing;
        ui.absorb.checked = absorbing;
        ui.absorb.disabled = !!p.ring;   // a ring has no edges
        if (!!p.ring !== periodic) {
            periodic = !!p.ring;
            if (periodic) setView(XMIN, XMAX); else setView(-40, 40);
        }
        buildPotential();
        if (relaunch) {
            // a single Gaussian takes the potential's suggested starting point
            if (prep.kind === "packet" && prep.shape === "gaussian" && p.packet) Object.assign(prep, p.packet);
            syncPrepUI();
            launch();
        }
    }

    function setMode(imag) {
        if (imag === imaginary) return;
        if (finished) launch();   // nothing left to relax: start again
        imaginary = imag;
        ui.modeReal.setAttribute("aria-pressed", String(!imag));
        ui.modeImag.setAttribute("aria-pressed", String(imag));
        root.classList.toggle("tdse-imaginary", imag);
        ui.imagNote.hidden = !imag;
        refresh();
        if (imag) {
            normalize();
            converged = false;
            lastE = Infinity;
        } else {
            // carry on in real time from wherever imaginary time got to
            const m = measure();
            E = m.E;
            driftDist = 0;
            absorbedLeft = absorbedRight = 0;
            time = 0;
            history.next = 0;
            resetClassical(m.x, 1);
        }
        validate();
        updateStatus();
        draw();
    }

    // ---- wiring: Hamiltonian ----
    ui.preset.addEventListener("change", () => selectPreset(ui.preset.value));
    ["p1", "p2"].forEach((which) =>
        ui[which].addEventListener("input", () => {
            ui[which + "out"].textContent = fmtStep(ui[which].value, ui[which].step);
            buildPotential();
            potentialChanged();
            draw();
        })
    );
    ui.gZero.addEventListener("click", () => {
        setInteraction(0);
        syncPrepUI();
        if (prep.kind === "packet" && (prep.shape === "sech" || prep.shape === "sech2")) launch();
        draw();
    });
    ui.g.addEventListener("input", () => {
        setInteraction(Number(ui.g.value));
        syncPrepUI();
        // a bright soliton's width is set by g: re-imprint it so it stays exact
        if (prep.kind === "packet" && (prep.shape === "sech" || prep.shape === "sech2")) launch();
        draw();
    });
    ui.absorb.addEventListener("change", () => {
        absorbing = ui.absorb.checked;
        buildAbsorber();
        syncPrepUI();   // closing the edges can make a custom potential bound
    });

    // ---- wiring: initial state ----
    // segmented choices behave like the selects they replace: a change event bubbles up
    function onSeg(seg, fn) {
        seg.addEventListener("click", (e) => {
            const b = e.target.closest("button");
            if (!b || b.disabled || b.getAttribute("aria-pressed") === "true") return;
            fn(b.dataset.v);
            syncPrepUI();
            launch();
            seg.dispatchEvent(new Event("change", { bubbles: true }));
        });
    }
    onSeg(ui.kindSeg, (v) => {
        prep.kind = v;
        if (v === "ground") prep.prepIn = "same";
    });
    onSeg(ui.shapeSeg, (v) => {
        prep.shape = v;
        if (v === "pair" || v === "sech2") prep.x0 = 0;
    });
    onSeg(ui.imprintSeg, (v) => {
        prep.imprint = v;
        if (v === "soliton" || v === "bump") prep.x0 = v === "soliton" ? -4 : 0;
    });
    ui.prepIn.addEventListener("change", () => {
        prep.prepIn = ui.prepIn.value;
        const spec = prepSpec();
        if (spec) {
            prep.pa = spec.p1 ? spec.p1.value : 0;
            prep.pb = spec.p2 ? spec.p2.value : 0;
        }
        syncPrepUI();
        launch();
    });
    ["prepA", "prepB"].forEach((which) => {
        ui[which].addEventListener("input", () => {
            prep[which === "prepA" ? "pa" : "pb"] = Number(ui[which].value);
            ui[which + "out"].textContent = fmtStep(ui[which].value, ui[which].step);
        });
        ui[which].addEventListener("change", () => launch());
    });
    ui.quench.addEventListener("change", () => {
        prep.sameG = !ui.quench.checked;
        if (!prep.sameG) prep.g0 = g;   // start from "no quench"
        syncPrepUI();
        launch();
    });
    // sliders: packets relaunch while dragging, states built by relaxation on release
    for (const key of Object.keys(fields)) {
        const input = ui[`in_${key}`];
        input.addEventListener("input", () => {
            prep[key] = Number(input.value);
            ui[`out_${key}`].textContent = fmtStep(prep[key], fields[key].step);
            if (!slowPrep()) launch();
        });
        input.addEventListener("change", () => {
            if (slowPrep()) launch();
        });
    }

    // ---- wiring: playback & display ----
    function updatePlayButton() {
        ui.play.innerHTML = playing
            ? '<i class="bi bi-pause-fill"></i><span>Pause</span>'
            : '<i class="bi bi-play-fill"></i><span>Play</span>';
    }

    ui.play.addEventListener("click", () => {
        playing = !playing;
        if (playing && finished) launch();   // nothing left to evolve: start over
        updatePlayButton();
        if (playing) loop();
    });
    ui.reset.addEventListener("click", () => launch());
    ui.modeReal.addEventListener("click", () => setMode(false));
    ui.modeImag.addEventListener("click", () => setMode(true));

    ui.draw.addEventListener("click", () => {
        drawMode = !drawMode;
        ui.draw.setAttribute("aria-pressed", String(drawMode));
        root.classList.toggle("tdse-drawing", drawMode);
        ui.draw.querySelector("span").textContent = drawMode ? "Done drawing" : "Draw your own";
        draw();
    });

    ["showRho", "showRe", "showIm", "showE", "showCl"].forEach((k) => ui[k].addEventListener("change", () => {
        layoutAxis();
        draw();
    }));
    ui.numDx.addEventListener("change", () => setNumerics(Number(ui.numDx.value), dt));
    ui.numDt.addEventListener("change", () => setNumerics(dx, Number(ui.numDt.value)));

    // Re/Im swing below zero and need room under the axis; otherwise use the full height
    function layoutAxis() {
        const signed = ui.showRe.checked || ui.showIm.checked;
        y0 = Math.round(signed ? h * 0.7 : h - 14) + 0.5;
    }

    // ------------------------------------------------------------------
    // demos (defined in tdse_demos.js): snap every setting to a scenario
    // ------------------------------------------------------------------

    const displayDefaults = { showRho: true, showRe: false, showIm: false, showE: true, showCl: true, showP: false, showPh: false, showHist: false };
    let activeDemo = -1, activeTopic = "tunnel";

    function showTopic(topic) {
        activeTopic = topic;
        ui.tabs.querySelectorAll("button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.topic === topic)));
        ui.chips.replaceChildren(...demos.map((d, i) => [d, i]).filter(([d]) => d.topic === topic).map(([d, i]) => {
            const b = document.createElement("button");
            b.type = "button";
            b.className = "tdse-chip";
            b.textContent = d.label;
            b.setAttribute("aria-pressed", String(i === activeDemo));
            b.addEventListener("click", () => applyDemo(i));
            return b;
        }));
    }

    // Tweaking a slider keeps the demo (most captions ask for exactly that) and only
    // notes that the settings changed; a new potential or a new kind of state ends it.
    let demoModified = false;
    function setActiveDemo(i) {
        activeDemo = i;
        demoModified = false;
        if (i >= 0 && demos[i].topic !== activeTopic) showTopic(demos[i].topic);
        ui.chips.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(b.textContent === (demos[i] || {}).label)));
        ui.demoCaption.textContent = i < 0
            ? "You have changed the setup. Pick a demo to start one of the guided examples again."
            : demos[i].caption;
        if (i >= 0 && demos[i].refs) {
            // the experiments (or, for the oldest, the papers) this demo follows
            const line = document.createElement("span");
            line.className = "tdse-demo-ref";
            line.append(demos[i].refs.length > 1 ? "Experiments: " : "Experiment: ");
            demos[i].refs.forEach(([text, href], n) => {
                if (n) line.append(", ");
                const a = document.createElement("a");
                a.href = href;
                a.textContent = text;
                a.target = "_blank";
                a.rel = "noopener";
                line.append(a);
            });
            ui.demoCaption.append(line);
        }
        ui.demoCaption.classList.toggle("tdse-demo-caption-custom", i < 0);
    }
    function markModified() {
        if (activeDemo < 0 || demoModified) return;
        demoModified = true;
        const note = document.createElement("span");
        note.className = "tdse-demo-mod";
        note.textContent = "Settings changed. Click the demo again to reset it.";
        ui.demoCaption.append(note);
    }

    for (const [topic, title] of topics) {
        const b = document.createElement("button");
        b.type = "button";
        b.setAttribute("role", "tab");
        b.dataset.topic = topic;
        b.textContent = title;
        b.addEventListener("click", () => showTopic(topic));
        ui.tabs.append(b);
    }

    function setSlider(which, v) {
        ui[which].value = v;
        ui[which + "out"].textContent = fmtStep(v, ui[which].step);
    }

    function applyDemo(i, autoplay = true) {
        const d = demos[i];
        setMode(false);
        setInteraction(d.g || 0);
        selectPreset(d.preset, false);
        if (d.p1 !== undefined) setSlider("p1", d.p1);
        if (d.p2 !== undefined) setSlider("p2", d.p2);
        buildPotential();
        prep = { ...defaultPrep, ...d.prep };
        for (const [k, v] of Object.entries({ ...displayDefaults, ...(d.show || {}) })) ui[k].checked = v;
        syncStrips();
        layoutAxis();
        ui.speed.value = d.speed || 10;
        history.dt = d.hist || 0.5;
        history.mode = d.histMode || "density";
        syncPrepUI();
        launch();
        if (d.imaginary) setMode(true);
        setActiveDemo(i);
        if (autoplay) {
            playing = true;
            updatePlayButton();
            loop();
        }
    }

    // any manual change means the demo's description no longer applies
    const endsDemo = (e) => e.target === ui.preset || ui.kindSeg.contains(e.target);
    ["input", "change"].forEach((type) => root.querySelector(".tdse-panel").addEventListener(type, (e) => (endsDemo(e) ? setActiveDemo(-1) : markModified())));

    // ------------------------------------------------------------------
    // accuracy badge
    // ------------------------------------------------------------------
    const statusIcons = {
        ok: "bi-check-circle-fill",
        info: "bi-hourglass-split",
        warn: "bi-exclamation-triangle-fill",
        bad: "bi-x-octagon-fill",
    };
    let shownLevel = null;

    function renderStatus() {
        const label = ui.status.querySelector("span");
        if (label.textContent !== status.text) label.textContent = status.text;
        ui.status.title = status.tip;
        if (status.level === shownLevel) return;
        // level changed: swap icon and pulse to draw the eye
        ui.status.dataset.level = status.level;
        ui.status.querySelector("i").className = "bi " + statusIcons[status.level];
        if (shownLevel !== null) {
            ui.status.classList.remove("tdse-pulse");
            void ui.status.offsetWidth;   // restart the animation
            ui.status.classList.add("tdse-pulse");
        }
        shownLevel = status.level;
    }

    // ------------------------------------------------------------------
    // drawing
    // ------------------------------------------------------------------
    // The palette comes from CSS custom properties on the widget (styles.scss), so
    // it follows the light/dark theme; it is read again whenever the theme flips.
    const css = getComputedStyle(root);
    const color = (name, fallback) => css.getPropertyValue(name).trim() || fallback;
    const rgbOf = (hex, fallback) => {
        const m = /^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(hex);
        return m ? m.slice(1).map((h) => parseInt(h, 16)) : fallback;
    };
    const C = {};
    let rhoRGB, dipRGB;                     // the same two colours as RGB, for the space-time strip
    function readColors() {
        Object.assign(C, {
            rho: color("--tdse-rho", "#5a2a9c"),
            rhoFill: color("--tdse-rho-fill", "rgba(90, 42, 156, 0.18)"),
            re: color("--tdse-re", "#2a2338"),
            im: color("--tdse-im", "#c2410c"),
            pot: color("--tdse-pot", "rgba(42, 35, 56, 0.09)"),
            potLine: color("--tdse-pot-line", "rgba(42, 35, 56, 0.45)"),
            axis: color("--tdse-axis", "rgba(42, 35, 56, 0.25)"),
            text: color("--tdse-text", "#75707e"),
            backing: color("--tdse-backing", "rgba(241, 235, 247, 0.92)"),
        });
        rhoRGB = rgbOf(C.rho, [90, 42, 156]);
        dipRGB = rgbOf(C.im, [194, 65, 12]);
    }
    readColors();

    let w = 0, h = 0, y0 = 0, eTop = 1;
    const top = 24;

    function fitCanvas(cv, cx) {
        const dpr = window.devicePixelRatio || 1;
        cv.width = Math.round(cv.clientWidth * dpr);
        cv.height = Math.round(cv.clientHeight * dpr);
        cx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function resize() {
        w = canvas.clientWidth;
        h = canvas.clientHeight;
        fitCanvas(canvas, ctx);
        strips.forEach((st) => !st.cv.hidden && fitCanvas(st.cv, st.cx));
        layoutAxis();
        draw();
    }

    const px = (xx) => ((xx - VMIN) / (VMAX - VMIN)) * w;
    const unpx = (p) => VMIN + (p / w) * (VMAX - VMIN);
    const yE = (e) => y0 - (e / eTop) * (y0 - top);
    const unyE = (y) => ((y0 - y) / (y0 - top)) * eTop;
    const yRho = (r) => y0 - (r / rho0) * (y0 - top) * 0.7;
    const yAmp = (a) => y0 - (a / Math.sqrt(rho0)) * (h - y0 - 6);
    const clampY = (y) => Math.max(2, Math.min(h - 2, y));

    function pathOf(fn) {
        ctx.beginPath();
        for (let j = j0; j <= j1; j++) {
            const X = px(x[j]), Y = clampY(fn(j));
            if (j === j0) ctx.moveTo(X, Y); else ctx.lineTo(X, Y);
        }
    }

    function draw() {
        if (!w || !onScreen()) return;
        renderStatus();
        eTop = presets[presetKey].scale;
        ctx.clearRect(0, 0, w, h);

        // potential
        ctx.beginPath();
        ctx.moveTo(0, y0);
        for (let j = j0; j <= j1; j++) ctx.lineTo(px(x[j]), Math.max(top, yE(V[j])));
        ctx.lineTo(w, y0);
        ctx.closePath();
        ctx.fillStyle = C.pot;
        ctx.fill();
        pathOf((j) => Math.max(top, yE(V[j])));
        ctx.strokeStyle = C.potLine;
        ctx.lineWidth = 1;
        ctx.stroke();

        // axis
        ctx.beginPath();
        ctx.moveTo(0, y0);
        ctx.lineTo(w, y0);
        ctx.strokeStyle = C.axis;
        ctx.stroke();

        // energy line
        const yEnergy = Math.max(top, yE(E));
        if (ui.showE.checked && !finished) {
            ctx.save();
            ctx.setLineDash([4, 4]);
            ctx.beginPath();
            ctx.moveTo(0, yEnergy);
            ctx.lineTo(w, yEnergy);
            ctx.strokeStyle = C.im;
            ctx.globalAlpha = 0.6;
            ctx.stroke();
            ctx.restore();
        }

        // |psi|^2
        if (ui.showRho.checked) {
            pathOf((j) => yRho(pR[j] * pR[j] + pI[j] * pI[j]));
            ctx.lineTo(w, y0);
            ctx.lineTo(0, y0);
            ctx.closePath();
            ctx.fillStyle = C.rhoFill;
            ctx.fill();
            pathOf((j) => yRho(pR[j] * pR[j] + pI[j] * pI[j]));
            ctx.strokeStyle = C.rho;
            ctx.lineWidth = 1.6;
            ctx.stroke();
        }
        if (ui.showRe.checked) {
            pathOf((j) => yAmp(pR[j]));
            ctx.strokeStyle = C.re;
            ctx.lineWidth = 1.1;
            ctx.stroke();
        }
        if (ui.showIm.checked) {
            pathOf((j) => yAmp(pI[j]));
            ctx.strokeStyle = C.im;
            ctx.lineWidth = 1.1;
            ctx.stroke();
        }

        // classical particle rides on the energy line
        if (ui.showCl.checked && classical.alive && classicalAllowed()) {
            ctx.beginPath();
            ctx.arc(px(classical.x), yEnergy, 5, 0, Math.PI * 2);
            ctx.fillStyle = "#ffffff";
            ctx.fill();
            ctx.lineWidth = 2;
            ctx.strokeStyle = C.im;
            ctx.stroke();
        }

        // readouts
        ctx.font = "13px 'Source Sans 3', sans-serif";
        ctx.fillStyle = C.text;
        ctx.textBaseline = "top";
        ctx.textAlign = "right";
        ctx.fillText(`${imaginary ? "τ" : "t"} = ${time.toFixed(1)}`, w - 8, 6);
        ctx.textAlign = "left";
        if (imaginary) {
            ctx.fillText(`Energy: ${E.toFixed(4)}`, 8, 6);
        } else if (presets[presetKey].readout) {
            ctx.fillText(presets[presetKey].readout(), 8, 6);
        } else if (fluxAllowed()) {
            let left = absorbedLeft, right = absorbedRight;
            for (let j = 0; j < N; j++) {
                const r = (pR[j] * pR[j] + pI[j] * pI[j]) * dx;
                if (x[j] < 0) left += r; else right += r;
            }
            const tot = left + right || 1;
            const trans = (prep.x0 < 0 ? right : left) / tot;
            ctx.fillText(`Transmitted: ${(100 * trans).toFixed(0)}%`, 8, 6);
        }
        if (drawMode) {
            ctx.textAlign = "center";
            ctx.fillStyle = C.rho;
            ctx.fillText("Drag to paint the potential, drag low to erase", w / 2, 6);
        } else if (finished) {
            // fade in and drift up slightly (ease-out)
            const t = Math.min(1, (performance.now() - finishedAt) / FADE_MS);
            const a = 1 - (1 - t) ** 3;
            ctx.save();
            ctx.globalAlpha = a;
            ctx.textAlign = "center";
            ctx.textBaseline = "middle";
            const lines = ["The wavepacket has left the boundaries.", "Click to place a new one, or press Play."];
            const yMid = (top + y0) / 2 + 8 * (1 - a);
            const oneLine = ctx.measureText(lines.join(" ")).width < w - 48;
            const shown = oneLine ? [lines.join(" ")] : lines;
            // soft backing so the message stays readable over barriers and walls
            const bw = Math.max(...shown.map((s) => ctx.measureText(s).width)) + 28;
            const bh = shown.length * 18 + 14;
            ctx.fillStyle = C.backing;
            ctx.beginPath();
            ctx.roundRect(w / 2 - bw / 2, yMid - bh / 2, bw, bh, 6);
            ctx.fill();
            ctx.fillStyle = C.text;
            shown.forEach((s, i) => ctx.fillText(s, w / 2, yMid + (i - (shown.length - 1) / 2) * 18));
            ctx.restore();
        }

        if (drag) drawPreview();
        strips.forEach((st) => {
            if (st.cv.hidden) return;
            const sw = st.cv.clientWidth, sh = st.cv.clientHeight;
            st.cx.clearRect(0, 0, sw, sh);
            st.draw(st.cx, sw, sh);
        });
    }

    // Space-time strip: one row of pixels per `history.dt` of simulated time, the
    // newest at the top (right under the plot it continues), older ones scrolling
    // down. Kept in an offscreen canvas that is shifted by a row as time goes on.
    // Mode "change" shows the density minus its value at t = 0 instead (purple where
    // denser, orange where thinner): small waves on a big condensate stand out.
    const history = { cv: document.createElement("canvas"), dt: 0.5, next: 0, mode: "density", ref: null };
    history.cx = history.cv.getContext("2d");

    // the density at t = 0, for the "change" mode
    function snapshotHistory() {
        history.ref = Float64Array.from(pR, (a, j) => a * a + pI[j] * pI[j]);
    }

    function resetHistory() {
        history.next = time;
        history.cv.width = ui.hcanvas.clientWidth || 1;
        history.cv.height = ui.hcanvas.clientHeight || 1;
        history.cx.clearRect(0, 0, history.cv.width, history.cv.height);
    }

    // After a light/dark flip the rows already recorded still carry the old theme's
    // colours. Every pixel is one of the two strip colours with its strength in the
    // alpha, so swap each to the matching new colour and keep the alpha.
    function recolorHistory(oldRho, oldDip) {
        const hw = history.cv.width, hh = history.cv.height;
        if (hw < 2 || hh < 2) return;
        if (String(oldRho) === String(rhoRGB) && String(oldDip) === String(dipRGB)) return;
        const img = history.cx.getImageData(0, 0, hw, hh), p = img.data;
        const dist = (i, c) => (p[i] - c[0]) ** 2 + (p[i + 1] - c[1]) ** 2 + (p[i + 2] - c[2]) ** 2;
        for (let i = 0; i < p.length; i += 4) {
            if (p[i + 3] === 0) continue;
            const c = dist(i, oldDip) < dist(i, oldRho) ? dipRGB : rhoRGB;
            p[i] = c[0]; p[i + 1] = c[1]; p[i + 2] = c[2];
        }
        history.cx.putImageData(img, 0, 0);
    }

    // called after every frame's steps: add the rows that have come due
    function recordHistory() {
        if (ui.hcanvas.hidden) return;
        const hw = history.cv.width, hh = history.cv.height;
        if (hw !== ui.hcanvas.clientWidth || hh !== ui.hcanvas.clientHeight) resetHistory();
        while (time >= history.next) {
            history.next += history.dt;
            history.cx.globalCompositeOperation = "copy";
            history.cx.drawImage(history.cv, 0, 1);
            history.cx.globalCompositeOperation = "source-over";
            const row = history.cx.createImageData(hw, 1);
            const change = history.mode === "change" && history.ref && history.ref.length === N;
            for (let X = 0; X < hw; X++) {
                const j = Math.max(0, Math.min(N - 1, Math.round((unpx((X / hw) * w) - XMIN) / dx - 1)));
                const r = pR[j] * pR[j] + pI[j] * pI[j];
                const d = change ? (r - history.ref[j]) / (0.25 * rho0) : r / rho0;
                const v = Math.min(1, Math.abs(d)) ** 0.6;
                row.data.set([...(d < 0 ? dipRGB : rhoRGB), Math.round(255 * v)], 4 * X);
            }
            history.cx.putImageData(row, 0, 0);
        }
    }

    function drawHistory(cx, sw, sh) {
        cx.drawImage(history.cv, 0, 0, sw, sh);
        cx.font = "11px 'Source Sans 3', sans-serif";
        cx.textBaseline = "top";
        cx.textAlign = "left";
        cx.fillStyle = C.text;
        cx.fillText(history.mode === "change"
            ? "Space–time, change since t = 0 (purple: denser, orange: thinner), newest at the top"
            : "Space–time: newest at the top", 8, 3);
        // time ticks down the right edge
        const span = sh * history.dt, every = [5, 10, 20, 50, 100, 200].find((e) => span / e <= 6) || 500;
        cx.textAlign = "right";
        cx.textBaseline = "middle";
        cx.strokeStyle = C.axis;
        for (let tt = every; tt < span; tt += every) {
            const y = Math.round(tt / history.dt) + 0.5;
            cx.beginPath();
            cx.moveTo(sw - 6, y);
            cx.lineTo(sw, y);
            cx.stroke();
            cx.fillText(`−${tt}`, sw - 9, y);
        }
    }

    // Phase strip: the phase of psi on a fixed axis from 0 to 2 pi, wrapping around:
    // a winding current is a sawtooth, wrapping once per unit of circulation. The phase
    // is measured from the cloud's left edge (the overall phase just rotates in time),
    // which sits at pi, in the middle: a dark soliton's step then drops towards 0 or
    // climbs towards 2 pi, so its direction stays visible. Drawn only where there is
    // density; where the phase wraps, the line leaves one edge and comes back at the other.
    let phaseU = null;
    function drawPhase(cx, sw, sh) {
        if (!phaseU || phaseU.length !== N) phaseU = new Float64Array(N);
        const u = phaseU, TWO_PI = 2 * Math.PI;
        let peak = 0;
        for (let j = j0; j <= j1; j++) peak = Math.max(peak, pR[j] * pR[j] + pI[j] * pI[j]);
        const dens = (j) => pR[j] * pR[j] + pI[j] * pI[j];
        let ref = j0;
        while (ref < j1 && dens(ref) < 0.1 * peak) ref++;
        // phase increment between neighbours j-1 -> j
        const inc = (j) => Math.atan2(pR[j - 1] * pI[j] - pI[j - 1] * pR[j], pR[j - 1] * pR[j] + pI[j - 1] * pI[j]);
        // Unwrapped first (each step between neighbours taken as the smallest one, so a
        // soliton's step shows the direction it moves in), then wrapped while drawing.
        u[ref] = Math.PI;
        for (let j = ref + 1; j <= j1; j++) u[j] = u[j - 1] + inc(j);
        for (let j = ref - 1; j >= j0; j--) u[j] = u[j + 1] - inc(j + 1);

        const pad = 14, Y = (th) => sh - pad - (th / TWO_PI) * (sh - 2 * pad);
        const turns = (v) => Math.floor(v / TWO_PI);
        cx.font = "11px 'Source Sans 3', sans-serif";
        cx.strokeStyle = C.axis;
        for (let m = 0; m <= 2; m++) {
            const y = Math.round(Y(m * Math.PI)) + 0.5;
            cx.setLineDash(m === 1 ? [2, 3] : []);
            cx.beginPath();
            cx.moveTo(0, y);
            cx.lineTo(sw, y);
            cx.stroke();
        }
        cx.setLineDash([]);
        cx.fillStyle = C.text;
        cx.textAlign = "left";
        cx.textBaseline = "top";
        cx.fillText("Phase of ψ", 8, 3);

        cx.strokeStyle = C.re;
        cx.lineWidth = 1.6;
        cx.beginPath();
        let prev = -1;   // the last point drawn, or -1 after a gap in the density
        for (let j = j0; j <= j1; j++) {
            if (dens(j) < 0.02 * peak) {
                prev = -1;
                continue;
            }
            const X = px(x[j]), v = u[j], w = v - turns(v) * TWO_PI;
            if (prev < 0) {
                cx.moveTo(X, Y(w));
            } else {
                const v0 = u[prev], n0 = turns(v0), n1 = turns(v);
                if (n0 !== n1) {
                    // crossed a multiple of 2 pi: run to that edge and come back at the other
                    const up = n1 > n0, edge = (up ? n1 : n0) * TWO_PI;
                    const Xm = px(x[prev]) + ((edge - v0) / (v - v0)) * (X - px(x[prev]));
                    cx.lineTo(Xm, Y(up ? TWO_PI : 0));
                    cx.moveTo(Xm, Y(up ? 0 : TWO_PI));
                }
                cx.lineTo(X, Y(w));
            }
            prev = j;
        }
        cx.stroke();

        // tick labels on a soft backing, so they stay readable over the curve
        cx.textAlign = "right";
        cx.textBaseline = "middle";
        ["0", "π", "2π"].forEach((t, m) => {
            const y = Math.round(Y(m * Math.PI)) + 0.5;
            const tw = cx.measureText(t).width;
            cx.fillStyle = C.backing;
            cx.fillRect(sw - tw - 8, y - 7, tw + 6, 14);
            cx.fillStyle = C.text;
            cx.fillText(t, sw - 4, y);
        });
    }

    // Momentum strip: |phi(k)|^2 over the preset's range of k. On a ring only
    // k = 2 pi m / L exist, so it is drawn as bars, one per winding number m.
    let pPeak = 0;
    const kBuf = { key: "", ks: null, P: null, P0: null };   // P0: the distribution at t = 0
    function drawMomentum(kctx, kw, kh) {
        const K = presets[presetKey].kRange || 5;
        const dk = periodic ? (2 * Math.PI) / RING_L : (2 * K) / 600;
        const key = `${K}|${periodic}|${N}`;
        if (kBuf.key !== key) {
            const n = Math.floor(K / dk);
            kBuf.ks = Float64Array.from({ length: 2 * n + 1 }, (_, i) => (i - n) * dk);
            kBuf.P = new Float64Array(2 * n + 1);
            kBuf.P0 = null;
            kBuf.key = key;
        }
        momentumDistribution(kBuf.ks, kBuf.P);
        if (!kBuf.P0) kBuf.P0 = Float64Array.from(kBuf.P);
        pPeak = Math.max(pPeak, ...kBuf.P);
        const base = kh - 16, X = (k) => ((k + K) / (2 * K)) * kw, Y = (v) => base - (v / (pPeak || 1)) * (base - 16);
        // axis and ticks
        kctx.strokeStyle = C.axis;
        kctx.beginPath();
        kctx.moveTo(0, base + 0.5);
        kctx.lineTo(kw, base + 0.5);
        kctx.stroke();
        kctx.font = "11px 'Source Sans 3', sans-serif";
        kctx.fillStyle = C.text;
        kctx.textAlign = "center";
        kctx.textBaseline = "top";
        const tick = (k, label) => {
            if (X(k) > 10 && X(k) < kw - 10) kctx.fillText(label, X(k), base + 3);
        };
        const byWinding = periodic && K <= 1;
        if (byWinding) {
            for (let m = -10; m <= 10; m += 2) tick((2 * Math.PI * m) / RING_L, String(m));
        } else {
            const stepK = K > 3 ? 1 : 0.5;
            for (let k = -Math.floor(K / stepK) * stepK; k <= K + 1e-9; k += stepK) tick(k, stepK === 1 ? String(Math.round(k)) : k.toFixed(1));
        }
        kctx.textAlign = "left";
        kctx.fillText(byWinding ? "Momentum k = 2πm/L, by winding m" : "Momentum distribution, k (dashed: at t = 0)", 8, 3);
        // where it started, for comparison
        if (!periodic && time > 0) {
            kctx.save();
            kctx.setLineDash([3, 3]);
            kctx.strokeStyle = C.potLine;
            kctx.beginPath();
            kBuf.ks.forEach((k, i) => (i ? kctx.lineTo(X(k), Y(kBuf.P0[i])) : kctx.moveTo(X(k), Y(kBuf.P0[i]))));
            kctx.stroke();
            kctx.restore();
        }
        // the distribution
        kctx.fillStyle = C.rhoFill;
        kctx.strokeStyle = C.rho;
        if (periodic) {
            const bw = Math.max(3, X(dk) - X(0) - 3);
            kBuf.ks.forEach((k, i) => kctx.fillRect(X(k) - bw / 2, Y(kBuf.P[i]), bw, base - Y(kBuf.P[i])));
        } else {
            kctx.beginPath();
            kctx.moveTo(X(kBuf.ks[0]), base);
            kBuf.ks.forEach((k, i) => kctx.lineTo(X(k), Y(kBuf.P[i])));
            kctx.lineTo(X(kBuf.ks[kBuf.ks.length - 1]), base);
            kctx.fill();
            kctx.lineWidth = 1.4;
            kctx.stroke();
        }
    }

    // ------------------------------------------------------------------
    // pointer: place / throw a packet (or place an imprint), or paint the potential
    // ------------------------------------------------------------------
    let drag = null;
    let painting = false, paintLast = 0, stirring = false;

    // what a drag sets: k for packets, v/c for a dark soliton, nothing for a bump
    function dragKey() {
        if (prep.kind === "packet") return "k";
        if (prep.kind === "ground" && prep.imprint === "soliton") return "v";
        return null;
    }
    function dragValue() {
        const key = dragKey(), f = fields[key];
        const units = unpx(drag.x) - unpx(drag.x0);
        const raw = key === "v" ? units / 15 : units / 3;
        const v = Math.round(raw / f.step) * f.step;
        return Math.max(f.min, Math.min(f.max, Number(v.toFixed(2))));
    }

    function drawPreview() {
        const X0 = drag.x0;
        ctx.save();
        ctx.setLineDash([3, 3]);
        ctx.strokeStyle = C.rho;
        ctx.beginPath();
        if (prep.kind === "packet") {
            // the packet's shape at the chosen spot
            const shape = { ...prep, x0: unpx(X0), k: 0 }, prof = [];
            let peak = 0;
            for (let X = 0; X <= w; X += 2) {
                const [a, b] = packetAmp(unpx(X), shape);
                prof.push(a * a + b * b);
                peak = Math.max(peak, a * a + b * b);
            }
            prof.forEach((r, i) => {
                const Y = y0 - (r / (peak || 1)) * (y0 - top) * 0.5;
                if (i === 0) ctx.moveTo(0, Y); else ctx.lineTo(i * 2, Y);
            });
        } else {
            // an imprint: mark where it goes
            ctx.moveTo(X0, top);
            ctx.lineTo(X0, y0);
        }
        ctx.stroke();
        ctx.restore();

        const key = dragKey(), val = key ? dragValue() : 0;
        if (val !== 0) {
            const Y = y0 - (y0 - top) * 0.6, X1 = drag.x;
            ctx.beginPath();
            ctx.moveTo(X0, Y);
            ctx.lineTo(X1, Y);
            ctx.lineTo(X1 - Math.sign(val) * 7, Y - 4);
            ctx.moveTo(X1, Y);
            ctx.lineTo(X1 - Math.sign(val) * 7, Y + 4);
            ctx.strokeStyle = C.rho;
            ctx.lineWidth = 1.5;
            ctx.stroke();
            ctx.fillStyle = C.rho;
            ctx.textAlign = "center";
            ctx.textBaseline = "bottom";
            ctx.fillText(key === "v" ? `v = ${val.toFixed(2)}c` : `k = ${val.toFixed(1)}`, (X0 + X1) / 2, Y - 6);
        }
    }

    function paint(fromX, toX, y) {
        const v = Math.max(0, unyE(y));
        if (presetKey !== "custom") {
            presets.custom.scale = presets[presetKey].scale;   // keep the current view
            presets.custom.absorbing = absorbing;               // and the boundary condition
            selectPreset("custom", false);
            syncPrepUI();
            setActiveDemo(-1);
        }
        const a = Math.max(VMIN, Math.min(fromX, toX) - 0.5), b = Math.min(VMAX, Math.max(fromX, toX) + 0.5);
        for (let j = j0; j <= j1; j++) if (x[j] >= a && x[j] <= b) V[j] = v < 0.05 ? 0 : v;
        potVersion++;
        refresh();
        potentialChanged();
        draw();
    }

    function localXY(e) {
        const r = canvas.getBoundingClientRect();
        return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    // the plot handles its own drags: no text selection or native drag-and-drop ghost image
    canvas.addEventListener("dragstart", (e) => e.preventDefault());

    canvas.addEventListener("pointerdown", (e) => {
        if (e.pointerType === "mouse" && e.button !== 0) return;
        e.preventDefault();
        window.getSelection()?.removeAllRanges();
        const p = localXY(e);
        try {
            canvas.setPointerCapture(e.pointerId);
        } catch (_) { /* synthetic events */ }
        if (drawMode) {
            drag = null;
            paintLast = unpx(p.x);
            paint(paintLast, paintLast, p.y);
            painting = true;
            return;
        }
        if (presets[presetKey].stir) {
            stirring = true;
            stir.target = Math.max(-32, Math.min(32, unpx(p.x)));
            return;
        }
        if (!placeable()) return;   // nothing to place: position fixed by the potential
        drag = { x0: p.x, x: p.x };
        draw();
    });

    canvas.addEventListener("pointermove", (e) => {
        const p = localXY(e);
        if (stirring) {
            stir.target = Math.max(-32, Math.min(32, unpx(p.x)));
        } else if (painting) {
            const xx = unpx(p.x);
            paint(paintLast, xx, p.y);
            paintLast = xx;
        } else if (drag) {
            drag.x = p.x;
            draw();
        }
    });

    canvas.addEventListener("pointerup", () => {
        stirring = false;
        if (painting) {
            painting = false;
            return;
        }
        if (!drag) return;
        const moved = Math.abs(drag.x - drag.x0) > 6;
        prep.x0 = unpx(drag.x0);
        const key = dragKey();
        if (key && moved) prep[key] = dragValue();
        drag = null;
        syncPrepUI();
        markModified();
        launch();
    });
    canvas.addEventListener("pointercancel", () => {
        drag = null;
        painting = false;
        stirring = false;
        draw();
    });

    // ------------------------------------------------------------------
    // main loop (only runs while visible)
    // ------------------------------------------------------------------
    let visible = false, running = false, frameCount = 0, sinceStatus = 0;

    function loop() {
        if (running || !playing || !onScreen()) return;
        running = true;
        requestAnimationFrame(function frame() {
            if (!playing || !onScreen()) {
                running = false;
                return;
            }
            // same simulated time per frame whatever dt is
            const steps = Math.max(1, Math.round((Number(ui.speed.value) * DEFAULT_DT) / dt));
            for (let s = 0; s < steps; s++) {
                step();
                stepClassical();
            }
            time += steps * dt;
            sinceStatus += steps * dt;
            recordHistory();

            if (imaginary) {
                const m = measure();
                E = m.E;
                // grow the density scale smoothly as the state sharpens
                let peak = 0;
                for (let j = j0; j <= j1; j++) peak = Math.max(peak, pR[j] * pR[j] + pI[j] * pI[j]);
                if (peak > rho0) rho0 = peak;
                converged = Math.abs(lastE - E) / (steps * dt) < 1e-5;
                lastE = E;
                updateStatus(m);
            } else if (presets[presetKey].timeDep && frameCount % 6 !== 5) {
                // a time-dependent potential does work on the system: track its energy
                // every frame, or a fast drive makes the energy line jump
                frameCount++;
                E = measure().E;
            } else if (++frameCount % 6 === 0) {
                const m = measure();
                if (presets[presetKey].timeDep) E = m.E;
                const ok = updateStatus(m, sinceStatus);
                sinceStatus = 0;
                // stop once essentially all probability has been absorbed at the edges
                if (absorbing && m.norm < 2e-3) {
                    finished = true;
                    finishedAt = window.matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : performance.now();
                    requestAnimationFrame(fadeIn);
                }
                if (!ok || finished) {
                    playing = false;
                    updatePlayButton();
                }
            }
            draw();
            requestAnimationFrame(frame);
        });
    }

    // redraw until the end message has fully faded in (the main loop is stopped by then)
    function fadeIn() {
        if (!finished || !onScreen()) return;
        draw();
        if (performance.now() - finishedAt < FADE_MS) requestAnimationFrame(fadeIn);
    }

    // Only simulate and render while at least 10% of the widget is on screen
    // and the tab is in the foreground.
    function onScreen() {
        return visible && !document.hidden;
    }

    function wake() {
        if (!onScreen()) return;
        draw();
        loop();
    }

    new IntersectionObserver((entries) => {
        visible = entries[0].intersectionRatio >= 0.1;
        wake();
    }, { threshold: [0, 0.1] }).observe(root);

    document.addEventListener("visibilitychange", wake);

    new ResizeObserver(resize).observe(canvas);

    // light/dark toggle: Quarto swaps the stylesheet and the body's class
    new MutationObserver(() => requestAnimationFrame(() => {
        const oldRho = rhoRGB, oldDip = dipRGB;
        readColors();
        recolorHistory(oldRho, oldDip);
        draw();
    })).observe(document.body, { attributes: true, attributeFilter: ["class"] });

    if (window.MathJax && window.MathJax.typesetPromise) window.MathJax.typesetPromise([root]);

    updatePlayButton();
    updateNumericsNote();
    prep = { ...defaultPrep };
    if (demos.length) {
        showTopic(demos[0].topic);
        applyDemo(0, playing);   // open on the first demo (paused if the visitor prefers reduced motion)
    } else {
        root.querySelector(".tdse-demos").hidden = true;
        selectPreset(presetKey);
    }
    resize();
})();
