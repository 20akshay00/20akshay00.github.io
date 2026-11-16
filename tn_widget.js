// Compressing a many-body state into a matrix product state.
//
// The exact ground state of a small Bose-Hubbard chain (6 bosons on 6 sites,
// 462 amplitudes, found by Lanczos) or a random state with the same particle
// number is cut into an MPS by successive SVDs (Schmidt decompositions), keeping
// at most D singular values at every cut. The widget steps through that
// construction, and plots every amplitude of the compressed state against the
// exact one, next to the simpler alternative of keeping the largest amplitudes.
(function () {
    const root = document.getElementById("tn-widget");
    if (!root) return;

    // ------------------------------------------------------------------
    // the many-body basis: L sites holding NB bosons
    // ------------------------------------------------------------------
    const L = 6, NB = 6;

    // basis[s][M] lists the occupations of s sites holding M bosons, with the
    // first site outermost: the configuration (n, rest) is found at
    // offset[s][M][n] + (index of rest in basis[s-1][M-n]).
    const basis = [[]], offset = [[]];
    for (let M = 0; M <= NB; M++) basis[0].push(M === 0 ? [[]] : []);
    for (let s = 1; s <= L; s++) {
        basis.push([]);
        offset.push([]);
        for (let M = 0; M <= NB; M++) {
            const list = [], off = [];
            for (let n = 0; n <= M; n++) {
                off.push(list.length);
                for (const rest of basis[s - 1][M - n]) list.push([n, ...rest]);
            }
            basis[s].push(list);
            offset[s].push(off);
        }
    }
    const configs = basis[L][NB], DIM = configs.length;

    // ------------------------------------------------------------------
    // Bose-Hubbard Hamiltonian (open chain, J = 1):
    //   H = -sum_i (b+_i b_i+1 + h.c.) + U/2 sum_i n_i (n_i - 1)
    // stored as the hopping matrix elements plus the diagonal interaction
    // ------------------------------------------------------------------
    const keyOf = (c) => c.reduce((k, n) => k * (NB + 1) + n, 0);
    const indexOf = new Map(configs.map((c, i) => [keyOf(c), i]));
    const pairs = new Float64Array(DIM);         // sum_i n_i (n_i - 1) / 2
    const hop = { from: [], to: [], amp: [] };
    configs.forEach((c, i) => {
        pairs[i] = c.reduce((e, n) => e + (n * (n - 1)) / 2, 0);
        for (let b = 0; b + 1 < L; b++) {
            if (!c[b]) continue;
            const d = c.slice();
            d[b]--;
            d[b + 1]++;
            const j = indexOf.get(keyOf(d)), amp = -Math.sqrt(c[b] * (c[b + 1] + 1));
            hop.from.push(i, j);
            hop.to.push(j, i);
            hop.amp.push(amp, amp);
        }
    });
    hop.from = Int32Array.from(hop.from);
    hop.to = Int32Array.from(hop.to);
    hop.amp = Float64Array.from(hop.amp);

    function applyH(U, v, out) {
        for (let i = 0; i < DIM; i++) out[i] = U * pairs[i] * v[i];
        for (let k = 0; k < hop.amp.length; k++) out[hop.to[k]] += hop.amp[k] * v[hop.from[k]];
    }

    // ------------------------------------------------------------------
    // linear algebra
    // ------------------------------------------------------------------
    const dot = (a, b) => {
        let s = 0;
        for (let i = 0; i < a.length; i++) s += a[i] * b[i];
        return s;
    };
    const axpy = (a, x, y) => {
        for (let i = 0; i < y.length; i++) y[i] += a * x[i];
    };

    // Eigen-decomposition of a small symmetric matrix A (n x n, row-major, which
    // is overwritten) by cyclic Jacobi rotations. Returns the eigenvalues and the
    // eigenvectors as the columns of V.
    function eigh(A, n) {
        const V = new Float64Array(n * n);
        for (let i = 0; i < n; i++) V[i * n + i] = 1;
        for (let sweep = 0; sweep < 60; sweep++) {
            let off = 0, diag = 0;
            for (let p = 0; p < n; p++) {
                diag += A[p * n + p] ** 2;
                for (let q = p + 1; q < n; q++) off += A[p * n + q] ** 2;
            }
            if (off <= 1e-30 * diag || off === 0) break;
            for (let p = 0; p < n; p++) {
                for (let q = p + 1; q < n; q++) {
                    const apq = A[p * n + q];
                    if (Math.abs(apq) < 1e-300) continue;
                    const theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
                    const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
                    const c = 1 / Math.sqrt(t * t + 1), s = t * c;
                    for (let k = 0; k < n; k++) {
                        const akp = A[k * n + p], akq = A[k * n + q];
                        A[k * n + p] = c * akp - s * akq;
                        A[k * n + q] = s * akp + c * akq;
                    }
                    for (let k = 0; k < n; k++) {
                        const apk = A[p * n + k], aqk = A[q * n + k];
                        A[p * n + k] = c * apk - s * aqk;
                        A[q * n + k] = s * apk + c * aqk;
                    }
                    for (let k = 0; k < n; k++) {
                        const vkp = V[k * n + p], vkq = V[k * n + q];
                        V[k * n + p] = c * vkp - s * vkq;
                        V[k * n + q] = s * vkp + c * vkq;
                    }
                }
            }
        }
        return { values: Array.from({ length: n }, (_, i) => A[i * n + i]), V };
    }

    // Ground state by Lanczos with full reorthogonalisation, started from the
    // previous ground state (so dragging the U slider converges in a few steps).
    function groundState(U, start) {
        const basisV = [], a = [], b = [];
        let v = Float64Array.from(start);
        axpy(1 / Math.sqrt(dot(v, v)) - 1, v, v);
        const w = new Float64Array(DIM);
        let best = null;
        for (let j = 0; j < 150; j++) {
            basisV.push(v);
            applyH(U, v, w);
            a.push(dot(w, v));
            for (let pass = 0; pass < 2; pass++) for (const u of basisV) axpy(-dot(w, u), u, w);
            const beta = Math.sqrt(dot(w, w));
            const m = basisV.length;
            if (m % 5 === 0 || beta < 1e-10 || m === DIM) {
                const T = new Float64Array(m * m);
                for (let i = 0; i < m; i++) {
                    T[i * m + i] = a[i];
                    if (i + 1 < m) T[i * m + i + 1] = T[(i + 1) * m + i] = b[i];
                }
                const { values, V } = eigh(T, m);
                let lo = 0;
                for (let i = 1; i < m; i++) if (values[i] < values[lo]) lo = i;
                best = { E: values[lo], y: Array.from({ length: m }, (_, i) => V[i * m + lo]) };
                if (beta * Math.abs(best.y[m - 1]) < 1e-9 || beta < 1e-10) break;
            }
            b.push(beta);
            v = Float64Array.from(w, (x) => x / beta);
        }
        const psi = new Float64Array(DIM);
        best.y.forEach((c, i) => axpy(c, basisV[i], psi));
        const norm = Math.sqrt(dot(psi, psi));
        for (let i = 0; i < DIM; i++) psi[i] /= norm;
        return { psi, E: best.E };
    }

    // a random state with the same particle number (Gaussian amplitudes)
    function randomState(seed) {
        let t = seed >>> 0;
        const rand = () => {
            t = (t + 0x6d2b79f5) >>> 0;
            let r = Math.imul(t ^ (t >>> 15), 1 | t);
            r = (r + Math.imul(r ^ (r >>> 7), 61 | r)) ^ r;
            return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
        };
        const psi = Float64Array.from({ length: DIM }, () => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand()));
        const norm = Math.sqrt(dot(psi, psi));
        return psi.map((x) => x / norm);
    }

    // ------------------------------------------------------------------
    // Compression into a matrix product state, cutting the chain from left to
    // right. At each cut the state is written as sum_a s_a |left_a>|right_a>
    // (a Schmidt decomposition, done block by block in the number of bosons on
    // the left) and only the D largest s_a are kept. The kept left states are
    // the MPS tensors; what remains of the state, R[a] = s_a |right_a>, moves on
    // to the next cut. Each step projects the state, so after the last site the
    // squared norm of what is left is the overlap |<psi|psi_D>|^2 with the
    // (normalised) compressed state.
    //
    // Returns, per bond, the squared Schmidt values seen there and how many were
    // kept, the overlap, how many numbers the MPS stores, and expand(), which
    // multiplies the tensors back out into the 6435 amplitudes of psi_D.
    // ------------------------------------------------------------------
    function compress(psi, D) {
        let q = [0], R = [psi];                 // left boson number and remainder per kept state
        const bonds = [], steps = [];
        let stored = 0;
        for (let k = 0; k < L; k++) {
            const s = L - k;                   // sites still to the right, including site k
            const blocks = new Map();          // left boson number Q -> rows (a, n)
            q.forEach((qa, a) => {
                for (let n = 0; n <= NB - qa; n++) {
                    const Q = qa + n;
                    if (!blocks.has(Q)) blocks.set(Q, []);
                    blocks.get(Q).push([a, n]);
                }
            });
            const cand = [];
            for (const [Q, rows] of blocks) {
                const cols = basis[s - 1][NB - Q].length, r = rows.length;
                if (!cols) continue;
                const M = rows.map(([a, n]) => {
                    const o = offset[s][NB - q[a]][n];
                    return R[a].subarray(o, o + cols);
                });
                const G = new Float64Array(r * r);
                for (let i = 0; i < r; i++) for (let j = i; j < r; j++) G[i * r + j] = G[j * r + i] = dot(M[i], M[j]);
                const { values, V } = eigh(G, r);
                values.forEach((wt, e) => {
                    if (wt < 1e-15) return;
                    const u = Float64Array.from(rows, (_, i) => V[i * r + e]);
                    const vec = new Float64Array(cols);
                    for (let i = 0; i < r; i++) axpy(u[i], M[i], vec);
                    cand.push({ wt, Q, vec, rows, u });
                });
            }
            cand.sort((u, v) => v.wt - u.wt);
            const keep = cand.slice(0, D);
            for (const c of keep) stored += c.rows.length;
            if (k + 1 < L) bonds.push({ weights: cand.map((c) => c.wt), kept: keep.length });
            steps.push({ s, q, keep });
            q = keep.map((c) => c.Q);
            R = keep.map((c) => c.vec);
        }
        const overlap = R.reduce((acc, r) => acc + r[0] * r[0], 0);

        // back from the last cut to the first: each remainder is rebuilt from the
        // kept states, which projects psi onto what the MPS can describe
        function expand() {
            let next = R;
            for (let k = L - 1; k >= 0; k--) {
                const { s, q: qPrev, keep } = steps[k];
                const prev = qPrev.map((qa) => new Float64Array(basis[s][NB - qa].length));
                keep.forEach((c, b) => {
                    const cols = next[b].length;
                    c.rows.forEach(([a, n], i) => {
                        const o = offset[s][NB - qPrev[a]][n];
                        axpy(c.u[i], next[b], prev[a].subarray(o, o + cols));
                    });
                });
                next = prev;
            }
            return next[0];
        }
        return { bonds, overlap, stored, expand };
    }

    // ---- end of physics

    // ------------------------------------------------------------------
    // markup
    // ------------------------------------------------------------------
    const D_MAX = 30;                       // the largest bond dimension any state here needs
    const fmtInt = (n) => n.toLocaleString("en-US");
    const sub = (k) => (L > 9 ? `{${k}}` : k);

    // The diagram is a grid: site columns (boxes, legs, labels) alternate with
    // bond columns. While stepping through the construction, the part not yet
    // cut is one wide box laid over the remaining sites.
    const cells = [];
    for (let k = 1; k <= L; k++) {
        const c = 2 * k - 1;
        cells.push(`<div class="tn-box" data-site="${k}" style="grid-column:${c}">\\(C^{(${k})}\\)</div>`);
        cells.push(`<div class="tn-leg" style="grid-column:${c}"></div>`);
        cells.push(`<div class="tn-n" style="grid-column:${c}">\\(n_${sub(k)}\\)</div>`);
        if (k < L) cells.push(`<div class="tn-bond" style="grid-column:${c + 1}"><span data-bond="${k - 1}"></span></div>`);
    }
    const columns = Array.from({ length: 2 * L - 1 }, (_, i) => (i % 2 ? "minmax(0.6rem, 1fr)" : "var(--tn-site)")).join(" ");

    root.classList.add("tdse", "tn");
    root.innerHTML = `
        <div class="tdse-stage">
            <div class="tn-mps" role="img" aria-label="The amplitudes C as a matrix product state: one tensor per site, with the bond dimensions written on the bonds">
                <div class="tn-lhs">\\(C_{n_1 \\dots n_${L}} =\\)</div>
                <div class="tn-chain" style="grid-template-columns:${columns}">
                    ${cells.join("\n                    ")}
                    <div class="tn-box tn-rest" data-ui="rest"><span data-ui="restFull">\\(C\\)</span><span data-ui="restSV" hidden>\\(S V^\\dagger\\)</span></div>
                </div>
            </div>
            <div class="tn-build">
                <div class="tn-build-bar">
                    <span class="tn-build-title">How it is built</span>
                    <div class="tdse-seg" role="group" aria-label="Step of the construction" data-ui="steps">
                        ${Array.from({ length: L }, (_, k) => `<button type="button" data-step="${k}" aria-pressed="false">${k ? `Cut ${k}` : "Start"}</button>`).join("")}
                    </div>
                    <button type="button" class="tdse-btn" data-ui="next"><span>Next</span><i class="bi bi-chevron-right" aria-hidden="true"></i></button>
                </div>
                <p class="tdse-note tn-stable" data-ui="buildNote"></p>
            </div>
        </div>
        <details class="tdse-controls tn-why">
        <summary class="tdse-controls-toggle"><i class="bi bi-lightbulb" aria-hidden="true"></i><span>Why does it work?</span><small>how few amplitudes matter, and what a small D keeps</small></summary>
        <div class="tdse-stage tn-why-stage">
            <p class="tn-amp-cap">Every amplitude <span class="math inline">\\(C_{n_1 \\dots n_${L}} = \\langle n_1 \\dots n_${L} | \\psi \\rangle\\)</span>, largest first: exact (grey) and from the MPS (teal)</p>
            <canvas class="tn-amp" aria-label="All amplitudes of the exact state, largest first, with those of the compressed state on top"></canvas>
        </div>
        <div class="tdse-panel tn-panel">
            <div class="tdse-group" role="group" aria-labelledby="tn-g-state">
                <div class="tdse-group-title" id="tn-g-state">The state</div>
                <p class="tdse-eq"><span class="math inline">\\(H = -J\\sum_i \\big(b_i^\\dagger b_{i+1} + \\mathrm{h.c.}\\big) + \\tfrac{U}{2}\\sum_i n_i(n_i-1)\\)</span></p>
                <p class="tdse-eq-sub">${L} bosons on ${L} sites: <span class="math inline">\\(|\\psi\\rangle = \\sum C_{n_1 \\dots n_${L}} |n_1 \\dots n_${L}\\rangle\\)</span> over ${fmtInt(DIM)} configurations</p>
                <div class="tdse-seg tdse-seg-full" role="group" aria-label="State" data-ui="kindSeg">
                    <button type="button" data-kind="ground" aria-pressed="true">Ground state</button>
                    <button type="button" data-kind="random" aria-pressed="false">Random state</button>
                </div>
                <label class="tdse-field" data-ui="uWrap"><span>Interaction <span class="math inline">\\(U/J\\)</span></span><input type="range" min="0" max="20" step="0.1" value="2" data-ui="U"><output data-ui="uOut">2.0</output></label>
                <p class="tdse-note tn-stable" data-ui="stateNote"></p>
            </div>
            <div class="tdse-group" role="group" aria-labelledby="tn-g-mps">
                <div class="tdse-group-title" id="tn-g-mps">Compression</div>
                <label class="tdse-field"><span>Bond dim. <span class="math inline">\\(D\\)</span></span><input type="range" min="1" max="${D_MAX}" step="1" value="${D_MAX}" data-ui="D"><output data-ui="dOut"></output></label>
                <p class="tn-params">Parameters in the MPS: <b data-ui="stored"></b></p>
                <div class="tn-score">
                    <div><div class="tn-score-v" data-ui="overlap"></div><div class="tn-score-l">Overlap of the MPS with the exact state</div></div>
                    <div class="tn-score-alt"><div class="tn-score-v" data-ui="topk"></div><div class="tn-score-l">Keeping the <span data-ui="topkN"></span> largest amplitudes instead</div></div>
                </div>
                <p class="tdse-note tn-stable" data-ui="mpsNote"></p>
            </div>
        </div>
        </details>
        <div class="tdse-credit">
            <svg class="tdse-claude" viewBox="0 0 24 24" aria-hidden="true"><path d="m4.7144 15.9555 4.7174-2.6471.079-.2307-.079-.1275h-.2307l-.7893-.0486-2.6956-.0729-2.3375-.0971-2.2646-.1214-.5707-.1215-.5343-.7042.0546-.3522.4797-.3218.686.0608 1.5179.1032 2.2767.1578 1.6514.0972 2.4468.255h.3886l.0546-.1579-.1336-.0971-.1032-.0972L6.973 9.8356l-2.55-1.6879-1.3356-.9714-.7225-.4918-.3643-.4614-.1578-1.0078.6557-.7225.8803.0607.2246.0607.8925.686 1.9064 1.4754 2.4893 1.8336.3643.3035.1457-.1032.0182-.0728-.164-.2733-1.3539-2.4467-1.445-2.4893-.6435-1.032-.17-.6194c-.0607-.255-.1032-.4674-.1032-.7285L6.287.1335 6.6997 0l.9957.1336.419.3642.6192 1.4147 1.0018 2.2282 1.5543 3.0296.4553.8985.2429.8318.091.255h.1579v-.1457l.1275-1.706.2368-2.0947.2307-2.6957.0789-.7589.3764-.9107.7468-.4918.5828.2793.4797.686-.0668.4433-.2853 1.8517-.5586 2.9021-.3643 1.9429h.2125l.2429-.2429.9835-1.3053 1.6514-2.0643.7286-.8196.85-.9046.5464-.4311h1.0321l.759 1.1293-.34 1.1657-1.0625 1.3478-.8804 1.1414-1.2628 1.7-.7893 1.36.0729.1093.1882-.0183 2.8535-.607 1.5421-.2794 1.8396-.3157.8318.3886.091.3946-.3278.8075-1.967.4857-2.3072.4614-3.4364.8136-.0425.0304.0486.0607 1.5482.1457.6618.0364h1.621l3.0175.2247.7892.522.4736.6376-.079.4857-1.2142.6193-1.6393-.3886-3.825-.9107-1.3113-.3279h-.1822v.1093l1.0929 1.0686 2.0035 1.8092 2.5075 2.3314.1275.5768-.3218.4554-.34-.0486-2.2039-1.6575-.85-.7468-1.9246-1.621h-.1275v.17l.4432.6496 2.3436 3.5214.1214 1.0807-.17.3521-.6071.2125-.6679-.1214-1.3721-1.9246L14.38 17.959l-1.1414-1.9428-.1397.079-.674 7.2552-.3156.3703-.7286.2793-.6071-.4614-.3218-.7468.3218-1.4753.3886-1.9246.3157-1.53.2853-1.9004.17-.6314-.0121-.0425-.1397.0182-1.4328 1.9672-2.1796 2.9446-1.7243 1.8456-.4128.164-.7164-.3704.0667-.6618.4008-.5889 2.386-3.0357 1.4389-1.882.929-1.0868-.0062-.1579h-.0546l-6.3385 4.1164-1.1293.1457-.4857-.4554.0608-.7467.2307-.2429 1.9064-1.3114Z"/></svg>
            <span>Made with the assistance of Claude Opus 5.5</span>
        </div>`;

    const ui = {};
    root.querySelectorAll("[data-ui]").forEach((el) => (ui[el.dataset.ui] = el));
    const bondLabels = Array.from(root.querySelectorAll("[data-bond]"));
    const siteBoxes = Array.from(root.querySelectorAll("[data-site]"));
    const amp = { cv: root.querySelector(".tn-amp") };
    amp.cx = amp.cv.getContext("2d");

    // the palette follows the light/dark theme (CSS custom properties, styles.scss)
    const css = getComputedStyle(root);
    const color = (name, fallback) => css.getPropertyValue(name).trim() || fallback;
    const C = {};
    function readColors() {
        Object.assign(C, {
            teal: color("--tdse-rho", "#006d77"),
            ink: color("--tdse-re", "#2c3e50"),
            axis: color("--tdse-axis", "rgba(44, 62, 80, 0.25)"),
        });
    }
    readColors();

    // ------------------------------------------------------------------
    // state
    // ------------------------------------------------------------------
    // build: how many cuts of the construction are shown (L - 1: the whole MPS)
    const S = { kind: "ground", U: 2, D: D_MAX, build: L - 1, psi: null, order: null, sq: null, mps: null, psiD: null };
    const groundCache = new Map();          // U (to 0.1) -> ground state, so revisiting is instant
    let lastGround = new Float64Array(DIM).fill(1);
    const random = randomState(20240607);

    function prepare() {
        if (S.kind === "random") {
            S.psi = random;
        } else {
            const key = S.U.toFixed(1);
            if (!groundCache.has(key)) {
                lastGround = groundState(S.U, lastGround).psi;
                groundCache.set(key, lastGround);
            }
            S.psi = groundCache.get(key);
        }
        // configurations ranked by |amplitude|, and the weight of the largest ones
        S.order = Array.from({ length: DIM }, (_, i) => i).sort((a, b) => Math.abs(S.psi[b]) - Math.abs(S.psi[a]));
        S.sq = new Float64Array(DIM);
        S.order.reduce((acc, i, r) => (S.sq[r] = acc + S.psi[i] ** 2), 0);
    }

    // ------------------------------------------------------------------
    // the diagram, at some step of its construction
    // ------------------------------------------------------------------
    function drawChain() {
        const k = S.build;                  // cuts made so far; sites k+1..L are still one tensor
        const merged = k < L - 1;
        siteBoxes.forEach((box, i) => (box.style.visibility = merged && i >= k ? "hidden" : ""));
        bondLabels.forEach((lab, i) => {
            const shown = i < k;
            lab.parentElement.style.visibility = shown || !merged ? "" : "hidden";
            lab.textContent = String(S.mps.bonds[i].kept);
            lab.classList.toggle("tn-cut", S.mps.bonds[i].kept < S.mps.bonds[i].weights.length);
        });
        ui.rest.hidden = !merged;
        ui.rest.style.gridColumn = `${2 * k + 1} / ${2 * L}`;
        ui.restFull.hidden = k > 0;
        ui.restSV.hidden = k === 0;
        ui.steps.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(S.stepped && Number(b.dataset.step) === k)));
        const again = S.stepped && k === L - 1;
        ui.next.querySelector("span").textContent = !S.stepped ? "Step through" : again ? "From the start" : "Next";
        ui.next.querySelector("i").className = `bi ${again ? "bi-arrow-counterclockwise" : "bi-chevron-right"}`;
        const versions = [buildNote(L - 1, true), ...Array.from({ length: L }, (_, j) => buildNote(j, false))];
        setNote(ui.buildNote, versions, k === L - 1 && !S.stepped ? 0 : k + 1);
    }

    function buildNote(k, intro) {
        const last = L - 1;
        if (intro) {
            return `One SVD at each cut, from left to right, keeping at most D singular values each time: those are the numbers on the bonds.`;
        }
        if (k === 0) {
            return `<b>Start.</b> Take all ${fmtInt(DIM)} amplitudes as a single tensor <i>C</i> with one leg per site.`;
        }
        const kept = S.mps.bonds[k - 1].kept, all = S.mps.bonds[k - 1].weights.length;
        const keep = kept < all ? `keep the ${kept} largest of its ${all} singular values` : `keep all ${all} nonzero singular values`;
        const rows = k === 1 ? "leg n<sub>1</sub>" : `incoming bond and leg n<sub>${k}</sub>`;
        return `<b>Cut ${k} of ${last}.</b> Read the uncut tensor as a matrix, its ${rows} against legs n<sub>${k + 1}</sub>…n<sub>${L}</sub>, and take its SVD, <i>U S V</i><sup>†</sup>. We ${keep}: <i>U</i> becomes <i>C</i><sup>(${k})</sup>, and ${k === last ? `<i>S V</i><sup>†</sup>, with a single leg left, is <i>C</i><sup>(${L})</sup>. Done.` : `<i>S V</i><sup>†</sup> is cut next.`}`;
    }

    // ------------------------------------------------------------------
    // the amplitudes, largest first: exact (grey) and from the MPS (teal)
    // ------------------------------------------------------------------
    function fit(c) {
        const dpr = window.devicePixelRatio || 1;
        c.w = c.cv.clientWidth;
        c.h = c.cv.clientHeight;
        c.cv.width = Math.round(c.w * dpr);
        c.cv.height = Math.round(c.h * dpr);
        c.cx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function drawAmps() {
        const { cx, w, h } = amp;
        cx.clearRect(0, 0, w, h);
        const left = 8, right = w - 8, top = 6, bottom = h - 6;
        const ex = S.order.map((i) => S.psi[i]), mp = S.order.map((i) => S.psiD[i]);
        let hi = 0, lo = 0;
        for (const v of ex) {
            hi = Math.max(hi, v);
            lo = Math.min(lo, v);
        }
        const span = hi - lo || 1;
        const X = (i) => left + ((i + 0.5) / DIM) * (right - left);
        const Y = (v) => bottom - ((Math.max(lo, Math.min(hi, v)) - lo) / span) * (bottom - top);

        cx.strokeStyle = C.axis;
        cx.lineWidth = 1;
        cx.beginPath();
        cx.moveTo(left, Math.round(Y(0)) + 0.5);
        cx.lineTo(right, Math.round(Y(0)) + 0.5);
        cx.stroke();

        const dot = (v, i, r) => {
            cx.beginPath();
            cx.arc(X(i), Y(v), r, 0, 2 * Math.PI);
            cx.fill();
        };
        cx.fillStyle = C.ink;
        cx.globalAlpha = 0.22;
        ex.forEach((v, i) => dot(v, i, 3.2));
        cx.fillStyle = C.teal;
        cx.globalAlpha = 0.9;
        mp.forEach((v, i) => dot(v, i, 1.4));
        cx.globalAlpha = 1;
    }

    // ------------------------------------------------------------------
    // readouts and notes
    // ------------------------------------------------------------------
    const sup = (n) => String(n).replace(/[-0-9]/g, (ch) => "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"["-0123456789".indexOf(ch)]);
    function fmtOverlap(F) {
        const err = 1 - F;
        if (err < 1e-12) return "1 (exact)";
        if (err > 0.01) return F.toFixed(3);
        const e = Math.floor(Math.log10(err)), m = err / 10 ** e;
        return `1 − ${m.toFixed(1)} × 10${sup(e)}`;
    }

    // A note whose text changes keeps the height of its longest version, so the
    // page below does not jump: all versions share one grid cell, one is visible.
    function setNote(el, versions, active) {
        if (el.children.length !== versions.length) el.innerHTML = versions.map(() => "<span></span>").join("");
        versions.forEach((html, i) => {
            const span = el.children[i];
            if (span.innerHTML !== html) span.innerHTML = html;
            span.style.visibility = i === active ? "" : "hidden";
            span.setAttribute("aria-hidden", String(i !== active));
        });
    }

    const STATE_NOTES = [
        "Every amplitude drawn at random: a typical state of the Hilbert space. Its weight is spread over all configurations, and it is nearly as entangled as a state can be.",
        "Superfluid: the atoms hop freely, so the weight is spread over hundreds of configurations. The phase is gapless (critical), so the entanglement does grow with the length of the chain, but only logarithmically: a modest D still captures it.",
        "Near the superfluid to Mott insulator transition (at U/J ≈ 3.3 in a long chain).",
        "Mott insulator: a gapped phase. The repulsion pins one atom to each site, so a few configurations carry almost all the weight, and the entanglement across a cut stays the same however long the chain (an area law).",
    ];
    const stateNote = () => (S.kind === "random" ? 0 : S.U < 2.5 ? 1 : S.U <= 4.5 ? 2 : 3);

    // what the comparison of the two cards shows
    function mpsNotes() {
        return [
            `Nothing is discarded yet, so the MPS reproduces every amplitude, with ${fmtInt(S.mps.stored)} parameters for ${fmtInt(DIM)} amplitudes. The saving comes from truncating: lower D.`,
            "No compression works for a random state, and the MPS does even worse than keeping the largest amplitudes: the entanglement across each cut is exactly what it cannot afford.",
            `With D = ${S.D} the bonds are too narrow for the entanglement across the cuts, and the MPS misses part of the state. Raise D.`,
            "Same number of parameters, very different results. The weight is spread out, so keeping the largest amplitudes fails; but little entanglement crosses each cut, and that is all an MPS needs.",
            "Both do well here: a few configurations carry almost all the weight, so even keeping the largest amplitudes works. The MPS stands out in the superfluid, where the weight is spread out.",
            `With ${fmtInt(S.mps.stored)} parameters both are close to exact. The difference shows at smaller D.`,
        ];
    }
    function mpsNote(topk) {
        if (S.mps.bonds.every((b) => b.kept === b.weights.length)) return 0;
        if (S.kind === "random") return 1;
        if (S.mps.overlap < 0.99) return 2;
        if (S.mps.overlap > topk + 0.005) return 3;
        return S.kind === "ground" && S.U > 4.5 ? 4 : 5;
    }

    function render() {
        S.mps = compress(S.psi, S.D);
        const v = S.mps.expand(), norm = Math.sqrt(dot(v, v));
        S.psiD = v.map((x) => x / norm);
        const topk = S.sq[Math.min(DIM, S.mps.stored) - 1];
        ui.dOut.textContent = String(S.D);
        ui.stored.textContent = fmtInt(S.mps.stored);
        ui.topkN.textContent = fmtInt(Math.min(DIM, S.mps.stored));
        ui.overlap.textContent = fmtOverlap(S.mps.overlap);
        ui.topk.textContent = fmtOverlap(S.mps.stored >= DIM ? 1 : topk);
        setNote(ui.stateNote, STATE_NOTES, stateNote());
        setNote(ui.mpsNote, mpsNotes(), mpsNote(topk));
        drawChain();
        drawAmps();
    }

    // ------------------------------------------------------------------
    // controls (work is done at most once per frame while sliders move)
    // ------------------------------------------------------------------
    let queued = false, dirtyState = true, started = false;
    function schedule(stateChanged) {
        dirtyState = dirtyState || stateChanged;
        if (queued || !started) return;
        queued = true;
        requestAnimationFrame(() => {
            queued = false;
            if (dirtyState) prepare();
            dirtyState = false;
            render();
        });
    }

    ui.kindSeg.addEventListener("click", (ev) => {
        const btn = ev.target.closest("button");
        if (!btn) return;
        S.kind = btn.dataset.kind;
        ui.kindSeg.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(b === btn)));
        ui.U.disabled = S.kind === "random";
        ui.uWrap.classList.toggle("tn-off", S.kind === "random");
        schedule(true);
    });
    ui.U.addEventListener("input", () => {
        S.U = Number(ui.U.value);
        ui.uOut.textContent = S.U.toFixed(1);
        schedule(true);
    });
    ui.D.addEventListener("input", () => {
        S.D = Number(ui.D.value);
        ui.dOut.textContent = String(S.D);
        schedule(false);
    });
    function showStep(k) {
        S.stepped = true;
        S.build = k;
        if (S.mps) drawChain();
    }
    ui.steps.addEventListener("click", (ev) => {
        const btn = ev.target.closest("button");
        if (btn) showStep(Number(btn.dataset.step));
    });
    // Next: the first press of a visitor starts from the beginning
    ui.next.addEventListener("click", () => showStep(!S.stepped || S.build === L - 1 ? 0 : S.build + 1));

    // the first ground state takes a moment: compute it once the widget is near the screen
    new IntersectionObserver((entries, obs) => {
        if (!entries[0].isIntersecting) return;
        obs.disconnect();
        started = true;
        schedule(true);
    }, { rootMargin: "300px" }).observe(root);

    new ResizeObserver(() => {
        fit(amp);
        if (S.psiD) drawAmps();
    }).observe(root);

    // light/dark toggle: Quarto swaps the stylesheet and the body's class
    new MutationObserver(() => requestAnimationFrame(() => {
        readColors();
        if (S.psiD) drawAmps();
    })).observe(document.body, { attributes: true, attributeFilter: ["class"] });

    ui.dOut.textContent = String(S.D);
    if (window.MathJax && window.MathJax.typesetPromise) window.MathJax.typesetPromise([root]);
})();
