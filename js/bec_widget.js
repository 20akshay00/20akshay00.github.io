// Raw absorption images of an oscillating condensate, played as a movie over
// the hold time: either the mean of 20 runs or, for each frame, one run picked
// at random.
//
// data/bec_shots.bin (deflated bytes): for each of the 15 hold times the mean of all
// 20 shots, then SHOTS single shots (shot-major), each frame H x W pixels of
// optical density, scaled so that 255 is the 99.9th percentile of the shots.
(function () {
    const root = document.getElementById("bec-widget");
    if (!root) return;
    const SRC = new URL("../data/bec_shots.bin", document.currentScript.src).href;

    const HOLD = [0, 0.2, 0.3, 0.5, 0.6, 0.8, 1.0, 1.1, 1.3, 1.5, 1.6, 1.8, 1.9, 2.1, 2.3];   // ms
    const T = HOLD.length, H = 96, W = 47, F = H * W, SHOTS = 5;
    const PX = 10.33;                       // oscillator lengths per pixel (both axes)
    const PERIOD = 1000 / 440;              // trap period, ms (fit to the cloud centre: 441.5 Hz)
    const FRAME_MS = 260;                   // time per hold time while playing

    // viridis, as in the lab's own plots
    const VIRIDIS = ["#440154", "#472d7b", "#3b528b", "#2c728e", "#21918c", "#28ae80", "#5ec962", "#addc30", "#fde725"]
        .map((h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)));
    const LUT = new Uint8ClampedArray(256 * 3);
    for (let v = 0; v < 256; v++) {
        const s = (v / 255) * (VIRIDIS.length - 1), i = Math.min(VIRIDIS.length - 2, Math.floor(s)), f = s - i;
        for (let c = 0; c < 3; c++) LUT[3 * v + c] = VIRIDIS[i][c] + f * (VIRIDIS[i + 1][c] - VIRIDIS[i][c]);
    }

    root.classList.add("tdse", "bec");
    root.innerHTML = `
        <figure>
            <canvas class="bec-frame" aria-label="Absorption image of the cloud at the current hold time"></canvas>
            <div class="bec-controls">
                <button type="button" class="tdse-btn" data-ui="play"></button>
                <div class="tdse-seg" role="group" aria-label="Which images" data-ui="mode">
                    <button type="button" data-mode="mean" aria-pressed="true">Mean</button>
                    <button type="button" data-mode="shot" aria-pressed="false">Single shot</button>
                </div>
            </div>
            <figcaption>Absorption images of an oscillating condensate after 12 ms of free flight (normalised optical density). Each frame is a new run.</figcaption>
        </figure>`;

    const ui = {};
    root.querySelectorAll("[data-ui]").forEach((el) => (ui[el.dataset.ui] = el));
    const cv = root.querySelector(".bec-frame"), cx = cv.getContext("2d");
    // axis text colour, following the light/dark theme
    const readText = () => getComputedStyle(root).getPropertyValue("--tdse-text").trim() || "#5b6b78";
    let TEXT = readText();

    // ------------------------------------------------------------------
    // state: t the hold time shown, shot -1 for the mean or the run shown
    // ------------------------------------------------------------------
    const S = { t: 0, mode: "mean", shot: -1, playing: !window.matchMedia("(prefers-reduced-motion: reduce)").matches };
    let tiles = null;                       // one small canvas per frame, in viridis, x horizontal

    function makeTiles(bytes) {
        tiles = [];
        for (let k = 0; k < (1 + SHOTS) * T; k++) {
            const tile = document.createElement("canvas");
            tile.width = H;
            tile.height = W;
            const img = tile.getContext("2d").createImageData(H, W);
            // turned on its side: data row r (x, most negative first) runs left to right
            for (let r = 0; r < H; r++) {
                for (let c = 0; c < W; c++) {
                    const v = bytes[k * F + r * W + c];
                    img.data.set([LUT[3 * v], LUT[3 * v + 1], LUT[3 * v + 2], 255], 4 * (c * H + r));
                }
            }
            tile.getContext("2d").putImageData(img, 0, 0);
            tiles.push(tile);
        }
    }

    // a new run for every frame in single-shot mode
    function pickShot() {
        S.shot = S.mode === "mean" ? -1 : Math.floor(Math.random() * SHOTS);
    }

    // ------------------------------------------------------------------
    // drawing: the hold time, the frame with its x axis, and a colour bar
    // ------------------------------------------------------------------
    const LABEL = 18, AXIS = 34, BAR = 9, GAP = 8, BAR_TEXT = 14;
    let w = 0, h = 0;
    function fit() {
        const dpr = window.devicePixelRatio || 1;
        w = cv.clientWidth;
        h = LABEL + ((w - GAP - BAR - BAR_TEXT) * W) / H + AXIS;   // as tall as the frame needs
        cv.style.height = `${Math.round(h)}px`;
        cv.width = Math.round(w * dpr);
        cv.height = Math.round(h * dpr);
        cx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function draw() {
        cx.clearRect(0, 0, w, h);
        if (!tiles) return;
        const fx = 0, fy = LABEL, fw = w - GAP - BAR - BAR_TEXT, fh = (fw * W) / H;
        cx.imageSmoothingEnabled = false;
        cx.drawImage(tiles[S.shot < 0 ? S.t : T * (1 + S.shot) + S.t], fx, fy, fw, fh);

        cx.font = "12px 'Source Sans 3', sans-serif";
        cx.fillStyle = TEXT;
        cx.strokeStyle = TEXT;
        cx.lineWidth = 1;

        // hold time above the frame
        cx.textAlign = "left";
        cx.textBaseline = "bottom";
        cx.fillText(`Hold time: ${(HOLD[S.t] / PERIOD).toFixed(2)} trap periods`, fx, fy - 4);

        // x axis: the centre of the frame is x = 0
        cx.textAlign = "center";
        cx.textBaseline = "top";
        for (const xv of [-400, 0, 400]) {
            const x = Math.round(fx + (H / 2 + xv / PX) * (fw / H)) + 0.5;
            cx.beginPath();
            cx.moveTo(x, fy + fh);
            cx.lineTo(x, fy + fh + 4);
            cx.stroke();
            cx.fillText(String(xv), x, fy + fh + 5);
        }
        cx.fillText("x (ℓₕ)", fx + fw / 2, fy + fh + 19);

        // colour bar, optical density from 0 to 1
        const bx = fx + fw + GAP;
        for (let i = 0; i < fh; i++) {
            const v = Math.round(255 * (1 - i / fh));
            cx.fillStyle = `rgb(${LUT[3 * v]}, ${LUT[3 * v + 1]}, ${LUT[3 * v + 2]})`;
            cx.fillRect(bx, fy + i, BAR, 1.5);
        }
        cx.fillStyle = TEXT;
        cx.textAlign = "left";
        cx.textBaseline = "top";
        cx.fillText("1", bx + BAR + 3, fy);
        cx.textBaseline = "bottom";
        cx.fillText("0", bx + BAR + 3, fy + fh);
    }

    // ------------------------------------------------------------------
    // controls
    // ------------------------------------------------------------------
    function updatePlay() {
        ui.play.innerHTML = S.playing
            ? '<i class="bi bi-pause-fill" aria-hidden="true"></i><span>Pause</span>'
            : '<i class="bi bi-play-fill" aria-hidden="true"></i><span>Play</span>';
    }
    ui.play.addEventListener("click", () => {
        S.playing = !S.playing;
        updatePlay();
        tick();
    });
    ui.mode.addEventListener("click", (ev) => {
        const btn = ev.target.closest("button");
        if (!btn) return;
        S.mode = btn.dataset.mode;
        ui.mode.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(b === btn)));
        pickShot();
        draw();
    });

    // step through the hold times while playing and on screen
    let visible = false, timer = null;
    function tick() {
        clearTimeout(timer);
        if (!S.playing || !visible || document.hidden || !tiles) return;
        timer = setTimeout(() => {
            S.t = (S.t + 1) % T;
            pickShot();
            draw();
            tick();
        }, FRAME_MS);
    }

    // ------------------------------------------------------------------
    // loading: fetch the images once the section is near the screen
    // ------------------------------------------------------------------
    async function load() {
        try {
            const res = await fetch(SRC);
            if (!res.ok) throw new Error(res.status);
            const bytes = new Uint8Array(await new Response(res.body.pipeThrough(new DecompressionStream("deflate"))).arrayBuffer());
            if (bytes.length !== (1 + SHOTS) * T * F) throw new Error("unexpected size");
            makeTiles(bytes);
            draw();
            tick();
        } catch (err) {
            root.querySelector("figcaption").textContent = "The images could not be loaded.";
        }
    }

    let loading = false;
    new IntersectionObserver((entries) => {
        visible = entries[0].isIntersecting;
        if (visible && !loading) {
            loading = true;
            load();
        }
        tick();
    }, { rootMargin: "200px" }).observe(root);
    document.addEventListener("visibilitychange", tick);

    new ResizeObserver(() => {
        fit();
        draw();
    }).observe(cv);

    // light/dark toggle: Quarto swaps the stylesheet and the body's class
    new MutationObserver(() => requestAnimationFrame(() => {
        TEXT = readText();
        draw();
    })).observe(document.body, { attributes: true, attributeFilter: ["class"] });

    updatePlay();
})();
