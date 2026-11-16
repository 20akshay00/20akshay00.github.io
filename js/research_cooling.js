// Research page: the page cools down as you read.
//
// The scroll position sets a temperature, from 1 uK at the top of the page to
// 10 nK at the bottom (on a log scale). The page background shifts from a faint
// warm tint to a cool one, and a thermometer in the corner doubles as a scrollbar:
// it shows T and, below Tc, the ideal-gas condensate fraction N0/N = 1 - (T/Tc)^3.
// Click or drag on it to move through the page.
(function () {
    const T_TOP = 1000, T_BOTTOM = 10, T_C = 150;   // nK
    const PALETTE = {
        light: { warm: [246, 233, 240], cold: [230, 234, 248] },
        dark: { warm: [33, 19, 29], cold: [17, 21, 40] },
    };

    const meter = document.createElement("div");
    meter.className = "gas-meter";
    meter.innerHTML = '<div class="gas-tube" role="scrollbar" aria-orientation="vertical" aria-label="Page position, as a temperature" tabindex="-1"><div class="gas-fill"></div><div class="gas-tc"></div></div><div class="gas-read"></div>';
    document.body.appendChild(meter);
    const tube = meter.querySelector(".gas-tube"), fill = meter.querySelector(".gas-fill"),
        tcMark = meter.querySelector(".gas-tc"), read = meter.querySelector(".gas-read");

    const maxScroll = () => Math.max(1, document.documentElement.scrollHeight - innerHeight);
    const level = (T) => Math.log(T / T_BOTTOM) / Math.log(T_TOP / T_BOTTOM);   // 1 at the top, 0 at the bottom
    const mix = (a, b, p) => a.map((v, i) => Math.round(v + (b[i] - v) * p));

    function update() {
        const p = Math.min(1, Math.max(0, scrollY / maxScroll()));
        const T = T_TOP * Math.pow(T_BOTTOM / T_TOP, p);
        const f = T < T_C ? 1 - Math.pow(T / T_C, 3) : 0;
        const pal = document.body.classList.contains("quarto-dark") ? PALETTE.dark : PALETTE.light;
        document.body.style.backgroundColor = `rgb(${mix(pal.warm, pal.cold, p).join(", ")})`;
        fill.style.height = `${100 * level(T)}%`;
        tcMark.style.bottom = `${100 * level(T_C)}%`;
        meter.classList.toggle("bec", f > 0);
        tube.setAttribute("aria-valuenow", Math.round(100 * p));
        read.innerHTML = `T = ${Math.round(T)} nK<br>${f > 0 ? `N<sub>0</sub>/N = ${Math.round(100 * f)}%` : "thermal gas"}`;
    }

    // the thermometer as a scrollbar: hot (top of the tube) is the top of the page
    function scrollFromPointer(e) {
        const r = tube.getBoundingClientRect();
        const p = 1 - Math.min(1, Math.max(0, (e.clientY - r.top) / r.height));
        scrollTo({ top: p * maxScroll(), behavior: "instant" });
    }
    tube.addEventListener("pointerdown", (e) => {
        tube.setPointerCapture(e.pointerId);
        scrollFromPointer(e);
    });
    tube.addEventListener("pointermove", (e) => tube.hasPointerCapture(e.pointerId) && scrollFromPointer(e));

    // scroll events can fire many times per frame: update at most once per frame
    let queued = false;
    const schedule = () => {
        if (queued) return;
        queued = true;
        requestAnimationFrame(() => { queued = false; update(); });
    };
    addEventListener("scroll", schedule, { passive: true });
    addEventListener("resize", schedule);
    new MutationObserver(update).observe(document.body, { attributes: true, attributeFilter: ["class"] });
    update();
})();
