// Gallery: the arcade cabinet (games), the Blender render window (3D models) and
// the Aseprite editor (pixel art). All are built by the listing templates in
// _templates/; this only switches which game / render / piece is showing.
(function () {
    const still = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

    // ---- the arcade cabinet ----
    document.querySelectorAll("[data-arcade]").forEach((cab) => {
        const slides = [...cab.querySelectorAll(".arcade-slide")];
        const signs = [...cab.querySelectorAll(".arcade-sign")];
        const abouts = [...cab.querySelectorAll(".arcade-about")];
        const picks = [...cab.querySelectorAll(".arcade-pick")];
        const go = cab.querySelector(".arcade-go");
        const count = cab.querySelector(".arcade-count");
        const pad = (n) => String(n).padStart(2, "0");
        let at = 0, attract = null, touched = false;

        function show(i) {
            at = (i + slides.length) % slides.length;
            [slides, signs, abouts, picks].forEach((list) => list.forEach((el, k) => el.classList.toggle("is-on", k === at)));
            go.href = abouts[at].dataset.href;
            go.setAttribute("aria-label", `Start ${abouts[at].dataset.title}`);
            count.textContent = `${pad(at + 1)} / ${pad(slides.length)}`;
        }

        // a cabinet button sinks for a moment when pressed
        function press(el) {
            el.classList.add("is-pressed");
            setTimeout(() => el.classList.remove("is-pressed"), 140);
        }

        // attract mode: the cabinet cycles through its games while it is on screen,
        // until someone plays with it
        function run(on) {
            clearInterval(attract);
            attract = on && !touched && !still ? setInterval(() => show(at + 1), 4500) : null;
        }
        function touch() { touched = true; run(false); }
        new IntersectionObserver(([e]) => run(e.isIntersecting)).observe(cab);

        cab.querySelectorAll(".arcade-hit[data-step]").forEach((b) =>
            b.addEventListener("click", () => { touch(); press(b); show(at + Number(b.dataset.step)); }));
        go.addEventListener("click", () => press(go));
        picks.forEach((p) => p.addEventListener("click", (e) => {
            e.preventDefault();
            touch();
            show(Number(p.dataset.index));
        }));
        cab.addEventListener("keydown", (e) => {
            if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
                e.preventDefault();
                touch();
                show(at + (e.key === "ArrowLeft" ? -1 : 1));
            }
        });
        cab.addEventListener("pointerdown", touch);
    });

    // ---- the Blender render window ----
    // Picking a model in the Outliner (or the arrows) shows its render.
    document.querySelectorAll("[data-blender]").forEach((win) => {
        const renders = [...win.querySelectorAll(".bl-render")];
        const objs = [...win.querySelectorAll(".bl-obj")];
        const props = [...win.querySelectorAll(".bl-props")];
        const name = win.querySelector(".bl-head-name");
        const count = win.querySelector(".bl-count");
        let at = 0;

        function show(i) {
            at = (i + renders.length) % renders.length;
            [renders, objs, props].forEach((list) => list.forEach((el, k) => el.classList.toggle("is-on", k === at)));
            name.textContent = objs[at].textContent.trim();
            count.textContent = `${at + 1} / ${renders.length}`;
        }

        objs.forEach((o) => o.addEventListener("click", () => show(Number(o.dataset.index))));
        win.querySelectorAll(".bl-step button").forEach((b) => b.addEventListener("click", () => show(at + Number(b.dataset.step))));
        win.addEventListener("keydown", (e) => {
            if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
                e.preventDefault();
                show(at + (e.key === "ArrowLeft" ? -1 : 1));
            }
        });
    });

    // ---- the Aseprite editor ----
    document.querySelectorAll("[data-ase]").forEach((ed) => {
        const arts = [...ed.querySelectorAll(".ase-art")];
        const props = [...ed.querySelectorAll(".ase-props")];
        const frames = [...ed.querySelectorAll(".ase-frame")];
        const file = ed.querySelector(".ase-file");
        const play = ed.querySelector(".ase-play");
        let at = 0, timer = null;

        function show(i) {
            at = (i + arts.length) % arts.length;
            [arts, props, frames].forEach((list) => list.forEach((el, k) => el.classList.toggle("is-on", k === at)));
            file.textContent = arts[at].dataset.file;
        }

        function stop() {
            clearInterval(timer);
            timer = null;
            play.innerHTML = "&#9654;";
            play.classList.remove("is-on");
        }

        frames.forEach((f) => f.addEventListener("click", () => { stop(); show(Number(f.dataset.index)); }));
        ed.querySelectorAll(".ase-step").forEach((b) =>
            b.addEventListener("click", () => { stop(); show(at + Number(b.dataset.step)); }));
        ed.addEventListener("keydown", (e) => {
            if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
                e.preventDefault();
                stop();
                show(at + (e.key === "ArrowLeft" ? -1 : 1));
            }
        });
        play.addEventListener("click", () => {
            if (timer) return stop();
            play.innerHTML = "&#10074;&#10074;";
            play.classList.add("is-on");
            show(at + 1);
            timer = setInterval(() => show(at + 1), 1400);
        });
    });
})();
