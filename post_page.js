// Tidbit and blog posts: put the title block on the page itself (the notebook page
// for tidbits, the airmail letter for blog posts; see .tidbit-post and .blog-post in
// styles.scss), so the page opens with its heading instead of starting mid-text.
(function () {
    const header = document.getElementById("title-block-header");
    const main = document.querySelector("main.content");
    if (!header || !main) return;
    // the banner layout classes position things on the page grid; on the page it is plain flow
    header.querySelectorAll(".page-columns, .page-full, .column-body").forEach((el) =>
        el.classList.remove("page-columns", "page-full", "column-body"));
    header.classList.remove("page-columns", "page-full");
    header.classList.add("post-head");
    // the letter's postmark carries the date
    const date = header.querySelector(".date");
    if (date) header.dataset.date = date.textContent.trim();
    main.prepend(header);
})();

// Figures that are dark ink on white paper (hand-drawn sketches, line plots) get the
// class "ink", which dark mode uses to invert them (styles.scss). An image counts as
// ink if its corners are near-white and only a small share of it is coloured, so
// photos, screenshots and heatmaps are left alone.
(function () {
    const SIZE = 48;
    const canvas = document.createElement("canvas");
    canvas.width = canvas.height = SIZE;
    const ctx = canvas.getContext("2d", { willReadFrequently: true });

    function classify(img) {
        try {
            ctx.clearRect(0, 0, SIZE, SIZE);
            ctx.drawImage(img, 0, 0, SIZE, SIZE);
            const d = ctx.getImageData(0, 0, SIZE, SIZE).data;
            const at = (x, y) => 4 * (y * SIZE + x);
            const white = (i) => d[i + 3] > 200 && d[i] > 225 && d[i + 1] > 225 && d[i + 2] > 225;
            const corners = [at(1, 1), at(SIZE - 2, 1), at(1, SIZE - 2), at(SIZE - 2, SIZE - 2)];
            if (!corners.every(white)) return;
            let coloured = 0, light = 0;
            for (let i = 0; i < d.length; i += 4) {
                const max = Math.max(d[i], d[i + 1], d[i + 2]), min = Math.min(d[i], d[i + 1], d[i + 2]);
                if (max - min > 60) coloured++;
                if (min > 200) light++;
            }
            const n = SIZE * SIZE;
            if (coloured / n < 0.15 && light / n > 0.6) img.classList.add("ink");
        } catch (e) {
            // an image from another site cannot be read back: leave it as it is
        }
    }

    document.querySelectorAll("main.content figure img, main.content .cell-output-display img").forEach((img) => {
        if (img.complete && img.naturalWidth) classify(img);
        else img.addEventListener("load", () => classify(img), { once: true });
    });
})();

// Inline maths never line-breaks, so on a narrow screen a long inline formula can be
// wider than the page and drag the whole page sideways. Any that is wider than its
// paragraph gets a line of its own that scrolls sideways (.math-too-wide in
// theme/_tidbits.scss); the rest stay inline. Checked once MathJax has typeset, and
// again when the window is resized.
(function () {
    function fit() {
        document.querySelectorAll("main.content span.math.inline").forEach((m) => {
            m.classList.remove("math-too-wide");
            const column = m.closest("p, li, td") || m.parentElement;
            if (m.getBoundingClientRect().width > column.clientWidth) m.classList.add("math-too-wide");
        });
    }
    function whenTypeset(run) {
        const mj = window.MathJax;
        if (mj && mj.startup && mj.startup.promise) mj.startup.promise.then(run);
        else window.addEventListener("load", run);
    }
    whenTypeset(fit);
    let pending = 0;
    window.addEventListener("resize", () => { clearTimeout(pending); pending = setTimeout(fit, 150); });
})();
