// The robot on the landing page: it strolls along the divider above the links,
// stops at each end for a moment (its idle animation), turns round and walks back.
// Sprite sheets: 32 x 32 frames side by side, drawn at 2x; it faces right.
// The script sits inside the page text, before the divider exists: wait for the page.
document.addEventListener("DOMContentLoaded", function () {
    const about = document.querySelector(".quarto-about-jolla");
    const line = about && about.querySelector(".about-sep");
    if (!line) return;

    const SCALE = 2, SPEED = 38, PAUSE = 1.4;   // px per s, s at each end
    const FEET = 26;                             // its feet are on row 25 of 32

    // load and decode both sheets now: the idle one is otherwise first fetched when the
    // robot reaches the end, and it blinks out while that loads. Holding the images
    // keeps them decoded, so later swaps are instant too.
    const sheets = ["robot_walk", "robot_idle"].map((n) => {
        const img = new Image();
        img.src = `/images/sprites/${n}.png`;
        img.decode().catch(() => {});
        return img;
    });

    const robot = document.createElement("div");
    robot.sheets = sheets;
    robot.className = "robot-walker";
    robot.setAttribute("aria-hidden", "true");
    about.style.position = "relative";
    about.appendChild(robot);

    const still = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    let x = 0, dir = 1, rest = 0, t = 0, last = 0, visible = false, raf = 0;

    // the divider's span, in the about box's coordinates
    function span() {
        const a = about.getBoundingClientRect(), r = line.getBoundingClientRect();
        return { left: r.left - a.left, right: r.right - a.left - 32 * SCALE, top: r.top - a.top };
    }

    function step(dt) {
        const s = span();
        t += dt;
        if (rest > 0) {
            rest -= dt;
        } else if (!still) {
            x += dir * SPEED * dt;
            if (x > s.right - s.left || x < 0) {
                x = Math.max(0, Math.min(s.right - s.left, x));
                dir = -dir;
                rest = PAUSE;
            }
        }
        const walking = rest <= 0 && !still;
        const sheet = walking ? "walk" : "idle";
        const frame = walking ? Math.floor(t * 11) % 10 : Math.floor(t * 2) % 2;
        robot.dataset.sheet = sheet;
        robot.style.backgroundPosition = `${-frame * 32 * SCALE}px 0`;
        robot.style.left = `${s.left + x}px`;
        robot.style.top = `${s.top - FEET * SCALE}px`;
        // face the way it is going (or, while resting, the way it will go next)
        robot.style.transform = `scaleX(${dir})`;
    }

    function loop(now) {
        step(Math.min(0.05, (now - last) / 1000 || 0));
        last = now;
        raf = visible ? requestAnimationFrame(loop) : 0;
    }

    new IntersectionObserver((entries) => {
        visible = entries[0].isIntersecting;
        if (visible && !raf) {
            last = performance.now();
            raf = requestAnimationFrame(loop);
        }
    }).observe(line);
    step(0);
});
