const canvas = document.getElementById('atom-canvas');
const ctx = canvas.getContext('2d');

let width, height, isMobile;
let atoms = [];
// The number of atoms follows the window's area (600 on a typical laptop screen),
// so a big monitor is as dense as a small one; capped, since every pair interacts.
const atomDensity = 600 / (1440 * 900);
const minAtoms = 250, maxAtoms = 1000;

// Atoms fade from slate (hot gas) to purple (condensate) as they cool. One palette
// per site theme (light, dark), picked from the class Quarto puts on the body.
// The canvas is transparent: the page's notebook paper shows through it.
const PALETTES = {
    light: { hot: [42, 35, 56, 0.17], cold: [90, 42, 156, 0.45] },
    dark: { hot: [226, 216, 240, 0.16], cold: [180, 140, 240, 0.5] },
};
const colorBins = 6;
let binColors;
function readPalette() {
    const p = document.body.classList.contains('quarto-dark') ? PALETTES.dark : PALETTES.light;
    binColors = Array.from({ length: colorBins }, (_, k) => {
        const cold = k / (colorBins - 1);
        const c = p.hot.map((h, i) => h + (p.cold[i] - h) * cold);
        return `rgba(${c[0] | 0}, ${c[1] | 0}, ${c[2] | 0}, ${c[3].toFixed(3)})`;
    });
}
readPalette();
new MutationObserver(() => {
    readPalette();
    if (reducedMotion) draw();   // the still frame is not redrawn otherwise
}).observe(document.body, { attributes: true, attributeFilter: ['class'] });

// With "reduce motion" set in the OS, the atoms are drawn once and stay put.
const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

const interactionRadius = 50;
const coreRadius = 9;

// --- the cursor ---
// TRAP_ON = false: the cursor simply pushes atoms away (radius 150 px, as before).
// TRAP_ON = true: the cursor is an optical trap. Held still, it cools nearby atoms (molasses friction) and gently confines
// them, so they settle into a single cloud. The cloud keeps a finite size
// because cold atoms still repel softly (like a Thomas-Fermi condensate).
// Moved quickly, it stirs and heats them instead.
const TRAP_ON = false;
const pushRadius = 150;
const trapRadius = 170;
const trapStiffness = 0.0006;
const molasses = 0.08;      // max velocity damping per frame at the centre
const trapCooling = 0.04;   // max temperature drop per frame at the centre
const ambientRate = 0.004;  // how fast atoms relax back to the hot gas (~4s)
const zeroPoint = 0.03;     // residual jiggle at T = 0
const coldRepulsion = 0.015; // interaction strength in the condensate: higher = bigger, looser cloud

let mouse = { x: -1000, y: -1000, speed: 0, active: false };

function init() {
    width = canvas.width = window.innerWidth;
    height = canvas.height = window.innerHeight;
    isMobile = width < 768;

    const currentNumAtoms = Math.round(Math.min(maxAtoms, Math.max(minAtoms, width * height * atomDensity)));

    atoms = [];
    for (let i = 0; i < currentNumAtoms; i++) {
        atoms.push({
            x: Math.random() * width,
            y: Math.random() * height,
            vx: (Math.random() - 0.5) * 2,
            vy: (Math.random() - 0.5) * 2,
            radius: Math.random() * 2 + 1,
            T: 1
        });
    }
}

function moveTrap(x, y) {
    if (mouse.active) {
        mouse.speed = Math.max(mouse.speed, Math.hypot(x - mouse.x, y - mouse.y));
    }
    mouse.x = x;
    mouse.y = y;
    mouse.active = true;
}

window.addEventListener('mousemove', (e) => moveTrap(e.clientX, e.clientY));
document.addEventListener('mouseleave', () => { mouse.active = false; mouse.x = mouse.y = -1000; });

window.addEventListener('touchstart', (e) => {
    mouse.active = false;
    moveTrap(e.touches[0].clientX, e.touches[0].clientY);
}, { passive: true });
window.addEventListener('touchmove', (e) => moveTrap(e.touches[0].clientX, e.touches[0].clientY), { passive: true });
window.addEventListener('touchend', () => { mouse.active = false; mouse.x = mouse.y = -1000; });

// The physics advances in fixed steps of 1/60 s, however often the screen refreshes:
// on a 144 Hz monitor a step per frame would run everything 2.4 times as fast.
const STEP = 1000 / 60;
let lastTime = 0, pending = 0;

// --- neighbour grid ---
// Checking every pair of atoms costs n^2 (half a million pairs per step for 1000
// atoms). Instead the atoms are sorted into square cells a bit wider than the
// interaction radius (the margin covers how far atoms move during a step), and
// each atom only looks at its own cell and the eight around it.
const CELL = interactionRadius + 10;
let gridCols = 0, gridRows = 0, cellHead = new Int32Array(0), cellNext = new Int32Array(0),
    cellX = new Int32Array(0), cellY = new Int32Array(0);

function buildGrid() {
    const n = atoms.length;
    gridCols = Math.max(1, Math.ceil(width / CELL));
    gridRows = Math.max(1, Math.ceil(height / CELL));
    if (cellHead.length !== gridCols * gridRows) cellHead = new Int32Array(gridCols * gridRows);
    if (cellNext.length !== n) {
        cellNext = new Int32Array(n);
        cellX = new Int32Array(n);
        cellY = new Int32Array(n);
    }
    cellHead.fill(-1);
    for (let i = 0; i < n; i++) {
        // atoms that drift past the edge share the edge cells
        const gx = Math.min(gridCols - 1, Math.max(0, Math.floor(atoms[i].x / CELL)));
        const gy = Math.min(gridRows - 1, Math.max(0, Math.floor(atoms[i].y / CELL)));
        cellX[i] = gx;
        cellY[i] = gy;
        const c = gy * gridCols + gx;
        cellNext[i] = cellHead[c];
        cellHead[c] = i;
    }
}

function animate(now) {
    pending += lastTime ? Math.min(now - lastTime, 100) : STEP;
    lastTime = now;
    while (pending >= STEP) {
        step();
        pending -= STEP;
    }
    draw();
    requestAnimationFrame(animate);
}

function step() {
    // 1 when the cursor rests, 0 when it is moving fast
    const stillness = Math.max(0, 1 - mouse.speed / 6);
    const stir = Math.min(1, mouse.speed / 25);
    mouse.speed *= 0.85;

    const R2 = interactionRadius * interactionRadius;
    buildGrid();

    for (let i = 0; i < atoms.length; i++) {
        let a = atoms[i];

        if (a.x < 50) a.vx += 0.05;
        if (a.x > width - 50) a.vx -= 0.05;
        if (a.y < 50) a.vy += 0.05;
        if (a.y > height - 50) a.vy -= 0.05;

        // only atoms in this cell and the eight around it can be close enough to interact;
        // each pair is handled once, by its lower index, as in a plain double loop
        const ci = cellX[i], cj = cellY[i];
        for (let gy = Math.max(0, cj - 1); gy <= Math.min(gridRows - 1, cj + 1); gy++)
        for (let gx = Math.max(0, ci - 1); gx <= Math.min(gridCols - 1, ci + 1); gx++)
        for (let j = cellHead[gy * gridCols + gx]; j !== -1; j = cellNext[j]) {
            if (j <= i) continue;
            let b = atoms[j];
            let dx = b.x - a.x;
            let dy = b.y - a.y;
            let distSq = dx * dx + dy * dy;

            if (distSq < R2 && distSq > 0) {
                let dist = Math.sqrt(distSq);
                let T = 0.5 * (a.T + b.T);
                // soft repulsion: strong in the hot gas, weaker (but never zero)
                // when cold, which sets the size of the condensed cloud
                let force = (0.15 * T + coldRepulsion * (1 - T)) * (interactionRadius - dist) / interactionRadius;
                // hard core keeps atoms from sitting on top of each other
                if (dist < coreRadius) force += 0.3 * (coreRadius - dist) / coreRadius;

                let fx = (dx / dist) * force;
                let fy = (dy / dist) * force;

                a.vx -= fx;
                a.vy -= fy;
                b.vx += fx;
                b.vy += fy;
            }
        }

        if (mouse.active && !TRAP_ON) {
            let dx = mouse.x - a.x;
            let dy = mouse.y - a.y;
            let dist = Math.sqrt(dx * dx + dy * dy);
            if (dist < pushRadius && dist > 0) {
                let force = (pushRadius - dist) / pushRadius;
                a.vx -= (dx / dist) * force * 0.6;
                a.vy -= (dy / dist) * force * 0.6;
            }
        } else if (mouse.active) {
            let dx = mouse.x - a.x;
            let dy = mouse.y - a.y;
            let dist = Math.sqrt(dx * dx + dy * dy);
            if (dist < trapRadius && dist > 0) {
                let w = (trapRadius - dist) / trapRadius;
                // gentle confinement + molasses friction + cooling while the cursor rests
                a.vx += dx * trapStiffness * stillness;
                a.vy += dy * trapStiffness * stillness;
                let friction = 1 - molasses * w * stillness;
                a.vx *= friction;
                a.vy *= friction;
                a.T -= trapCooling * w * stillness;
                // stirring + heating while it moves
                a.vx -= (dx / dist) * w * 0.6 * stir;
                a.vy -= (dy / dist) * w * 0.6 * stir;
                a.T += 0.05 * w * stir;
            }
        }

        // relax back towards the hot ambient gas
        a.T += (1 - a.T) * ambientRate;
        a.T = Math.min(1, Math.max(0, a.T));

        // thermal kicks
        let noise = 0.25 * a.T + zeroPoint;
        a.vx += (Math.random() - 0.5) * noise;
        a.vy += (Math.random() - 0.5) * noise;

        let damping = 0.99 - 0.04 * (1 - a.T);
        a.vx *= damping;
        a.vy *= damping;

        a.x += a.vx;
        a.y += a.vy;
    }
}

function draw() {
    ctx.clearRect(0, 0, width, height);
    // draw in a few colour bins so each bin is a single path
    for (let k = 0; k < colorBins; k++) {
        ctx.fillStyle = binColors[k];
        ctx.beginPath();
        for (const a of atoms) {
            if (Math.round((1 - a.T) * (colorBins - 1)) !== k) continue;
            ctx.moveTo(a.x + a.radius, a.y);
            ctx.arc(a.x, a.y, a.radius, 0, Math.PI * 2);
        }
        ctx.fill();
    }
}


let lastWidth = window.innerWidth;

window.addEventListener('resize', () => {
    const currentWidth = window.innerWidth;
    const currentHeight = window.innerHeight;

    if (isMobile) {
        if (Math.abs(currentWidth - lastWidth) > 50) {
            lastWidth = currentWidth;
            init();
        } else {
            width = canvas.width = currentWidth;
            height = canvas.height = currentHeight;
        }
    } else {
        lastWidth = currentWidth;
        init();
    }
    if (reducedMotion) draw();
});

init();
if (reducedMotion) {
    // a still scatter of atoms: the look of the page without the motion
    draw();
} else {
    requestAnimationFrame(animate);
}
