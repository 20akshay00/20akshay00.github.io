const canvas = document.getElementById('atom-canvas');
const ctx = canvas.getContext('2d');

let width, height, isMobile;
let atoms = [];
const numAtoms = 600;

const bgColor = '#fcfcf9';
const gridColor = "rgba(44, 62, 80, 0.075)";
// atoms fade from slate (hot gas) to teal (condensate) as they cool
const hotColor = [44, 62, 80, 0.175];
const coldColor = [0, 109, 119, 0.45];
const colorBins = 6;

const gridSize = 100;
const interactionRadius = 50;
const coreRadius = 9;

// --- the cursor is an optical trap ---
// Held still, it cools nearby atoms (molasses friction) and gently confines
// them, so they settle into a single cloud. The cloud keeps a finite size
// because cold atoms still repel softly (like a Thomas-Fermi condensate).
// Moved quickly, it stirs and heats them instead.
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

    const currentNumAtoms = isMobile ? 250 : numAtoms;

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

function drawGrid(angle = Math.PI / 180 * 30) {
    const diag = Math.sqrt(width ** 2 + height ** 2);
    ctx.save();
    ctx.translate(width / 2, height / 2);
    ctx.rotate(angle);
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let x = -diag / 2; x <= diag / 2; x += gridSize) {
        ctx.moveTo(x, -diag / 2);
        ctx.lineTo(x, diag / 2);
    }
    for (let y = -diag / 2; y <= diag / 2; y += gridSize) {
        ctx.moveTo(-diag / 2, y);
        ctx.lineTo(diag / 2, y);
    }
    ctx.stroke();
    ctx.restore();
}

const binColors = Array.from({ length: colorBins }, (_, k) => {
    const cold = k / (colorBins - 1);
    const c = hotColor.map((h, i) => h + (coldColor[i] - h) * cold);
    return `rgba(${c[0] | 0}, ${c[1] | 0}, ${c[2] | 0}, ${c[3].toFixed(3)})`;
});

function animate() {
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, width, height);

    drawGrid();

    // 1 when the cursor rests, 0 when it is moving fast
    const stillness = Math.max(0, 1 - mouse.speed / 6);
    const stir = Math.min(1, mouse.speed / 25);
    mouse.speed *= 0.85;

    const R2 = interactionRadius * interactionRadius;

    for (let i = 0; i < atoms.length; i++) {
        let a = atoms[i];

        if (a.x < 50) a.vx += 0.05;
        if (a.x > width - 50) a.vx -= 0.05;
        if (a.y < 50) a.vy += 0.05;
        if (a.y > height - 50) a.vy -= 0.05;

        for (let j = i + 1; j < atoms.length; j++) {
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

        if (mouse.active) {
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

    requestAnimationFrame(animate);
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
});

init();
animate();
