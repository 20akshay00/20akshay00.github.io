const canvas = document.getElementById('atom-canvas');
const ctx = canvas.getContext('2d');

let width, height, isMobile;
let atoms = [];
const numAtoms = 600;
const atomColor = "rgba(44, 62, 80, 0.175)";
const gridColor = "rgba(44, 62, 80, 0.075)";
const mouseRadius = 150;
const gridSize = 100;
const interactionRadius = 50;

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
            radius: Math.random() * 2 + 1
        });
    }
}

let mouse = { x: -1000, y: -1000 };
window.addEventListener('mousemove', (e) => {
    if (!isMobile) {
        mouse.x = e.clientX;
        mouse.y = e.clientY;
    }
});

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

function animate() {
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, width, height);

    drawGrid();

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

            if (distSq < interactionRadius * interactionRadius) {
                let dist = Math.sqrt(distSq);
                let force = (interactionRadius - dist) / interactionRadius;
                let angle = Math.atan2(dy, dx);
                let fx = Math.cos(angle) * force * 0.15;
                let fy = Math.sin(angle) * force * 0.15;

                a.vx -= fx;
                a.vy -= fy;
                b.vx += fx;
                b.vy += fy;
            }
        }

        if (!isMobile) {
            let dxMouse = mouse.x - a.x;
            let dyMouse = mouse.y - a.y;
            let distMouse = Math.sqrt(dxMouse * dxMouse + dyMouse * dyMouse);
            if (distMouse < mouseRadius) {
                let force = (mouseRadius - distMouse) / mouseRadius;
                let angle = Math.atan2(dyMouse, dxMouse);
                a.vx -= Math.cos(angle) * force * 0.6;
                a.vy -= Math.sin(angle) * force * 0.6;
            }
        } else {
            a.vx += (Math.random() - 0.5) * 0.1;
            a.vy += (Math.random() - 0.5) * 0.1;
        }

        a.vx *= 0.99;
        a.vy *= 0.99;

        a.x += a.vx;
        a.y += a.vy;

        ctx.fillStyle = atomColor;
        ctx.beginPath();
        ctx.arc(a.x, a.y, a.radius, 0, Math.PI * 2);
        ctx.fill();
    }

    requestAnimationFrame(animate);
}

window.addEventListener('resize', init);
init();
animate();