// Guided demos for the interactive Schrodinger / Gross-Pitaevskii widget
// (tdse_widget.js). Each demo snaps the widget to a setup and shows a caption.
// This file is data only: edit the text freely. Captions can be wrapped over
// several lines; line breaks and indentation are collapsed into single spaces.
//
// Fields of a demo:
//   label     the chip's name
//   topic     which tab it is under (one of the keys in `topics`)
//   refs      [text, link] pairs for the experiments it follows (optional)
//   preset    the potential (a key of `presets` in tdse_widget.js)
//   p1, p2    the potential's two sliders; g the interaction strength
//   speed     simulation steps per frame (default 10); imaginary: start in imaginary time
//   prep      the initial state: P({...}) a wavepacket, G({...}) a ground state,
//             or { kind: "leftwell" }; the keys match the initial-state panel
//   show      which displays to switch on (showP momentum, showPh phase,
//             showHist space-time, showRe/showIm, showE energy, showCl classical)
//   hist      simulated time per row of the space-time strip (default 0.5)
//   histMode  "change" shows the density minus its value at t = 0
//   caption   what to look for
//
// Several captions quote numbers the simulation produces (percentages, widths,
// periods): if you change a demo's parameters, those numbers change too.
window.TDSE_DEMOS = (() => {
    // the tabs: [key, title]
    const topics = [
        ["tunnel", "Tunnelling"],
        ["waves", "Waves"],
        ["lattice", "Lattices"],
        ["bec", "Condensates"],
        ["solitons", "Solitons"],
        ["superfluid", "Superfluidity"],
        ["driven", "Driven"],
    ];

    // a wavepacket / a ground state, with the usual defaults
    const P = (o) => ({ kind: "packet", shape: "gaussian", ...o });
    const G = (o) => ({ kind: "ground", prepIn: "same", imprint: "none", sameG: true, ...o });

    const demos = [
        // ---- tunnelling
        {
            label: "Tunnelling", topic: "tunnel", preset: "barrier", p1: 2.4, p2: 0.8,
            prep: P({ x0: -15, k: 1.8, sigma: 3 }), show: { showHist: true }, hist: 0.2,
            caption: `
                When the wavepacket has less energy than the barrier, the classical particle always
                bounces back, but the quantum particle has a probability of tunneling through.
            `,
        },

        {
            label: "Resonant tunnelling", topic: "tunnel", preset: "double", p1: 4, p2: 3,
            prep: P({ x0: -25, k: 2.43, sigma: 7 }), show: { showHist: true }, hist: 0.3,
            caption: `
                At this energy, one of these barriers on its own lets through about 44%. The pair
                lets through about 70%: waves bouncing between the barriers interfere and cancel
                most of the reflection. Change the momentum to 1.95 and only about 5% gets through.
            `,
        },

        {
            label: "Reflection off a cliff", topic: "tunnel", preset: "step", p1: 3, speed: 20,
            prep: P({ x0: 18, k: -0.5, sigma: 5 }), show: { showHist: true }, hist: 0.4,
            caption: `
                The packet runs down a sudden drop. A classical particle would just speed up and
                carry on, but almost half of the quantum probability density bounces back off the edge of the cliff.
            `,
        },

        {
            label: "Ramsauer–Townsend effect", topic: "tunnel", preset: "well", p1: 2, p2: 3, speed: 20,
            prep: P({ x0: -24, k: 0.85, sigma: 7 }), show: { showHist: true }, hist: 0.5,
            caption: `
                A slow particle passing over a dip. Classically it just speeds up and carries on,
                but a quantum particle is partly reflected by the dip. At this energy the wave fits the
                well exactly and 99% passes through. Change the momentum to 0.3 and only about 37%
                does.
            `,
        },

        {
            label: "Double-well tunnelling", topic: "tunnel", preset: "doublewell", p1: 0.3, p2: 8, speed: 30,
            refs: [["Albiez et al., PRL 2005", "https://arxiv.org/abs/cond-mat/0411757"]],
            prep: { kind: "leftwell" }, show: { showHist: true }, hist: 2,
            caption: `
                A particle starting in one well tunnels to the other and back again, even though it
                never has enough energy to climb over the barrier.
            `,
        },

        // ---- waves
        {
            label: "Spreading", topic: "waves", preset: "free",
            prep: P({ x0: -12, k: 1, sigma: 0.7 }), show: { showHist: true, showP: true }, hist: 0.2,
            caption: `
                A tightly localised packet contains a wide range of momenta, so it spreads out as it
                moves: the more precisely a particle is placed, the faster it disperses. The
                momentum distribution below never changes: spreading is just the different momenta
                running apart.
            `,
        },

        {
            label: "Interference", topic: "waves", preset: "free",
            prep: P({ shape: "pair", x0: 0, k: 2, sigma: 1.5, d: 20, phi: 0 }), show: { showHist: true }, hist: 0.1,
            caption: `
                Two packets pass through each other and form fringes spaced by π/k, as in
                matter-wave interference experiments with condensates.
            `,
        },

        {
            label: "Talbot effect", topic: "waves", preset: "pulse", p1: 1, speed: 5,
            refs: [["Chapman et al., PRA 1995", "https://doi.org/10.1103/PhysRevA.51.R14"], ["Deng et al., PRL 1999", "https://doi.org/10.1103/PhysRevLett.83.5407"]],
            prep: P({ x0: 0, k: 0, sigma: 15 }), show: { showHist: true, showE: false }, hist: 0.2,
            caption: `
                A light grating flashed on a wide, resting wave only ripples its phase, which the
                density cannot see. As the wave evolves, stripes appear, sharpest at t = π, fade to
                flat at 2π, come back shifted by half a period at 3π, and are gone again at 4π, the
                Talbot time d²/π for grating period d = 2π. The wave keeps re-imaging the grating;
                the space–time strip shows the resulting 'Talbot carpet'. It fades after a few
                cycles as the diffracted parts fly apart.
            `,
        },

        {
            label: "Coherent state", topic: "waves", preset: "harmonic", p1: 0.2,
            prep: P({ x0: -10, k: 0, sigma: 1.58 }), show: { showHist: true }, hist: 0.4,
            caption: `
                In a harmonic trap, a packet of exactly the right width sloshes back and forth
                without changing shape, in step with the classical particle.
            `,
        },

        {
            label: "Anderson localisation", topic: "waves", preset: "speckle", p1: 1, p2: 1, speed: 30,
            refs: [["Billy et al., Nature 2008", "https://arxiv.org/abs/0804.1621"]],
            prep: P({ x0: -1.6, k: 0, sigma: 3 }), show: { showHist: true }, hist: 1,
            caption: `
                In a random speckle potential the packet spreads a little and then simply stops:
                waves scattered along all the possible paths interfere and cancel the spreading. A
                small, more energetic part escapes; without the disorder, all of it would keep
                spreading.
            `,
        },

        // ---- lattices
        {
            label: "Bloch oscillations", topic: "lattice", preset: "lattice", p1: 1.5, p2: 0.05, speed: 30,
            refs: [["Ben Dahan et al., PRL 1996", "https://doi.org/10.1103/PhysRevLett.76.4508"], ["Gustavsson et al., PRL 2008", "https://arxiv.org/abs/0710.5083"]],
            prep: P({ x0: 8, k: 0, sigma: 5 }), show: { showCl: false, showP: true, showHist: true }, hist: 1,
            caption: `
                The tilt pushes the particle downhill with a constant force, yet instead of sliding
                away it swings back and forth over many lattice sites and returns every 2π/(F·d) ≈
                63 time units. The momentum strip shows why: the force changes the momentum
                steadily, but whenever it reaches the edge of the Brillouin zone (k = ±π/2), the
                lattice Bragg-reflects it to the opposite edge, so it never runs away.
            `,
        },

        {
            label: "Landau–Zener tunnelling", topic: "lattice", preset: "lattice", p1: 0.7, p2: 0.08, speed: 30,
            refs: [["Anderson & Kasevich, Science 1998", "https://doi.org/10.1126/science.282.5394.1686"]],
            prep: P({ x0: 15, k: 0, sigma: 5 }), show: { showHist: true, showP: true, showCl: false }, hist: 1,
            caption: `
                In a shallower lattice, the particle still Bloch-oscillates, but at each turnaround
                part of it tunnels into a higher energy band and slides away downhill. It leaks out
                in regular pulses (streaks in the space–time strip), the principle behind the first
                pulsed atom laser (1998).
            `,
        },

        {
            label: "Conveyor belt", topic: "lattice", preset: "belt", p1: 14, p2: 60, speed: 20,
            refs: [["Schrader et al., 2001", "https://arxiv.org/abs/quant-ph/0107029"]],
            prep: P({ x0: -15, k: 0, sigma: 4 }), show: { showHist: true }, hist: 0.7,
            caption: `
                A deep optical lattice is slid 30 units to the right and stopped. Each atom sits in
                a well and oscillates there far faster than the lattice moves, so it simply rides
                along, like on a conveyor belt; most of them arrive. Lower the depth to 2 and
                restart: the wells are now too shallow to hold anything, and the lattice slides past
                underneath and leaves the atoms behind.
            `,
        },

        // ---- condensates
        {
            label: "Faraday waves", topic: "bec", preset: "faraday", p1: 0.4, p2: 2, g: 20, speed: 30,
            refs: [["Engels et al., PRL 2007", "https://arxiv.org/abs/cond-mat/0701028"]],
            prep: G({}), show: { showHist: true, showP: true }, hist: 1.5, histMode: "change",
            caption: `
                The interaction strength pulses up and down. For a while nothing seems to happen,
                but ripples that oscillate at half the drive frequency are amplified with every
                cycle. After t ≈ 120 a regular pattern of stripes appears, spaced by ≈ 5.2: the
                wavelength of sound waves at that frequency. These are Faraday waves, first seen on
                vibrated liquids and in a condensate in 2007. Set Ω to 1.2 and restart: the stripes
                come out wider (≈ 7.5).
            `,
        },

        {
            label: "Ground state", topic: "bec", preset: "harmonic", p1: 0.2, g: 10, imaginary: true, speed: 30,
            prep: P({ x0: -8, k: 0, sigma: 1.6 }),
            caption: `
                Imaginary time relaxes any state into the ground state. With repulsive interactions,
                the condensate spreads into a broad, flat-topped Thomas–Fermi profile.
            `,
        },

        {
            label: "Time of flight", topic: "bec", preset: "free", g: 20,
            refs: [["Anderson et al., Science 1995", "https://doi.org/10.1126/science.269.5221.198"]],
            prep: G({ prepIn: "harmonic", pa: 0.1 }), show: { showHist: true }, hist: 0.3,
            caption: `
                A condensate prepared in a trap is released into free space and flies apart, pushed
                outwards by its own interactions. Set g to 0 and restart: it starts out smaller and
                expands far more slowly. Photographing condensates after such a “time of flight” is
                how most experiments see them.
            `,
        },

        {
            label: "Two condensates interfere", topic: "bec", preset: "free", g: 1,
            refs: [["Andrews et al., Science 1997", "https://doi.org/10.1126/science.275.5300.637"]],
            prep: G({ prepIn: "twotraps", pa: 1, pb: 30, zoom: 0.2 }), show: { showHist: true }, hist: 0.3,
            caption: `
                Two condensates prepared in separate traps are released. Where they overlap, they
                form interference fringes: the 1997 experiment that showed condensates are coherent
                matter waves.
            `,
        },

        {
            label: "Speed of sound", topic: "bec", preset: "box", p1: 60, g: 20,
            refs: [["Andrews et al., PRL 1997", "https://doi.org/10.1103/PhysRevLett.79.553"]],
            prep: G({ imprint: "bump", x0: 0, eps: 0.15, bw: 4 }), show: { showHist: true }, hist: 0.5, histMode: "change",
            caption: `
                A small bump on a uniform condensate splits into two pulses that travel at the speed
                of sound, c = √(g n) ≈ 0.6 here. Increase g and restart: the condensate gets stiffer
                and the sound gets faster. In the space–time strip, which shows how the density has
                changed since the start, the pulses are straight lines whose slope is the speed of
                sound.
            `,
        },

        {
            label: "Shock waves", topic: "bec", preset: "box", p1: 60, g: 20,
            refs: [["Hoefer et al., PRA 2006", "https://arxiv.org/abs/cond-mat/0603389"]],
            prep: G({ imprint: "bump", x0: 0, eps: 2, bw: 2 }), show: { showHist: true }, hist: 0.5, histMode: "change",
            caption: `
                A large bump does not split into smooth sound pulses: their fronts steepen and break
                up into trains of ripples, a dispersive shock wave. Lower the height ε towards 0.15
                and restart to get ordinary sound again.
            `,
        },

        {
            label: "Breathing mode", topic: "bec", preset: "harmonic", p1: 0.1, g: 20, speed: 20,
            refs: [["Moritz et al., PRL 2003", "https://arxiv.org/abs/cond-mat/0307607"]],
            prep: G({ prepIn: "harmonic", pa: 0.154 }), show: { showHist: true }, hist: 1,
            caption: `
                A condensate prepared in a slightly tighter trap breathes in and out once released
                into this one. With interactions it breathes at √3·ω (one breath every ≈ 36 time
                units); set g to 0 and restart, and it breathes at 2ω (every ≈ 31). Experiments use
                exactly this to measure interactions.
            `,
        },

        {
            label: "Self-trapping", topic: "bec", preset: "doublewell", p1: 0.3, p2: 8, g: 0.4, speed: 30,
            refs: [["Albiez et al., PRL 2005", "https://arxiv.org/abs/cond-mat/0411757"]],
            prep: { kind: "leftwell" }, show: { showHist: true }, hist: 2,
            caption: `
                The double well again, now with interacting bosons. Tunnelling is still possible,
                but above a threshold the interactions lock the atoms in the well they started in.
                Set g to 0.1 and restart to see them tunnel back and forth again.
            `,
        },

        // ---- solitons
        {
            label: "Dark soliton", topic: "solitons", preset: "harmonic", p1: 0.1, g: 20, speed: 20,
            refs: [["Burger et al., PRL 1999", "https://arxiv.org/abs/cond-mat/9910487"]],
            prep: G({ imprint: "soliton", x0: -3.6, v: 0 }), show: { showPh: true, showHist: true }, hist: 1,
            caption: `
                A dip in the condensate with a phase step of about π across it keeps its shape and
                oscillates through the cloud, more slowly than the cloud itself would. In the phase
                strip the step points down while the soliton moves right and up while it moves left,
                passing through π each time it turns around.
            `,
        },

        {
            label: "Grey soliton", topic: "solitons", preset: "box", p1: 30, g: 20, speed: 20,
            refs: [["Burger et al., PRL 1999", "https://arxiv.org/abs/cond-mat/9910487"]],
            prep: G({ imprint: "soliton", x0: -5, v: 0.5 }), show: { showHist: true, showPh: true }, hist: 1,
            caption: `
                In a flat condensate a black soliton stands still, but a shallower grey one moves,
                here at half the speed of sound. Its phase step is smaller than π: the shallower the
                dip, the faster it moves and the smaller the step. The step points down while it
                moves right and up while it moves left. It travels without spreading and bounces off
                the walls of the box.
            `,
        },

        {
            label: "Colliding condensates", topic: "solitons", preset: "harmonic", p1: 0.1, g: 20, speed: 20,
            refs: [["Weller et al., PRL 2008", "https://arxiv.org/abs/0803.4352"]],
            prep: P({ shape: "pair", x0: 0, k: 0, sigma: 2, d: 24, phi: 0 }), show: { showHist: true, showPh: true }, hist: 1,
            caption: `
                Two condensates fall towards each other in a trap, interfere where they overlap, and
                pass through each other. At first the fringes vanish again as the clouds separate,
                but each collision leaves a few dark solitons behind. After three or four passes the
                cloud is full of them and never separates cleanly again. Set g to 0 and restart: the
                fringes come and go with every pass, forever.
            `,
        },

        {
            label: "Bright soliton", topic: "solitons", preset: "free", g: -2,
            refs: [["Khaykovich et al., Science 2002", "https://arxiv.org/abs/cond-mat/0205378"], ["Strecker et al., Nature 2002", "https://arxiv.org/abs/cond-mat/0204532"]],
            prep: P({ shape: "sech", x0: -15, k: 1.5 }), show: { showHist: true }, hist: 0.3,
            caption: `
                With attractive interactions, the pull between bosons exactly cancels the spreading,
                so the packet travels without changing shape. Set g to 0 to watch it disperse.
            `,
        },

        {
            label: "Soliton collisions", topic: "solitons", preset: "free", g: -2,
            refs: [["Nguyen et al., Nat. Phys. 2014", "https://arxiv.org/abs/1407.5087"]],
            prep: P({ shape: "sech2", x0: 0, k: 0.3, d: 16, phi: 0 }), show: { showHist: true }, hist: 0.5,
            caption: `
                Two bright solitons heading for each other in phase merge briefly into one tall peak
                and pass straight through. Set the phase φ to π (out of phase): they repel and turn
                around without ever touching.
            `,
        },

        {
            label: "Soliton interferometer", topic: "solitons", preset: "trapbarrier", p1: 0.78, p2: 0, g: -1, speed: 30,
            refs: [["McDonald et al., PRL 2014", "https://arxiv.org/abs/1403.3485"]],
            prep: P({ shape: "sech", x0: -20, k: 0 }), show: { showHist: true }, hist: 1,
            caption: `
                A bright soliton slides down the trap into a narrow barrier that splits it in two.
                The halves swing out, come back after one trap period and meet at the barrier again,
                which recombines them. With both sides equal, 96% leaves to the right. Raise the
                right side by ε = 0.025 and restart: that half picks up an extra phase, and now 96%
                leaves to the left.
            `,
        },

        {
            label: "Soliton train", topic: "solitons", preset: "box", p1: 60, g: -6,
            refs: [["Strecker et al., Nature 2002", "https://arxiv.org/abs/cond-mat/0204532"]],
            prep: G({ sameG: false, g0: 20 }), show: { showHist: true }, hist: 0.5,
            caption: `
                A wide condensate is prepared with repulsive interactions (g₀ = 20), which are then
                suddenly made attractive (g = −6). The cloud is unstable: ripples grow and it breaks
                up into a train of bright solitons, as in experiments in 2002.
            `,
        },

        // ---- superfluidity
        {
            label: "Superfluid flow", topic: "superfluid", preset: "obstacle", p1: 0.3, p2: 0.3, g: 20, speed: 60,
            refs: [["Raman et al., PRL 1999", "https://arxiv.org/abs/cond-mat/9909109"]],
            prep: G({}), show: { showHist: true }, hist: 2.5, histMode: "change",
            caption: `
                An obstacle is dragged through the condensate at 0.3 times the speed of sound. The
                condensate simply flows around it: nothing is left behind and no energy is
                deposited, so there is no friction at all. This is superfluidity.
            `,
        },

        {
            label: "Critical velocity", topic: "superfluid", preset: "obstacle", p1: 0.8, p2: 0.3, g: 20, speed: 30,
            refs: [["Raman et al., PRL 1999", "https://arxiv.org/abs/cond-mat/9909109"]],
            prep: G({}), show: { showHist: true }, hist: 1, histMode: "change",
            caption: `
                The same obstacle at 0.8 times the speed of sound. Above a critical speed (here
                about half the speed of sound) superfluidity breaks down: the obstacle sheds dark
                solitons and sound in its wake and heats the condensate. Watch the energy line rise.
            `,
        },

        {
            label: "Persistent current", topic: "superfluid", preset: "ring", p1: 0.7, g: 20, speed: 30,
            refs: [["Ryu et al., PRL 2007", "https://arxiv.org/abs/0709.0012"], ["Ramanathan et al., PRL 2011", "https://arxiv.org/abs/1101.0019"]],
            prep: G({ imprint: "kick", wind: 2 }), show: { showPh: true, showP: true },
            caption: `
                A condensate flowing around a ring: its phase climbs steadily by 4π = 2 × 2π on the
                way around (w = 2, see the phase strip), steepest at the weak link where the flow is
                fastest. It flows past the weak link without slowing down, indefinitely. Its
                circulation can only take whole-number values, and can only change where the density
                drops to zero. Set w to 5 and restart: that flow is too fast for the weak link,
                which sheds solitons, the circulation jumps by whole steps (phase slips), and the
                current dies away.
            `,
        },

        {
            label: "Phase slips", topic: "superfluid", preset: "ringstir", p1: 0.6, p2: 0.8, g: 20, speed: 40,
            refs: [["Wright et al., PRL 2013", "https://arxiv.org/abs/1208.3608"]],
            prep: G({}), show: { showHist: true, showP: true }, hist: 3, histMode: "change",
            caption: `
                A weak link is dragged around the ring at 0.6 times the speed of sound. The
                condensate cannot speed up smoothly: whenever the flow through the link gets too
                fast, a phase slip lets the circulation jump by exactly one, 0, 1, 2, 3, 4. At t =
                300 the link is removed, and the current keeps going on its own. Set the speed to
                0.3 and restart: a superfluid ignores slow stirring and does not rotate at all.
            `,
        },

        {
            label: "Stir it yourself", topic: "superfluid", preset: "stir", p1: 0.4, g: 20, speed: 30,
            refs: [["Raman et al., PRL 1999", "https://arxiv.org/abs/cond-mat/9909109"]],
            prep: G({}), show: { showHist: true }, hist: 0.5, histMode: "change",
            caption: `
                Drag on the plot to move the obstacle through the condensate; its speed is shown at
                the top. Move it slowly and the fluid flows around it without resistance. Move it
                faster than about half the speed of sound and it starts to shed dark solitons and
                sound.
            `,
        },

        // ---- driven
        {
            label: "Kapitza–Dirac diffraction", topic: "driven", preset: "pulse", p1: 3, speed: 5,
            refs: [["Gould et al., PRL 1986", "https://doi.org/10.1103/PhysRevLett.56.827"]],
            prep: P({ x0: 0, k: 0, sigma: 2.5 }), show: { showHist: true, showP: true }, hist: 0.15,
            caption: `
                A standing wave of light flashed on for an instant ripples the phase of the wave.
                That diffracts it into packets moving at 0, ±1, ±2, … (multiples of the lattice's
                recoil), visible as separate peaks in the momentum distribution. Their weights are
                squared Bessel functions: 27% at rest, 31% each way at ±1, 5% at ±2. Lower the pulse
                area to 1.5 and restart: 75% stays at rest.
            `,
        },

        {
            label: "Delta-kick cooling", topic: "driven", preset: "lens", p1: 1,
            refs: [["Ammann & Christensen, PRL 1997", "https://doi.org/10.1103/PhysRevLett.78.2088"], ["Kovachy et al., PRL 2015", "https://arxiv.org/abs/1407.6995"]],
            prep: G({ prepIn: "harmonic", pa: 0.5 }), show: { showHist: true, showP: true }, hist: 0.2,
            caption: `
                A cloud released from a trap expands, each part moving at a speed proportional to
                its distance from the centre. At t = 12 a harmonic potential flashed on briefly
                slows every part by just that much, like a lens collimating light. The expansion
                stops and the momentum distribution narrows sixfold: the cloud is about 34 times
                colder. Experiments use this to reach picokelvin. Set the lens power to 0 (no lens)
                or 2 (focuses, then expands again) and restart.
            `,
        },

        {
            label: "Shortcut to adiabaticity", topic: "driven", preset: "ramp", p1: 10, p2: 1,
            refs: [["Schaff et al., 2010", "https://arxiv.org/abs/1006.1495"]],
            prep: G({}), show: { showHist: true }, hist: 1,
            caption: `
                The trap is opened tenfold in just T = 10. With a plain ramp the cloud would be left
                breathing, between widths 1.6 and 6.2 (set s to 0 and restart); even a ramp 20 times
                slower leaves it breathing. The designed ramp, which for a moment even turns the
                trap into a gentle hill, delivers it at rest, exactly at the new ground-state width
                of 3.16. The space–time strip shows the difference at a glance.
            `,
        },

        {
            label: "Frozen tunnelling", topic: "driven", preset: "shaken", p1: 0.029, speed: 30,
            refs: [["Kierig et al., PRL 2008", "https://arxiv.org/abs/0803.1406"]],
            prep: { kind: "leftwell" }, show: { showHist: true }, hist: 2,
            caption: `
                The double well again, now tilted back and forth. At this drive strength (A·d/Ω =
                2.4, d the distance between the wells, a zero of the Bessel function J₀) the
                tunnelling almost stops: at least three quarters of the particle stays in the well
                it started in. This is coherent destruction of tunnelling. Set A to 0.015 and
                restart: it tunnels again, more slowly than without the drive.
            `,
        },

        {
            label: "Adiabatic passage", topic: "driven", preset: "sweep", p1: 600, speed: 40,
            prep: G({}), show: { showHist: true }, hist: 5,
            caption: `
                The particle starts in the lower (left) well, then the tilt is slowly reversed.
                Changed slowly enough, the particle follows the lowest state and ends up in the
                other well (about 96%). Set the sweep time to 40 and restart: now the change is too
                fast to follow and about 90% stays behind (Landau–Zener).
            `,
        },

        {
            label: "Paul trap", topic: "driven", preset: "paul", p1: 0.4,
            refs: [["Paul, Rev. Mod. Phys. 1990", "https://doi.org/10.1103/RevModPhys.62.531"]],
            prep: P({ x0: 3, k: 0, sigma: 1.2 }), show: { showHist: true }, hist: 0.3,
            caption: `
                A particle on top of a hill should roll off. Shaking the hill's curvature up and
                down fast enough keeps it balanced on top: Kapitza's upside-down pendulum, and how
                Paul traps hold single ions. The fast wobble is the shaking itself; the slow swing
                is the particle oscillating in the effective trap the shaking creates. Set A to 0
                and it falls off.
            `,
        },

        {
            label: "Atom interferometer", topic: "driven", preset: "bragg", p1: 0,
            refs: [["Kasevich & Chu, PRL 1991", "https://doi.org/10.1103/PhysRevLett.67.181"], ["Torii et al., PRA 2000", "https://arxiv.org/abs/cond-mat/9908160"]],
            prep: P({ x0: -25, k: 0, sigma: 6 }), show: { showP: true },
            caption: `
                A short pulse of a moving light lattice can kick a resting atom to speed 2, or stop
                a moving one. A half-length pulse does this only halfway: the atom ends up both
                resting and moving. That is the first pulse, a beam splitter, and the two halves
                drift apart. The second, full pulse swaps them (a mirror): the moving half stops,
                the resting half starts moving, so they meet again. The third, half pulse recombines
                them, and how the two halves interfere decides where the atom goes. At φ = 0, about
                95% ends at rest; at φ = π (3.14), about 92% moves on. Anything that shifts the
                phase of one path relative to the other, such as gravity, shows up in this split,
                which is how atom interferometers measure it.
            `,
        },
    ];

    return { topics, demos };
})();
