using Plots, LinearAlgebra

"""
    init_cond(xgrid, μ, σ, p)
Generates a moving Gaussian wavepacket.
"""
function init_cond(xgrid, μ, σ, p)
    return @. Complex(1/(σ*sqrt(2*π))^0.5 * exp(p*im*xgrid - ((xgrid - μ)^2/(4*σ^2))))
end

"""
    solve_schrodinger(dx, xi, xf, m, dt, mass, packet_params, potential_fn)
Solves the TDSE using the Crank-Nicolson method.
"""
function solve_schrodinger(dx, xi, xf, m, dt, mass, packet_params, potential_fn)
    xgrid = collect(xi:dx:xf)
    n = length(xgrid)
    V_vec = potential_fn.(xgrid)

    # constructing the evolution matrices 
    α = im*dt/(4*mass*dx^2) * ones(Complex, n-1)
    ξ = 1 .+ (im * dt/2) * (1/(mass*dx^2) .+ V_vec)
    γ = 1 .- (im * dt/2) * (1/(mass*dx^2) .+ V_vec)

    M₁ = SymTridiagonal(ξ, -α)
    M₂ = SymTridiagonal(γ, α)

    ψ = zeros(Complex, (n, m))
    ψ[:, 1] = init_cond(xgrid, packet_params...)

    for j in 1:(m-1)
        ψ[:, j+1] = M₁ \ (M₂ * ψ[:, j])
    end

    return xgrid, V_vec, ψ
end

# simulation parameters
begin
    V_barrier(x) = (x > -5.0 && x < 5.0) ? 2.5 : 0.0

    nSteps = 2200
    dx, dt = 0.01, 0.025
    x_start, x_end = -40.0, 40.0
    mass = 5.0
    packet = (-25.0, 2.0, 6.0) # μ, σ, p
end

# run simulation
xgrid, V_vals, ψs = solve_schrodinger(dx, x_start, x_end, nSteps, dt, mass, packet, V_barrier)

# generate animation
anim = @animate for j = 1:10:nSteps
    plot(xgrid, abs2.(ψs[:, j]),
        ylim=(-0.05, 0.35),
        label="",
        xlabel="Position (x)",
        ylabel="Probability density (|ψ(x)|²)",
        lw=1.5,
        fill=(0, 0.15, :blue),
        framestyle=:box,
        title="")

    plot!(xgrid, V_vals,
        color=:black,
        lw=1,
        label="Potential V(x)",
        alpha=0.5)
end

mp4(anim, "./images/research/tdse_wavepacket.mp4", fps=30)