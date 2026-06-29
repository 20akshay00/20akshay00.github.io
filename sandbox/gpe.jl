using Plots, FFTW, Statistics, LinearAlgebra

function init_cond(xgrid, μ, σ, p)
    return @. Complex(1/(σ*sqrt(2*π))^0.5 * exp(p*im*xgrid - ((xgrid - μ)^2/(4*σ^2))))
end

function solve_gpe_imaginary_time(dx, x_max, dt, n_steps, mass, g, potential_fn)
    xgrid = collect((-x_max):dx:(x_max-dx))
    n = length(xgrid)
    V_vals = potential_fn.(xgrid)

    dk = 2π / (2 * x_max)
    kgrid = fftshift((-n÷2):(n÷2-1)) * dk

    ψs = zeros(ComplexF64, n, n_steps)
    μs = Float64[]

    psi = init_cond(xgrid, 3.0, 0.5, 0.0)
    psi /= sqrt(sum(abs2.(psi)) * dx)

    K_prop = exp.(-0.5 * dt * (1/mass) .* kgrid .^ 2)

    for j in 1:n_steps
        ψs[:, j] .= psi

        psi .*= exp.(-0.5 * dt .* (V_vals .+ g .* abs2.(psi)))

        psi_k = fft(psi)
        psi_k .*= K_prop
        psi = ifft(psi_k)

        psi .*= exp.(-0.5 * dt .* (V_vals .+ g .* abs2.(psi)))

        norm_factor = sqrt(sum(abs2.(psi)) * dx)
        psi /= norm_factor

        if j > 1
            push!(μs, -log(norm_factor) / dt)
        end
    end

    return xgrid, V_vals, ψs, μs
end

begin
    V_harmonic(x) = 0.5 * x^2
    dx, dt = 0.05, 0.002
    x_limit = 10.0
    n_steps = 2500
    mass = 1.0
    g = 20.0
end

xgrid, V_vals, ψs, μs = solve_gpe_imaginary_time(dx, x_limit, dt, n_steps, mass, g, V_harmonic)

anim = @animate for j = 1:5:n_steps
    p1 = plot(xgrid, abs2.(ψs[:, j]),
        fill=(0, 0.2, :purple),
        title="",
        lab = "",
        ylabel="Probability density (|ψ(x)|²)",
        xlabel = "Position (x)",
        ylim=(-0.05, 1.0),
        lw=2,
        color=:purple,
        framestyle=:box)

    plot!(p1, xgrid, V_vals, color=:black, alpha=0.3, label="Potential V(x)")

    current_idx = max(1, j-1)
    p2 = plot(1:current_idx, μs[1:current_idx],
        title="",
        xlabel="Iteration",
        ylabel="Chemical potential (μ)",
        label="",
        color=:red,
        lw=1.5,
        xlim=(0, n_steps),
        ylims=(4.8, 16.),
        framestyle=:box)

    plot(p1, p2, size=(900, 400), layout=(1, 2), margins=5Plots.mm)
end

mp4(anim, "../images/research/gpe_imag_time.mp4", fps=30)