using Plots, LaTeXStrings
pyplot()

using Random
nsteps = 20
Y = zeros(Int, nsteps + 1)
T = zeros(Float64, nsteps + 1)
for i = 1:nsteps
    T[i+1] = T[i] + randexp()
    Y[i+1] = Y[i] + 1
end

# plot(T,Y,linetype=:steppost, xlabel=L"t", ylabel=L"Y(t)", fmt=:svg, legend=false)
# plot(T,Y,linetype=:steppost, xlabel=L"t", ylabel=L"Y(t)", fmt=:pdf, legend=false)
# fname = string(pwd(),"/figs/unitpoisson.pdf")
# savefig(fname)

using DiffEqBiological, Latexify
rn = @reaction_network begin
   (k₊,k₋), A + B ↔ C
    α, 0 --> C
    β, A --> 0
    γ, B --> 0
end k₊ k₋ α β γ
latexify(rn; env=:chemical)

latexify(rn)

using OrdinaryDiffEq
u0 = [0., 50., 100.]
p  = [1., 1., 1., 1., 1.]
tspan = (0.,30.)
oprob = ODEProblem(rn, u0, tspan, p)
osol = solve(oprob, Tsit5())
p1 = plot(osol, fmt=:svg, linestyle=:dash, xlabel=L"t")
xlabel!(p1,L"t")

using DiffEqJump
u0 = [0, 50, 100]
dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(dprob, Direct(), rn)
jsol = solve(jprob, SSAStepper())
# p2 = plot(osol, fmt=:svg, linestyle=:dash, xlabel=L"t")
# plot!(p2, jsol, fmt=:svg, xlabel=L"t")
# plot(p1,p2, size=(800,400))

# fmt = :pdf
# p1 = plot(osol, fmt=fmt, linestyle=:dash, xlabel=L"t")
# p2 = plot(osol, fmt=fmt, linestyle=:dash, xlabel=L"t")
# plot!(p2, jsol, fmt=fmt, xlabel=L"t")
# plot(p1,p2, size=(800,400))
# fname = string(pwd(),"/abtoc.pdf")
# savefig(fname)