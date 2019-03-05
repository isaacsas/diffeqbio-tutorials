using DifferentialEquations, DiffEqBiological, Plots, Latexify;
pyplot();

repressilator = @reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end α K n δ γ β μ;

latexify(repressilator; env=:chemical)

latexify(repressilator)

speciesmap(repressilator)

paramsmap(repressilator)

# ODES

# parameters [α,K,n,δ,γ,β,μ]
p = (.5, 40, 2, log(2)/120, 5e-3, 20*log(2)/120, log(2)/60)

# initial condition [m₁,m₂,m₃,P₁,P₂,P₃]
u₀ = [0.,0.,0.,20.,0.,0.]

# time to solve over
tspan = (0., 10000.)

# create the ODEProblem we want to solve
oprob = ODEProblem(repressilator, u₀, tspan, p)

sol = solve(oprob, saveat=10.)
plot(sol)

# JUMPS

# first we redefine the initial condition to be integer valued
u₀ = [0,0,0,20,0,0]

# first we create a discrete problem to encode that our species are integer valued:
dprob = DiscreteProblem(repressilator, u₀, tspan, p)

# now we create a JumpProblem, and specify Gillespie's Direct Method as the solver:
jprob = JumpProblem(dprob, Direct(), repressilator, save_positions=(false,false))

# now let's solve and plot the jump process:
sol = solve(jprob, SSAStepper(), saveat=10.)
plot(sol)

# τ-leaping
lprob = JumpProblem(dprob, Direct(), repressilator.regular_jumps)
lsol = solve(lprob, SimpleTauLeaping(), dt=.1)
plot(lsol, plotdensity=1000)

# CLE SDE:
bdp = @reaction_network begin
    c₁, X --> 2X
    c₂, X --> 0
    c₃, 0 --> X
  end c₁ c₂ c₃
p = (1.0,2.0,50.)
u₀ = [5.]
tspan = (0.,4.)

# SDEProblem for CLE
sprob = SDEProblem(bdp, u₀, tspan, p)

# solve and plot
sol = solve(sprob, saveat=.004)
plot(sol)

# other stuff
latexify(repressilator.symjac)

# solve bcr example
using ReactionNetworkImporters, Sundials, DataFrames, CSVFiles
networkname = "testbcrbng"
tf = 10000.
datadir = joinpath(dirname(pathof(ReactionNetworkImporters)),"../data/bcr")
fname    = joinpath(datadir, "bcr.net")
gdatfile = joinpath(datadir, "bcr.gdat")
print("getting gdat file...")
gdatdf = DataFrame(load(File(format"CSV", gdatfile), header_exists=true, spacedelim=true) )
prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname); 
rnbng = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
addodes!(rnbng; build_jac=false, build_symfuncs=false)
boprob = ODEProblem(rnbng, u0, (0.,tf), p)
asykgroups = prnbng.groupstoids[:Activated_Syk]
asyksyms = findall(x -> x ∈ asykgroups, rnbng.syms_to_ints)
bsol = solve(boprob, CVODE_BDF(),dense=false, saveat=1., abstol=1e-8, reltol=1e-8)
basyk = sum(bsol[asykgroups,:], dims=1)
plot(gdatdf[:time][2:end], gdatdf[:Activated_Syk][2:end], xscale=:log10, label=:BNG)
plot!(bsol.t[2:end], basyk'[2:end], label=:DEBio, xscale=:log10, linestyle=:dash)
xlabel!("t")


# bifurcation plot:
using DiffEqBiological, Plots
rn = @reaction_network begin
         k1, 0 --> Y
         k2p, Y --> 0
         k2pp*P, Y --> 0
         mm(1-P,(k3p+k3pp*A),J3), 0 --> P
         mmr(P,k4*m/J4,J4), Y + P --> Y 
end k1 k2p k2pp k3p k3pp A J3 k4 m J4
p = [0.04,0.04,1.,1.,10.0,0.,0.04,35.,.1,.04]

pmap = Dict{Symbol,Float64}()
foreach(sym -> (sym != :m) && (pmap[sym] = p[rn.params_to_ints[sym]]), rn.params)
fix_parameters(rn;pmap...)

bif = bifurcations(rn,p,:m,(.01,.65))
bif_plot(bif)