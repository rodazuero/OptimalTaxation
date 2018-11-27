cd("C:\\Users\\d-wills\\Dropbox (Uniandes)\\1_5_TaxesAndInformality\\PlannerPblmCodes")

using Roots
using Plots
using DifferentialEquations


include("Definitions.jl")
#initialize init_parameters
pa = init_parameters();

#Partial equilibrium; define prices
ω=0.2;
λ=0.05;

uw0 = 0.01
param = (ω, pa)

# Mirrlees
include("Mirrlees.jl")
#1. Create solution object
u0 = [uw0,0.0,0.0,λ]
tspan = (exp(pa.μ_w-3*pa.σ2_w),exp(pa.μ_w+3*pa.σ2_w))
prob = ODEProblem(mirrlees!,u0,tspan,param)
solvec=Array{OrdinaryDiffEq.ODECompositeSolution}(undef,1)
solvec[1] = solve(prob)

#2. Find/ check bound where BVP is positive and negative
solvediffeq!(solvec, [80.5,0.0,0.0,λ], tspan, ω, pa)
solvediffeq!(solvec, [81.5,0.0,0.0,λ], tspan, ω, pa)

#3. Solve BVP
solvepartialeq!(solvec, [80.5, 81.5], u0, tspan, λ, ω, pa)

#4. Plot solution
θvec= tspan[1]:((tspan[2]-tspan[1])/99):tspan[2];
martaxes= Array{Float64}(undef,100);
taxes= Array{Float64}(undef,100);
bases= Array{Float64}(undef,100);
labor= Array{Float64}(undef,100);
mass=  Array{Float64}(undef,100);
i=1;
hw(θ) =pdf( Normal(pa.μ_w, pa.σ2_w), log(θ));
for θ in θvec
    martaxes[i], taxes[i], bases[i], labor[i] = recover(θ, solvec[1], λ, ω, pa);
    mass[i]= hw(θ);
    global i=i+1;
end

plot(bases,martaxes)
plot(bases, taxes)
