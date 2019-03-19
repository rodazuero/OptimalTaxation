cd("C:\\Users\\d-wills\\Dropbox (Uniandes)\\1_5_TaxesAndInformality\\PlannerPblmCodes")

using Roots
#using Plots
using DifferentialEquations


include("Definitions.jl")
include("Controls.jl")
include("States.jl")
include("MyFindZeros.jl")

#Parameters
pa = init_parameters();

#u0

uw0    = 0.025
μ0     = 0.0
e0     =  exp(pa.μ_e-3.0*pa.σ2_e)
ϕ_e0   = -0.002697
y_agg0 =  0.0
λ0     =  1.08
l_agg0 =  0.0
ω0     =  0.015

θ0= exp(pa.μ_w-3.0*pa.σ2_w)

n_FI = (λ0 * pa.α * e0 / ω0)^(1.0/(1-pa.α))
A0 = - e0*n_FI^pa.α + ω0/λ0*n_FI;
ϕ_e0 = (λ0*A0 - uw0^pa.ϕ + λ0*uw0)* pa.he(θ0, e0)

l_FI= (ω0*θ0/ (λ0*pa.χ))^(1.0/pa.ψ)

ϕ_e0 = (-0.1*λ0- uw0^pa.ϕ + λ0*uw0)* pa.he(θ0, e0)

init_controls( uw0, e0, ϕ_e0, λ0, ω0, pa);

 #1. Create solution object
 u0 = [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0]
 tspan = (θ0,exp(pa.μ_w+3.0*pa.σ2_w))
 alg= RK4();
 prob = ODEProblem(find_states!,u0,tspan,pa)
 solvec = solve(prob,alg, dt=0.01, adaptive= false)


pa.σ= 0.2
pa.α= 0.7
pa.β= 0.2
e0  = 2.5
λ0  = 1.0
ω0  = 0.15

plot  [(0.5*x^0.7 + 0.16/0.7*0.15*x  + 0.1*1.2)^1.2 ;  1.2*0.1/ 0.2 + 1.2/0.2 * 0.3/0.7 *0.15*x ] ; x from 0 to 40

θ=  1.332913470700001
e =0.15312099274670196; uw= 28.410744445290298; ϕ_e= -0.36956890658247793; μ = 5.267834902222252e-6;
    ss = State(e, uw, ϕ_e, μ, λ0, ω0);
