cd("C:\\Users\\d-wills\\Dropbox (Uniandes)\\1_5_TaxesAndInformality\\Codes\\OptimalTaxation\\Planner")

using Roots
#using Plots
#using DifferentialEquations


include("Definitions.jl")
include("Controls.jl")
include("States.jl")
include("MyFindZeros.jl")

#Parameters
pa = init_parameters();


#u0

uw0    = 0.00001
μ0     = 0.0
e0     =  pa.θ_e_lb;
ϕ_e0   = -0.1
y_agg0 =  0.0
λ0     =  1.00
l_agg0 =  0.0
ω0     =  0.055

θ0 = pa.θ_w_lb;
ss0 = State(e0, uw0, ϕ_e0, μ0, λ0, ω0);

h_e0= pa.he(θ0, e0)
h_w0= pa.hw(θ0, e0)

A0= ((uw0^pa.ϕ-λ0*uw0)*h_e0 + ϕ_e0) / (λ0*h_e0)
#Test  new_find_controls

new_find_controls( θ0, ss0, pa)



 #1. Create solution object
 u0 = [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0]
 tspan = (θ0,5.0)
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
