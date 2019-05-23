cd("C:\\Users\\d-wills\\Dropbox (Uniandes)\\1_5_TaxesAndInformality\\Codes\\OptimalTaxation\\Planner")

using Roots
#using Plots
#using DifferentialEquations


include("Definitions.jl")
include("Controls.jl")
include("States.jl")
#include("MyFindZeros.jl")
include("MyRungeKutta.jl")
include("BoundaryValueProblem.jl")

#Parameters
pa = init_parameters();


#u0

uw0    = 0.000000926495
μ0     = 0.0
e0     =  pa.θ_e_lb;
ϕ_e0   = -0.04
y_agg0 =  0.0
λ0     =  1.0
l_agg0 =  0.0
ω0     =  5.0

θ0 = pa.θ_w_lb;
ss0 = State(e0, uw0, ϕ_e0, μ0, λ0, ω0);

h_e0= pa.he(θ0, e0)
h_w0= pa.hw(θ0, e0)

A0= ((uw0^pa.ϕ-λ0*uw0)*h_e0 + ϕ_e0) / (λ0*h_e0)
#Test  new_find_controls
##new_find_controls( θ0, ss0, pa)

#Test states
##du0 = Array{Float64,1}(undef,8)
u0 = Array{Float64,1}(undef,8)
u0[1] = uw0;
u0[2] = μ0;
u0[3] = e0;
u0[4] = ϕ_e0;
u0[5] = y_agg0;
u0[6] = λ0;
u0[7] = l_agg0;
u0[8] = ω0;
##find_states!(du0,u0,pa,θ0)

#Test my runge kutta
Nspan = 500
xlb= log(pa.θ_w_lb);
xub = log(pa.θ_w_ub);
xstep = (xub - xlb)/(Nspan - 1);
xspan = xlb:xstep:xub;

solution = Array{Float64}(undef,Nspan,8);
my_runge_kutta!(solution,u0,xspan,xstep,pa)
solution[end, 3]


 final_mu([uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0] ,500, pa)
 final_mu([0.01, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0] ,500, pa)

find_uw0(u0, uw0, 0.01, Nspan, pa)

 #1. Create solution object
 #u0 = [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0]
 #tspan = (pa.θ_w_lb,10.0)
 #alg= RK4();
 #prob = ODEProblem(find_states!,u0,tspan,pa)
 #solvec = solve(prob,alg, dt=0.01, adaptive= false)
