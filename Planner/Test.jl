cd("C:\\Users\\dwills\\Dropbox\\1_5_TaxesAndInformality\\Codes\\OptimalTaxation\\Planner")
#cd("/Users/danielwillsr/Documents/Planner")

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


#Define initial state vector

uw0    = 0.01 #0.000000926495
μ0     = 0.0
e0     =  pa.θ_e_lb;
ϕ_e0   = -0.00009
y_agg0 =  0.0
λ0     =  1.0
l_agg0 =  0.0
ω0     =  10.0


#Test  new_find_controls
##new_find_controls( θ0, ss0, pa)

#Test states
##du0 = Array{Float64,1}(undef,8)
#u0 = Array{Float64,1}(undef,8)
#u0[1] = uw0;
#u0[2] = μ0;
#u0[3] = e0;
##u0[4] = ϕ_e0;
#u0[5] = y_agg0;
#u0[6] = λ0;
#u0[7] = l_agg0;
#u0[8] = ω0;
##find_states!(du0,u0,pa,θ0)

#Test my runge kutta
#Nspan = 500
#y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0]
#xlb= log(pa.θ_w_lb);
#xub = log(pa.θ_w_ub);
#xstep = (xub - xlb)/(Nspan - 1);
#xspan = xlb:xstep:xub;
#solution = Array{Float64}(undef,Nspan,8);
#my_runge_kutta!(solution,y0,xspan,xstep,pa)

# Solve \mu(theta_ub)=0 and plot the solution
Nspan = 100;
y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0];
uw0_low=10.0^-10;
uw0_up =1.0;
sol = find_uw0(y0, uw0_low, uw0_up, Nspan, pa);

xlb= log(pa.θ_w_lb);
xub = log(pa.θ_w_ub);
xstep = (xub - xlb)/(Nspan - 1);
xspan = xlb:xstep:xub;
θspan = exp.(xspan);

controls = Array{Float64}(undef,Nspan,4);
recover_controls!(controls, θspan, sol);
