#cd("C:\\Users\\dwills\\Dropbox\\1_5_TaxesAndInformality\\Codes\\OptimalTaxation\\Planner")
#cd("/Users/danielwillsr/Documents/Planner")
cd("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Planner")

using Roots
using NLopt
using Statistics
using PyPlot
using DataFrames
using CSV
using NLsolve
#using Plots
#using DifferentialEquations

#Old parameters
include("Definitions.jl") #Definitions2.jl
include("Controls.jl") #NewControls.jl
include("States.jl") #States_Ralwsian.jl
#include("MyFindZeros.jl")
include("MyRungeKutta.jl")
#include("NewMyRungeKutta.jl") #NewMyRungeKutta.jl
include("BoundaryValueProblem.jl")

#Define values for the model parameters
pa = init_parameters();

#Define initial state vector

uw0    = 1 #guess
μ0     =-1
e0     =  pa.θ_e_lb;
ϕ_e0   =  -3.2 #guess
#ϕ_e0   =  -1.233e-5 #Min- e=1.07
y_agg0 =  0.0
λ0     =  1 #guess
l_agg0 =  0.0
ω0     =  1.6 #guess=#
l_new0 =  0.0
y_new0 =  0.0
ystar0 = 0.0
lstar0 = 0.0 #n=#


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

Nspan = 500
y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0, 0, 0 ];
xlb= pa.θ_w_lb;
xub = pa.θ_w_ub;
xstep = (xub - xlb)/(Nspan - 1);
xspan = xlb:xstep:xub;
solution = Array{Float64}(undef,Nspan,12);
agg = Array{Float64}(undef,Nspan,4);
my_runge_kutta!(solution,y0,xspan,xstep,pa)
#my_runge_kutta!(solution,[uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0],xspan,xstep,pa)

solution[end,:]
using DelimitedFiles
writedlm("SolutionNew.csv",solution,';')

#=xlb= log(pa.θ_w_lb);
xub = log(pa.θ_w_ub);
xstep = (xub - xlb)/(Nspan - 1);
xspan = xlb:xstep:xub;
θspan = exp.(xspan);=#

θspan = Array{Float64,1}
θspan = collect(xspan)


controls = Array{Float64}(undef,Nspan,4);
recover_controls!(controls, θspan, solution);
controls

#FINAL BOUNDARY CONDITIONS

#Boundary for μ
    #Utilitarian planner
    #bound_μ= (pa.φ*solution[end,1]^(pa.φ-1)-solution[end,6])*integral
    #Ralwsian planner
    bound_μ= (pa.φ*solution[end,1]^(pa.φ-1)-solution[end,6])*integral
    dd_μ= solution[end,2]-bound_μ=
#Boundary for e
    bound_e= -((solution[end,1]^pa.ϕ+solution[end,6]*((solution[end,3]*controls[end,2]^pa.α)-((pa.β*(controls[end,1])^(1+pa.σ))/(1+pa.σ))-solution[end,1])-solution[end,8]*controls[end,2])*pa.he(xub, solution[end,3]))
    dd=bound_e-solution[end,4]


C= DataFrame(controls)
names!(C,[:z, :n, :l, :p] )
RungeKutta=hcat(DataFrame(solution),C, DataFrame(theta=θspan))
CSV.write("RungeKuttaNew.csv", RungeKutta)

# Test distance

#fixed_initial_states= [μ0, e0, y_agg0, l_agg0]ju
fixed_initial_states= [μ0, e0, y_agg0, l_agg0, y_new0, l_new0, ystar0, lstar0 ]
x = [uw0, ϕ_e0, λ0, ω0]
grad = [];

io=open("Results.txt", "w")
println(io, "uw0, phie_0, lambda0, omega0, mu, delta_e, Y, L, distance")
close(io);

dist= distance!(x, grad, fixed_initial_states, Nspan, pa, solution, io);
objective(x, grad) = distance!(x, grad, fixed_initial_states, Nspan, pa, solution, io);

# Shooting (documentacion en NLopt.jl)
n = length(x);

#Local optimization
local_opt = Opt(:LN_NELDERMEAD, n)
local_opt.xtol_rel = 1e-3
local_opt.min_objective = objective

#=(minf,minx,ret) = optimize(local_opt, [uw0, ϕ_e0, λ0, ω0])
numevals = local_opt.numevals # the qnumber of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")
println("got $minf at $minx after $numevals iterations (returned $ret)")=#


#Feasible region for global optimization

lb1 = [eps(), -0.1, 14, 8]
ub1 = [ 0.0015, 0.0, 16, 12]
lb2 = [eps(), -0.1, 14, 8.3]
ub2 = [ 0.0014, 0, 14.5, 10.5]
lb3 = [eps(), -0.1, 14, 8]
ub3 = [ 0.0014, 0, 31, 15.5]
lb4 = [eps(), -0.1, 12, 8]
ub4 = [ 0.001, 0, 16, 10]

lb6 = [eps(), -0.1, 13, eps()]
ub6 = [ 0.001, 0, 16, 1.0]

lbr2 = [eps(), -0.3, 10, 6]
ubr2 = [ 0.001, 0, 18, 12]


#Global optimization
opt = Opt(:G_MLSL_LDS, n)
opt.lower_bounds = lbr2
opt.upper_bounds = ubr2
opt.local_optimizer = local_opt::Opt
opt.xtol_rel = 1e-3
opt.min_objective = objective
opt.maxtime=3600

(minf,minx,ret) = optimize(opt, [uw0, ϕ_e0, λ0, ω0])
numevals = local_opt.numevals # the qnumber of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")
println("got $minf at $minx after $numevals iterations (returned $ret)")


#Recover the initial states values
uw0 = minx[1];
ϕ_e0= minx[2];
λ0 = minx[3];
ω0 = minx[4];

# Solve \mu(theta_ub)=0 and plot the solution
Nspan = 500
y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0, 0, 0 ];
xlb= log(pa.θ_w_lb);
xlb= pa.θ_w_lb;
xub = log(pa.θ_w_ub);
xub = pa.θ_w_ub;
xstep = (xub - xlb)/(Nspan - 1);
xspan = xlb:xstep:xub;
solution = Array{Float64}(undef,Nspan,12);
agg = Array{Float64}(undef,Nspan,4);
slope=((pa.θ_e_ub-pa.θ_e_lb)*(Nspan-1))/((log(pa.θ_w_ub)-log(pa.θ_w_lb))*Nspan)
my_runge_kutta!(solution,y0,xspan,xstep,pa)

θspan = Array{Float64,1}
θspan = collect(xspan)

controls = Array{Float64}(undef,Nspan,4);
recover_controls!(controls, θspan, solution);
controls
bound_e= - (solution[end,1]^pa.ϕ+solution[end,6]*((solution[end,3]*controls[end,2]^pa.α)-((pa.β*(controls[end,1])^(1+pa.σ))/(1+pa.σ))-solution[end,1])-solution[end,8]*controls[end,2])*pa.he(xub, solution[end,3])


##GRAPHS
pygui(true)
rc("font", family="serif")

#States plot
#=for i in 1:4
    x=i
    estados=subplot(42i)
    estados[1,1]=plot(θspan, solution[:,i])
end=#

#States
estados=subplot(421)
suptitle("Optimal states")
estados=plot(θspan[1:500], solution[1:500,1])
ylabel("u_w")
estados=subplot(422)
estados[1,2]=plot(θspan, solution[:,2])
plot(θspan,repeat([0],500),color="lime")
ylabel("μ")
estados=subplot(423)
estados[2,1]=plot(θspan, solution[:,3])
plot(θspan,repeat([pa.θ_e_ub],500),color="lime")
ylabel("e")
estados=subplot(424)
estados[2,2]=plot(θspan, solution[:,4])
#plot(θspan,repeat([bound_e],500),color="lime")
plot(θspan,repeat([0],500),color="lime")
ylabel("φe")
estados=subplot(425)
estados[3,1]=plot(θspan, solution[:,5])
ylabel("Y")
estados=subplot(426)
estados[3,2]=plot(θspan, solution[:,6])
ylabel("λ")
estados=subplot(427)
estados[4,1]=plot(θspan, solution[:,7])
plot(θspan,repeat([0],500),color="lime")
ylabel("L")
estados=subplot(428)
estados[4,2]=plot(θspan, solution[:,8])
ylabel("ω")

#Markets size
estados=subplot(121)
estados[1,2]=plot(θspan, solution[:,10])
ylabel("Y*")
estados=subplot(122)
estados[1,2]=plot(θspan, solution[:,9])
plot(θspan,repeat([0],500),color="lime")
ylabel("L*")

#Auxiliar states
estados=subplot(121)
estados[1,2]=plot(θspan, solution[:,11])
plot(θspan,repeat([0],500),color="lime")
ylabel("L*")
estados=subplot(122)
estados[1,2]=plot(θspan, solution[:,12])
plot(θspan,repeat([0.15],500),color="lime")
ylabel("Y*")

#Controls
ctrl=subplot(221)
suptitle("Optimal controls")
ctrl=plot(θspan, controls[:,1])
ylabel("z")
ctrl=subplot(222)
ctrl[1,2]=plot(θspan, controls[:,2])
ylabel("n")
ctrl=subplot(223)
ctrl[2,1]=plot(θspan, controls[:,3])
ylabel("l")
ctrl=subplot(224)
ctrl[2,2]=plot(θspan, controls[:,4])
ylabel("p")
