module MySimpleModel

using Roots
using NLopt
using Statistics
#using PyPlot
using DataFrames
using CSV
using NLsolve
using DifferentialEquations
#using Plots

#Global Problem
include("Definitions2.jl")
include("MJControlsAlgorithmGlobal3.jl")
include("States_NewRungeKutta.jl")
include("NewMyRungeKuttaReverse.jl")
#Entrepreneurs Problem
include("NewControlsAlgorithmE.jl")
include("States_NewRungeKuttaE.jl")
include("NewMyRungeKuttaEReverse.jl")

#Define values for the model parameters
pa = init_parameters();

#Entrepreneurs Problem
    #Define proportion of agents in global problem
    gp     =   0.993
    ue0    =   640.0
    μe0    =   0.0 - 1.0e-10
    ye0    =   0.0
    λe0    =   1.0
    le0    =   0.0
    ωe0    =   1.34297

    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb = pa.θ_e_ub - ((1.0-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    solutione = Array{Float64}(undef,Nspan,10);
    my_runge_kuttae_reverse!(solutione,y_end,espan,estep,pa,pa.θ_w_ub)

#Global Problem (Reverse)
    μ_end = solutione[1,2];
    e_end = elb;
    #Definitions for the functions of e and p:
    par   = -0.75;
    upper_cons = 0.0
    cons  = ((elb-(pa.θ_e_lb+upper_cons))/((pa.θ_w_ub-pa.θ_w_lb)^(1.0+par)));

    Nspan = 500
    y_end = [μ_end, e_end];
    xlb   = pa.θ_w_lb;
    xub   = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xub:-xstep:xlb;
    solution = Array{Float64}(undef,Nspan,2);
    my_runge_kutta_reverse!(solution,y_end,xspan,xstep,pa, par, cons)

    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,1);
    recover_controls!(controls, θspan, solution, par, cons);
    controls

    #Calculating the values of the observed μ and e (analitic solution):
    real = Array{Float64}(undef,Nspan,2);
    fill!(real, NaN);

    for i=1:Nspan
        real[i,2] = cons*(θspan[i]-pa.θ_w_lb)^(1.0+par)+(pa.θ_e_lb+upper_cons)
        real[i,1] = -1.0 + ((θspan[i]-pa.θ_w_lb)*(real[i,2]-pa.θ_e_lb))/((pa.θ_w_ub-pa.θ_w_lb)*(pa.θ_e_ub-pa.θ_e_lb))
    end

    fig, figura=plt.subplots(1,3)
    fig.suptitle("Our RK solutions")
        #μ:
    figura[1].plot(θspan[1:500], solution[1:500,2])
    figura[1].plot(θspan[1:500], real[:,2])
    figura[1].legend(["μ RK", "μ real"],loc="upper right")
    figura[1].set(ylabel="μ")
        #e:
    figura[2].plot(θspan[1:500], solution[1:500,2])
    figura[2].plot(θspan[1:500], real[:,2])
    figura[2].legend(["e RK", "e real"],loc="upper right")
    figura[2].set(ylabel="e")
        #p:
    figura[3].plot(θspan[1:500], controls[1:500,1])
    figura[3].set(ylabel="p")

    #Calculating the function through the Runge-Kutta package:
    parameters = [par,cons,pa.θ_w_lb,λe0]; #Parameters to solve the differencial equations

    initial_val = [μ_end, e_end]; #Initial values to solve the differencial equations
    θ_limits    = (pa.θ_w_ub,pa.θ_w_lb) #The θspan, to solve the equations

    prob = ODEProblem(my_runge_kutta_reverse_RKPack!,initial_val,θ_limits,parameters)
    alg = RK4()
    solution_RKPack = solve(prob,alg_hints=[:stiff])
    alg = Rosenbrock23()
    solution_RKPack = solve(prob,alg)
    #Saving the results in a matrix:
    sol_RKPackage = Array{Float64}(undef,Nspan,2);
    fill!(sol_RKPackage, NaN);
    for i = 1:Nspan
        (mu_sol,e_sol) = solution_RKPack(θspan[i]);
        sol_RKPackage[i,1] = mu_sol;
        sol_RKPackage[i,2] = e_sol;
    end
    diff_estimated = Array{Float64}(undef,Nspan,4);
    fill!(diff_estimated, NaN);
    diff_estimated[:,1] = real[:,1]-sol_RKPackage[:,1];
    diff_estimated[:,2] = real[:,1]-solution[:,1];
    diff_estimated[:,3] = real[:,2]-sol_RKPackage[:,2];
    diff_estimated[:,4] = real[:,2]-solution[:,2];

        #Viendo lo que generó el paquete:
        fig, figura=plt.subplots(1,2)
        fig.suptitle("Package RK solutions")
            #μ:
        figura[1].plot(θspan[1:500], sol_RKPackage[1:500,2])
        figura[1].plot(θspan[1:500], real[:,2])
        figura[1].legend(["μ RK", "μ real"],loc="upper right")
        figura[1].set(ylabel="μ")
            #e:
        figura[2].plot(θspan[1:500], sol_RKPackage[1:500,2])
        figura[2].plot(θspan[1:500], real[:,2])
        figura[2].legend(["e RK", "e real"],loc="upper right")
        figura[2].set(ylabel="e")

    #Viendo las diferencias con entre nuestra solución y el paquete:
    fig, figura_RK=plt.subplots(1,2)
    fig.suptitle("Runge Kutta Approximation Method")
            #μ:
        figura_RK[1].plot(θspan[:], diff_estimated[:,1])
        figura_RK[1].plot(θspan[:], diff_estimated[:,2])
        figura_RK[1].legend(["μ from package", "μ from our RK"],loc="upper right")
        figura_RK[1].set(ylabel="μ", xlabel="θw")

        figura_RK[2].plot(θspan[:], diff_estimated[:,3])
        figura_RK[2].plot(θspan[:], diff_estimated[:,4])
        figura_RK[2].legend(["e from package", "e from our RK"],loc="upper right")
        figura_RK[2].set(ylabel="e", xlabel="θw")

end  # module MySimpleModel
