cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\InformalityGlobalAndEntrepreneurs")
cd("C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\InformalityGlobalAndEntrepreneurs")

using Roots
using NLopt
using Statistics
#using PyPlot
using DataFrames
using CSV
using NLsolve
using Plots
#using DifferentialEquations

#Global Problem
include("Definitions.jl")

pa = init_parameters();

#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    gp     =   0.993

    ue    =   640.0 #guess
    μe    =   0.0 - 1.0e-10
    ye    =   0.0
    λe    =   1.0
    le    =   0.0
    ωfe   =   1.343 #guess
    lie   =   0.0
    ωie   =   1.343 #guess
    wie   =   1.343 #guess
    ϕwe   =   0.0
    ye_agg=   0.0
    le_agg=   0.0
    lie_agg=  0.0


    #Construct state object
    sse = StateE(ue, μe, ye_agg, λe, le_agg, ωfe, lie_agg, ωie, wie, ϕwe);

    θe = pa.θ_e_ub
    θ = pa.θ_w_ub

    #Recover ditributions
    h_e= pa.he(θ, θe);

    #Defining bounds we use in various cases (limits of z):
    n_full_info = ((sse.λe*pa.α*θe)/sse.ωfe)^(1.0/(1.0-pa.α));
    z_max   = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    z_lwbar = eps(); #The smallest possible z.
    z_1     = z_max;
    z_2     = -pa.σ*sse.μe*n_full_info^pa.α/(sse.λe*h_e);
    #Keep the lower number:
    z_upbar = min(z_1,z_2);
    z_upbar = z_2

    #With informality
        #Defining the functions we are using:
        fun_ni(n)   = ((pa.α*θe*n^(pa.α-1.0)-sse.wie)/pa.δ)^(1.0/pa.γ);
        n_opt(z)    = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
        fun_nz(z,n) = (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωfe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ)+
                       fun_ni(n)^(1.0-pa.γ)/n^(2.0-pa.α)*pa.α*(1.0-pa.α)*θe/(pa.δ*pa.γ)*(sse.λe*pa.δ*fun_ni(n)^pa.γ+sse.ωie-sse.ωfe)*h_e);
        fun_z(z)    = fun_nz(z,n_opt(z));

    #Without informality
        #Defining the functions we are using:
        n_optw(z)   = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
        fun_nzw(z,n) = (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωfe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ));
        fun_zw(z)  = fun_nzw(z,n_optw(z));


    #1.Grid
    Nz = 1000; #Number of Zs
    zstep = (z_upbar - z_lwbar)/(Nz-1);
    zgrid = range(z_lwbar, stop = z_upbar, length = Nz);
    #println("Number of z = ", Nz)

    grid_z = Array{Float64}(undef,Nz,6);
    fill!(grid_z,NaN);

    #1. Grid
    iteration=NaN
    iterationw=NaN

    for j = 1:Nz

        println(j)

        potential_z  = collect(zgrid)[j];
        potential_n  = n_opt(potential_z);
        potential_cpo = fun_z(potential_z);
        potential_zw  = collect(zgrid)[j];
        potential_nw  = n_optw(potential_zw);
        potential_cpow = fun_zw(potential_zw);

        grid_z[j,1] = potential_z
        grid_z[j,2] = potential_n
        grid_z[j,3] = potential_cpo
        grid_z[j,4] = potential_zw
        grid_z[j,5] = potential_nw
        grid_z[j,6] = potential_cpow

        if  j>1 && grid_z[j,3]*grid_z[j-1,3]<0
            iteration=j
        end

        if  j>1 && grid_z[j,6]*grid_z[j-1,6]<0
            iterationw=j
        end

    end


    if iteration=NaN
        iteration=Nz
    end

    if iterationw=NaN
        iterationw=Nz
    end


    new_lwbar=grid_z[iteration-1,1]
    new_upbar=grid_z[iteration,1]

    new_lwbarw=grid_z[iteration-1,4]
    new_upbarw=grid_z[iteration,4]


    if fun_z(new_lwbar)*fun_z(new_upbar)<0
        potential_z = find_zero(fun_z, (new_lwbar,new_upbar), Bisection());
    else
        error("Bisection: signs equal --> Cannot solve.")
    end


    if fun_zw(new_lwbarw)*fun_zw(new_upbarw)<0
        potential_z = find_zero(fun_zw, (new_lwbarw,new_upbarw), Bisection());
    else
        error("Bisection: signs equal --> Cannot solve.")
    end



            if  j>1 && potential_cpo*grid_z[(j-1),3]<0
                iteration=j
            end

            if  j>1 && potential_cpow*grid_z[(j-1),6]<0
                iterationw=j
            end
