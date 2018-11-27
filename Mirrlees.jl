function mirrlees!(du,u,param,θ)
    uw = u[1];
    μ = u[2];
    y = u[3];
    λ = u[4];
    (ω, pa)= param;

    μ_w = 1.7626;
    σ2_w = 1.0921;
    h_w =pdf( Normal(μ_w,σ2_w), log(θ));

    coeff = λ*pa.χ*h_w - μ*(1+pa.ψ)*pa.χ/θ;
    if coeff > 0
        l = (λ*ω*θ*h_w/coeff)^(1/pa.ψ);
    else
        l = Inf;
    end


    du[1] = pa.χ*(l^(1+pa.ψ))/θ;
    if uw > 0.0
        du[2] = (λ - pa.ϕ*uw^(pa.ϕ-1))*h_w;
    else
        println("pilas")
        du[2] = Inf;
    end
    du[3] = ω*θ*l - uw - (pa.χ*l^(1+pa.ψ))/(1+pa.ψ)
    du[4] = 0.0;


    nothing
end

function bc1!(residual, u,param,θ)
    residual[1] = u[1][2]  # the multiplier at θ_lb should be zero
    residual[2] = u[end][2]   # the multiplier at θ_ub should be zero
    residual[3] = u[1][3]  # Lower boundary for Y
    residual[4] = u[end][2] - pa.G   # Upper boundary for Y

end

function solvediffeq!(sol::Array, u0, tspan, ω, pa)
    param = (ω, pa)

    prob = ODEProblem(mirrlees!,u0,tspan,param)
    sol[1] = solve(prob)
    sol[1].u[end][2], sol[1].u[end][3]
end


function solvepartialeq!(sol::Array, bounds, u0, tspan, λ, ω, pa)

    param = (ω, pa);
    bvpobjective(x) = solvediffeq!(sol, [x,0.0,0.0,λ], tspan, ω, pa)[1];
    uw0=find_zero(bvpobjective, bounds);
    u0=[uw0-1000*eps(),0.0,0.0,λ];
    prob = ODEProblem(mirrlees!, u0, tspan, param)
    sol[1] = solve(prob)
    nothing
end

function recover(θ, sol::OrdinaryDiffEq.ODECompositeSolution, λ, ω, pa)
    uw = sol(θ)[1];
    μ  = sol(θ)[2];
    y  = sol(θ)[3];
    λ  = sol(θ)[4];

    h_w =pdf( Normal(pa.μ_w, pa.σ2_w), log(θ));

    coeff = λ*pa.χ*h_w - μ*(1+pa.ψ)*pa.χ/θ;
        if coeff > 0
            l = (λ*ω*θ*h_w/coeff)^(1/pa.ψ);
        else
            l = Inf;
        end

    marginaltax = (ω*θ-pa.χ*(l^pa.ψ))/ω*θ;
    tax = ω*θ*l - uw - (pa.χ*l^(1+pa.ψ))/(1+pa.ψ);
    base = ω*θ*l;
    labor= l;

    marginaltax, tax, base, labor
end
