

function dg_de(s,e,pa)
    x1 = log(s);
    x2 = log(e);

    M = (x1 - pa.μ_w)^2.0/pa.σ2_w - 2.0*pa.σ_we*(x1 - pa.μ_w)*(x2 - pa.μ_e)/(pa.σ2_w*pa.σ2_e)^0.5 + (x2 - pa.μ_e)^2.0/pa.σ2_e;
    dM_dx2 = 2.0*((x2 - pa.μ_e)/pa.σ2_e - pa.σ_we*(x1 - pa.μ_w)/((pa.σ2_w*pa.σ2_e)^0.5));
    dnorm_dx2 = 1.0/(2.0*π* (pa.σ2_w*pa.σ2_e)^0.5 * (1.0 - pa.σ_we^2.0)^0.5) * (1.0/s) * exp(-M/(2.0*(1.0 - pa.σ_we^2.0)))*( -1.0/(1.0 - pa.σ_we^2.0) )*dM_dx2;
    norm_biv =  1.0/(2.0*π* (pa.σ2_w*pa.σ2_e)^0.5 * (1.0 - pa.σ_we^2.0)^0.5) * (1.0/s) * exp(-M/(2.0*(1.0 - pa.σ_we^2.0)));


    dg_de= (dnorm_dx2 - norm_biv)/(e^2.0)
end

function integrate_dg_de(θ, e, pa)
    (N,)=size(pa.nodes);
    b = θ;


    integral=0.0;
    for i=1:N
        x = (pa.nodes[i] +1.0)*b/2.0;
        integral += dg_de(x,e,pa)*pa.weights[i];
    end
    integral = integral*b/2.0;
end

function dhe(s,e,pa)
    #dhe=-(pdf(dist_marginal_e, log(e))*(1+((ln(log(e))-pa.μ_e)/pa.σ2_e)))*cdf(dist_w_given_e(log(e)), log(s) )*(1/log(e))
    dhe=-(pdf(Normal(pa.μ_e,((pa.σ2_e)^0.5)), log(e))*(1/e)*(1+((log(e)-pa.μ_e)/pa.σ2_e)))*cdf(Normal(pa.μ_w + pa.σ_we*((pa.σ2_w)^0.5)/((pa.σ2_e)^0.5)*(log(e)-pa.μ_e),((1.0-pa.σ_we^2.0)*pa.σ2_w)^0.5), log(s))*(1/e)
end

function find_states!(du,u,pa,θ,ini)
    uw    = u[1];
    μ     = u[2];
    e     = u[3];
    ϕ_e   = u[4];
    y_agg = u[5];
    λ     = u[6];
    l_agg = u[7];
    ω     = u[8];
    l_new = u[9];
    y_new = u[10];

    #Construct state object
    #ss = State(e, uw, ϕ_e, μ, λ, ω);
    ss = State(ini[3], ini[1], ini[4], ini[2], ini[6], ini[8]);

    #Find optimal controls
    #(z, n, l, p) = new_find_controls( ini[11], ss, pa);
    (z, n, l, p, κnz) = new_find_controls( ini[11], ss, pa);
    println("z = ", z, "n = ", n, "l = ", l, "p = ", p)
    any(isnan,(z, n, l, p)) && error("Function find_states gets NaN controls")

    h_e=  pa.he( θ, e);
    h_w=  pa.hw( θ, e);
    h_tot= h_w + p*h_e;

    Vw = pa.indicator*uw^pa.ϕ + θ*l*ω - λ*uw - λ*pa.χ*l^(1.0+pa.ψ)/(1.0+pa.ψ);
    Ve = pa.indicator*uw^pa.ϕ + λ*e*n^pa.α - λ*pa.β*z^(1.0+pa.σ)/(1.0+pa.σ) - ω*n - λ*uw;
    #Non independent distributions
    #dhw_de = pa.gg( θ , e);
    #dhe_de=  integrate_dg_de(θ,e,pa); #Non independent distributions
    #=Independent distributions
    dhe_de= dhe(θ,e,pa)
    dhw_de= (pdf(Normal(pa.μ_e,((pa.σ2_e)^0.5)), log(e))*(1/e))*((pdf(Normal(pa.μ_w,((pa.σ2_w)^0.5)), log(θ)))*(1/θ))
    #dVe_de= λ*n^pa.α;=#
    #Uniform distribution
    ∂he_de = 0.0;
    ∂hw_de = pa.gg(θ,e);
    #∂hw_de = pa.gg(θ,e); #In the uniform case.

    du[1] = pa.χ*l^(1.0+pa.ψ)/θ;
    if uw > 0.0
        du[2] = (λ - pa.indicator*pa.ϕ*uw^(pa.ϕ-1.0))*h_tot;
    else
        du[2] = Inf;
    end
    du[3] = p;
    du[4] = -( Vw*∂hw_de + Ve*p*∂he_de + λ*n^pa.α*(p*h_e + κnz));
    du[5] = ( e*n^pa.α - pa.β*z^(1.0+pa.σ)/(1.0+pa.σ) )*p*h_e - uw*h_tot - pa.χ*l^(1.0+pa.ψ)/(1.0+pa.ψ)*h_w;
    du[6] = 0.0;
    du[7] = θ*l*h_w - (n-pa.ϵ)*p*h_e;
    du[8] = 0.0;
    du[9] = θ*l*h_w + (n-pa.ϵ)*p*h_e;
    du[10]= e*n^pa.α*p*h_e - pa.β*z^(1+pa.σ)/(1+pa.σ)*p*h_e; #Production - evasion

    println("ϕ_e' =", du[4], "u_w' =", du[1])

    nothing
end
