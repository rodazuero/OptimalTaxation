function find_states!(du::Array{Float64,1},u::Array{Float64,1},pa,θ::Float64,controls::Array{Float64,1})
    #Definitions we are using:
    uw     = u[1];
    μ      = u[2];
    e      = u[3];
    ϕe     = u[4];
    y_agg  = u[5];
    λ      = u[6];
    lf_agg = u[7];
    ωf     = u[8];
    li_agg = u[9];
    ωi     = u[10];
    wi     = u[11];
    ϕw     = u[12];

    lf_new = u[13];
    li_new = u[14];
    y_new  = u[15];

    z  = controls[1];
    n  = controls[2];
    ni = controls[3];
    l  = controls[4];
    li = controls[5];
    p  = controls[6];

    h_w = pa.hw(θ, e);
    h_e = pa.he(θ, e);
    hh  = pa.hh(θ, e, p);

    #Value of entrepreneurs and workers:
    Vw = pa.indicator*uw^pa.ϕ - λ*uw + θ*l*ωf + (ωi-ωf)*θ*li - λ*( pa.χ/(1.0+pa.ψ)*l^(1.0+pa.ψ) + pa.κ/(1.0+pa.ρ)*(θ*li)^(1.0+pa.ρ));
    Ve = pa.indicator*uw^pa.ϕ - λ*uw - ωf*n - (ωi-ωf)*ni + λ*( e*n^pa.α - pa.δ/(1.0+pa.γ)*ni^(1.0+pa.γ) - pa.β/(1.0+pa.σ)*z^(1.0+pa.σ));
    #Uniform distribution:
    ∂he_∂e = 0.0;
    ∂hw_∂e = pa.gg(θ,e);

    du[1] = pa.χ/θ*l^(1.0+pa.ψ);
    if uw > 0.0
        du[2] = (λ - pa.indicator*pa.ϕ*uw^(pa.ϕ-1.0))*hh;
    else
        du[2] = Inf;
    end
    du[3] = p;
    du[4] = - ( ∂he_∂e*Ve*p + ∂hw_∂e*Vw + λ*n^pa.α*(1.0-pa.α/pa.γ*ni/n)*p*h_e);
    du[5] = (e*n^pa.α - pa.δ/(1.0+pa.γ)*ni^(1.0+pa.γ) - pa.β/(1+pa.σ)*z^(1+pa.σ) - uw)*p*h_e -
            (uw + pa.χ/(1.0+pa.ψ)*l^(1.0+pa.ψ) + pa.κ/(1.0+pa.ρ)*(θ*li)^(1.0+pa.ρ))*h_w;
    du[6] = 0.0;
    du[7] = θ*(l-li)*h_w - (n-ni)*p*h_e;
    du[8] = 0.0;
    du[9] = θ*li*h_w - ni*p*h_e;
    du[10] = 0.0;
    du[11] = 0.0;
    du[12] = ss.λ/pa.ρ*(ni*p*h_e-li*θ*h_w);

    du[13] = θ*(l-li)*h_w + (n-ni)*p*h_e;
    du[14] = θ*li*h_w + ni*p*h_e;
    du[15] = ( e*n^pa.α - pa.β*z^(1+pa.σ)/(1+pa.σ) - pa.δ/(1.0+pa.γ)*ni^(1.0+pa.γ) )*p*h_e; #Production - evasion - informality;

    nothing
end
