function find_states!(du::Array{Float64,1},u::Array{Float64,1},pa,θ::Float64,ini::Array{Float64,1})
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

    #Construct state object
    ss = State(ini[1], ini[2], ini[3], ini[4], ini[5], ini[6], ini[7], ini[8],
               ini[9], ini[10], ini[11], ini[12]);

    #Find optimal controls
    (z, n, ni, l, li, p) = new_find_controls(θ, ini[end], sse, pa);
    any(isnan,(z, n, ni, l, li, p)) && error("Function find_states gets NaN controls.")

    h_w = pa.hw(θ, e);
    h_e = pa.he(θ, e);
    hh  = pa.hh(θ, e, p);

    #Value of entrepreneurs and workers:
    Vw = pa.indicator*uw^pa.ϕ - λ*uw + θ*l*ωf + (ωi-ωf)*θ*li - λ*( pa.χ/(1.0+pa.ψ)*l^(1.0+pa.ψ) + pa.κ/(1.0+pa.ρ)*(θ*li)^(1.0+pa.ρ));
    Ve = pa.indicator*uw^pa.ϕ - λ*uw - ωf*n - (ωi-ωf)*ni + λ*( e*n^pa.α - pa.δ/(1.0+pa.γ)*ni^(1.0+pa.γ) - pa.β/(1.0+pa.σ)*z^(1.0+pa.σ));
    ∂ni_∂e = pa.α/(pa.δ*pa.γ)*ni^(1.0-pa.γ)/n^(1.0-pa.α);
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
    du[4] = - ( ∂he_∂e*Ve*p + ∂hw_∂e*Vw + λ*n^pa.α*p*h_e - ∂ni_∂e*(λ*δ*ni^pa.γ + ωi - ωf)*p*h_e );
    du[5] = (e*n^pa.α - pa.δ/(1.0+pa.γ)*ni^(1.0+pa.γ) - pa.β/(1+pa.σ)*z^(1+pa.σ) - uw)*p*h_e -
            (uw + pa.χ/(1.0+pa.ψ)*l^(1.0+pa.ψ) + pa.κ/(1.0+pa.ρ)*(θ*li)^(1.0+pa.ρ))*h_w;
    du[6] = 0.0;
    du[7] = θ*(l-li)*h_w - (n-ni)*p*h_e;
    du[8] = 0.0;
    du[9] = θ*li*h_w - ni*p*h_e;
    du[10] = 0.0;
    du[11] = 0.0;
    du[12] = - (li^(1.0-pa.ρ)/(pa.ρ*pa.κ*θ^pa.ρ)*(ωi - ωf - λ*pa.κ*(θ*li)^pa.ρ)*θ*h_w +
                ni^(1.0-pa.γ)/(pa.δ*pa.γ)*(ωi - ωf + λ*pa.δ*ni^pa.γ)*p*h_e);

    du[13] = θ*(l-li)*h_w + (n-ni)*p*h_e;
    du[14] = θ*li*h_w + ni*p*h_e;
    du[15] = ( e*n^pa.α - pa.β*z^(1+pa.σ)/(1+pa.σ) - pa.δ/(1.0+pa.γ)*ni^(1.0+pa.γ) )*p*h_e; #Production - evasion - informality;

    nothing
end
