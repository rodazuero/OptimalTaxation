function find_statese!(du::Array{Float64,1},u::Array{Float64,1},pa,θ::Float64,θe::Float64,ini::Array{Float64,1})
    ue      = u[1];
    μe      = u[2];
    ye_agg  = u[3];
    λe      = u[4];
    le_agg  = u[5];
    ωfe     = u[6];
    lie_agg = u[7];
    ωie     = u[8];
    wie     = u[9];
    ϕie     = u[10];
    
    le_new  = u[11];
    lie_new = u[12];
    ye_new  = u[13];

    #Construct state object
    sse = StateE(ini[1], ini[2], ini[3], ini[4], ini[5], ini[6], ini[7], ini[8], ini[9], ini[10]);

    #Find optimal controls
    (ze, ne, nie) = new_find_controlse(θ, ini[end], sse, pa);
    any(isnan,(ze, ne, nie)) && error("Function find_statese gets NaN controls.")

    h_e = pa.he(θ, θe);

    du[1] = ne^pa.α*(1.0 - pa.β*ze^pa.σ)
    if ue > 0.0
        du[2] = (λe - pa.indicator*pa.ϕ*ue^(pa.ϕ-1.0))*h_e;
    else
        du[2] = Inf;
    end
    du[3] = (θe*ne^pa.α - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - pa.β/(1+pa.σ)*ze^(1+pa.σ) - ue)*h_e;
    du[4] = 0.0;
    du[5] = -(ne-nie)*h_e;
    du[6] = 0.0;
    du[7] = -nie*h_e;
    du[8] = 0.0;
    du[9] = 0.0;
    du[10] = 0.0;

    du[11] = (ne-nie)*h_e;
    du[12] = nie*h_e;
    du[13] = (θe*ne^pa.α - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - pa.β/(1+pa.σ)*ze^(1+pa.σ) - ue)*h_e;

    nothing
end
