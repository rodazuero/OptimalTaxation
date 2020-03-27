function find_statese!(du::Array{Float64,1},u::Array{Float64,1},pa,θe::Float64,controls::Array{Float64,1})
    #Definitions we are using:
    ue      = u[1];
    μe      = u[2];
    ye_agg  = u[3];
    λe      = u[4];
    lfe_agg = u[5];
    ωfe     = u[6];
    lie_agg = u[7];
    ωie     = u[8];
    wie     = u[9];
    ϕwe     = u[10];

    lfe_new = u[11];
    lie_new = u[12];
    ye_new  = u[13];

    ze  = controls[1];
    ne  = controls[2];
    nie = controls[3];

    h_e = pa.he(pa.θ_w_ub, θe);

    #The derivative vector:
    du[1] = ne^pa.α*(1.0 - pa.β*ze^pa.σ)
    println("ue' = ", du[1])
    if ue > 0.0
        du[2] = (λe - pa.indicator*pa.ϕ*ue^(pa.ϕ-1.0))*h_e;
    else
        du[2] = Inf;
    end
    #du[3] = (θe*ne^pa.α - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - pa.β/(1+pa.σ)*ze^(1+pa.σ) - ue)*h_e;
    du[3] = (θe*ne^pa.α - pa.β/(1+pa.σ)*ze^(1+pa.σ) - ue)*h_e;
    du[4] = 0.0;
    #du[5] = -(ne-nie)*h_e;
    du[5] = -(ne)*h_e;
    du[6] = 0.0;
    du[7] = -nie*h_e;
    du[8] = 0.0;
    du[9] = 0.0;
    du[10] = -1.0/(pa.δ*pa.γ)*nie^(1.0-pa.γ)*(λe*pa.δ*nie^pa.γ-ωfe+ωie)*h_e;

    du[11] = (ne-nie)*h_e;
    du[12] = nie*h_e;
    du[13] = ( θe*ne^pa.α - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - pa.β/(1+pa.σ)*ze^(1+pa.σ) )*h_e;

    #println(du)
    nothing
end
