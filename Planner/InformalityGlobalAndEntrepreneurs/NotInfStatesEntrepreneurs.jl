function NotInffind_statese!(du::Array{Float64},u,pa,θ,θe,ini)
    ue     = u[1];
    μe     = u[2];
    ye_agg = u[3];
    λe     = u[4];
    le_agg = u[5];
    ωe     = u[6];
    le_new = u[7];
    ye_new = u[8];

    #Construct state object
    #sse = State(ue, μe, λe, ωe);
    sse = NotInfStateE(ini[1], ini[2], ini[4], ini[6]);

    #Find optimal controls
    (ze, ne) = NotInfnew_find_controlse( θ, ini[9], sse, pa);
    any(isnan,(ze, ne)) && error("Function find_statese gets NaN controls")

    h_e=  pa.he(θ, θe);

    du[1] = ne^pa.α*(1.0 - pa.β*ze^pa.σ)
    if ue > 0.0
        du[2] = (λe - pa.indicator*pa.ϕ*ue^(pa.ϕ-1))*h_e;
    else
        du[2] = Inf;
    end
    du[3] = (θe*ne^pa.α- (pa.β/(1+pa.σ))*ze^(1+pa.σ) - ue)*h_e;
    du[4] = 0.0;
    du[5] = -ne*h_e;
    du[6] = 0.0;
    du[7] = ne*h_e;
    du[8] = (θe*ne^pa.α- (pa.β/(1+pa.σ))*ze^(1+pa.σ) - ue)*h_e;

    nothing
end
