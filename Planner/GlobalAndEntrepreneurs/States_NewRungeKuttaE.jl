

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

function find_statese!(du::Array{Float64},u,pa,θ,θe,ini)
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
    sse = StateE(ini[1], ini[2], ini[4], ini[6]);

    #Find optimal controls
    (ze, ne) = new_find_controlse( θ, ini[9], sse, pa);
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
