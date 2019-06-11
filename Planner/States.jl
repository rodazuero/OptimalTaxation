


function dg_de(s,e,pa)
    x1 = log(s);
    x2 = log(e);

    M = (x1 - pa.μ_w)^2.0/pa.σ2_w - 2.0*pa.σ_we*(x1 - pa.μ_w)*(x2 - pa.μ_e)/(pa.σ2_w*pa.σ2_e)^0.5 + (x2 - pa.μ_e)^2.0/pa.σ2_e;
    dM_dx2 = (x2 - pa.μ_e)/pa.σ2_e - pa.σ_we*(x1 - pa.μ_w)/((pa.σ2_w*pa.σ2_e)^0.5);
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

function find_states!(du,u,pa,θ)
    uw    = u[1];
    μ     = u[2];
    e     = u[3];
    ϕ_e   = u[4];
    y_agg = u[5];
    λ     = u[6];
    l_agg = u[7];
    ω     = u[8];

    #Construct state object
    ss = State(e, uw, ϕ_e, μ, λ, ω);
    #Find optimal controls
    (z, n, l, p) = new_find_controls( θ, ss, pa);

    h_e=  pa.he( θ, e);
    h_w=  pa.hw( θ, e);
    h_tot= h_w + p*h_e;

    Vw = uw^pa.ϕ + (θ*l*ω - λ*(uw + pa.χ*(l^(1+pa.ψ))/(1+pa.ψ) ));
    Ve = uw^pa.ϕ + λ*e*n^pa.α - λ*pa.β*(z^(1+pa.σ))/(1+pa.σ) - ω*n - λ*uw;
    dhw_de = pa.gg( θ , e);
    dhe_de=  integrate_dg_de(θ,e,pa);
    dVe_de= λ*n^pa.α;


    du[1] = pa.χ*(l^(1+pa.ψ))/θ
    if uw > 0.0
        du[2] = (λ - pa.ϕ*uw^(pa.ϕ-1))*h_tot;
    else
        du[2] = Inf;
    end
    du[3] = p;
    du[4] = -( Vw*dhw_de + dVe_de*p*h_e + Ve*p*dhe_de );
    du[5] = ( e*(n^pa.α) - pa.β*(z^(1+pa.σ) )/(1+pa.σ) )*p*h_e - uw*h_tot - pa.χ*(l^(1+pa.ψ))/(1+pa.ψ)*h_w;
    du[6] = 0.0;
    du[7] = θ*l*h_w - (n-pa.ϵ)*p*h_e;
    du[8] = 0.0;

    nothing
end
