


function dg_de(θ,e,pa)
    x1 = log(θ);
    x2 = log(e)
    z = (x1 - pa.μ_w)^2.0/pa.σ2_w - 2.0*pa.σ_we*(x1 - pa.μ_w)*(x2 - pa.μ_e)/(pa.σ2_w*pa.σ2_e)^0.5 + (x2 - pa.μ_e)^2.0/pa.σ2_e;
    dz_dx2 = 2*(x2 - pa.μ_e)/pa.σ2_e - 2*pa.σ_we*(x1 - pa.μ_w)/((pa.σ2_w*pa.σ2_e)^0.5);
    dnorm_dx2 = 1.0/(2.0*π* (pa.σ2_w*pa.σ2_e)^0.5 * (1.0 - pa.σ_we^2.0)^0.5) * exp(-z/(2.0*(1.0 - pa.σ_we^2.0)))*( -1.0/( 2.0*(1.0 - pa.σ_we^2.0)) )*dz_dx2;
    dnorm_dx2/e
end

function integrate_dg_de(θ, e, pa)
    (N,)=size(pa.nodes);
    integral=0;
    for i=1:N
        s = (pa.nodes[i] +1)*θ/2;
        integral += dg_de(s,e,pa)*pa.weights[i];
    end
    integral = integral*θ/2;
end

function find_states!(du,u,pa,θ)
    uw    = u[1];
    μ     = u[2];
    e     = u[3];
    ϕ_e   = u[4]
    y_agg = u[5];
    λ     = u[6];
    l_agg = u[7];
    ω     = u[8];

    ss = State(e, uw, ϕ_e, μ);

    ctrl = Array{Control}(undef,1);
    find_controls!(ctrl, θ, λ, ω, ss, pa);
    cc = ctrl[1];

    h_e= pa.he( log(θ), log(ss.e));
    h_w= pa.hw( log(θ), log(ss.e));
    h_tot= h_w + cc.p*h_e;

    Vw = uw + (θ*ss.l*ω - λ*(uw + pa.χ*(cc.l^(1+pa.ψ))/(1+pa.ψ) ));
    Ve = uw + λ*e*cc.n^pa.α - λ*pa.β*(z^(1+pa.σ))/(1+pa.σ) - ω*cc.n - λ*uw;
    dhw_de = pa.gg( log(θ) , log(e));
    dhe_de=  integrate_dg_de(θ,e,pa);
    dVe_de= λ*n^pa.α;



    du[1] = pa.χ*(cc.l^(1+pa.ψ))
    if uw > 0.0
        du[2] = (λ - pa.ϕ*uw^(pa.ϕ-1))*h_tot;
    else
        du[2] = Inf;
    end
    du[3] = cc.p;
    du[4] = -( Vw*dhw_de + dVe_de*p*h_e + Ve*p*dhe_de );
    du[5] = ( e*(cc.n^α) - pa.β*(cc.z^(1+pa.σ) )/(1+pa.σ) )*cc.p*h_e -uw*h_tot - pa.χ*(cc.l^(1+pa.ψ))/(1+pa.ψ)*h_w;
    du[6] = 0.0;
    du[7] = θ*cc.l*h_w - cc.n*cc.p*h_e;
    du[8] = 0.0;

    nothing
end
