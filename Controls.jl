


function find_controls!(ctrl::Array{Control}, θ, λ, ω, ss, pa)
    h_e= pa.he( log(θ), log(ss.e));
    h_w= pa.hw( log(θ), log(ss.e));

    objective(z,n,l,p,κ) = ( ss.μ*pa.χ*l^(1+pa.ϕ)/θ + λ*( ss.e*n^pa.α*p*h_e - pa.β*z^(1+pa.σ)/(1+pa.σ)*p*h_e  - pa.χ*l^(1+pa.ϕ)*h_w/(1+pa.ϕ) ) + ω*( θ*l*h_w - n*p*h_e ) + ss.ϕ_e*p + κ*( pa.χ*l^(1+pa.ϕ)/θ  - p*n^pa.α*(1-pa.β*z^pa.σ) ) );
    # The planner can always set everthing to zero.
    mmaximum = 0.0;

    # 1 - Look for candidate optimal controls with z=0
    cnst1= pa.α*( λ*ss.uw*h_e - ss.uw^pa.ϕ*h_e - ss.ϕ_e) / (ω*h_e*(1-pa.α));
    if cnst1 > 0
        n = cnst1;
        κ = λ*ss.e*pa.α*h_e - ω*n^(pa.α-1)*h_e;

        if κ < 0
            coeff = λ*pa.χ*h_w - (ss.μ + λ*κ)*(1+pa.ψ)*pa.χ/θ;
            if coeff > 0
                l = (ω*θ*h_w / coeff)^(1/pa.ψ);
                p = (pa.χ*l^(1+pa.ψ)) / ( θ*n^pa.α );
                mmaximum = objective(0.0,n,l,p,κ)
                ctrl[1] = Control(n, p, 0.0, l);

            else
                l = Inf;
                p = Inf;
                mmaximum = Inf;
                ctrl[1] = Control(n, p, 0.0, l);
            end
        end # If κ>0, there is no solution with p=0.

    else
        cnst2 = θ*λ*h_w/(1+pa.ψ) - ss.μ;
        if cnst2 < 0.0
            p = Inf;
            κ = cnst2;
            n = eps();
            l = Inf;
            mmaximum = Inf;
            ctrl[1] = Control(n, p, 0.0, l);
        end # Else, complementary slackness for p=0 doesn't hold.
    end

    # 2 - Look for interior candidates

    cnst3 = 1.0 - pa.α*(1+pa.σ);
    N(z) = (λ*z)/(pa.σ*(1+pa.σ)) * (1-pa.α - (pa.β*z^pa.σ)*( 1/(1+pa.σ) - pa.α ) ) - param1/h_e; #cf eq 21







    vz=find_z( θ, λ, ω, ss, pa)
    (Nz,) = size(vz);



    for i=1:Nz
        z = vz[i];
        NN = λ*z/pa.σ * (1-pa.α - (pa.β*z^pa.σ)*( 1/(1+pa.σ) - pa.α ) ) - param1/h_e; #cf eq 21
        n = (NN^(1/pa.α)) / ((1-pa.α)*λ*ss.e);
        κ = λ*z*h_e/(n^pa.α);
        coeff = (ss.μ + λ*κ)*(1+pa.ψ)*pa.χ/θ - λ*pa.χ*h_w;
        if coeff > 0
            l = (ω*θ*h_w / coeff)^(1/pa.ψ);
            p = θ*n^pa.α / (pa.χ*l^(1+pa.ψ));
        else
            l = Inf;
            p = Inf;
        end

        candidate = objective(z,n,l,p,κ);
        if candidate > mmaximum
            mmaximum = candidate;
            ctrl[1] = Control(n, p, z, l);
        end
    end

    nothing
end
