


function find_controls( θ, ss, pa)
    h_e= pa.he(θ, ss.e);
    h_w= pa.hw(θ, ss.e);

    nn = NaN;
    zz = NaN;
    ll = NaN;
    pp = NaN;

    A= ((ss.uw^pa.ϕ-ss.λ*ss.uw)*h_e + ss.ϕ_e) / (ss.λ*h_e);

    aux(z,n) =  ss.λ*pa.χ*h_w - (ss.μ + z*ss.λ*h_e / (pa.σ* n^pa.α) )*(1+pa.ψ)*pa.χ/θ;
    l(coef) = coef > 0.0 ? (ss.ω*θ*h_w/coef)^(1.0/pa.ψ) : Inf;
    lz_extremum(coef) = coef > 0.0 ? (ss.ω*θ*h_w/(coef*pa.χ))^(1.0/pa.ψ) : Inf;
    p(l,z,n) = pa.χ*l^(1.0+pa.ψ) / (θ*n^pa.α*(1.0-pa.β*z^pa.σ));
    z_aux(n) = z(n,ss,pa, θ);
    z_ub(n) =  z(n,ss,pa, θ) - (1.0 / pa.β)^(1.0/pa.σ);
    zprime(n) = pa.α*pa.σ*ss.e*n^(pa.α-1.0) - ((pa.α*(1.0+pa.σ)-1.0)/pa.α)*ss.ω/ss.λ;
    eq6_aux(n) = eq6(n,ss,pa, θ)

#Compute the bounds for n such that z(n) is positive and less beta z ^sigma < 1.


    nub= ( (1.0/pa.β)^((1.0+pa.σ)/pa.σ)*(pa.β/(1.0+pa.σ)) - A )*(ss.λ*pa.α /(ss.ω*(1.0-pa.α))) ;
    auxz1 = pa.χ*h_w-(1+pa.ψ)*(ss.μ + ss.λ/pa.σ*( (1/pa.β)^(1/pa.σ) - 1e-10)/nub^pa.α*h_e  )/θ - (1 + pa.ψ)/θ * (ss.λ*ss.e- ss.ω/pa.α *nub^(1-pa.α))*h_e;

    if nub < 0.0
        return     nn, zz, ll, pp;
    elseif A > 0.0 # z(0)= (1+sigma)*A > 0 <=> A > 0
        nlb= 0.0;
    else
        nlb = BigFloat(- A*ss.λ*pa.α / (ss.ω*(1.0-pa.α)));
        auxz0 = pa.χ*h_w-(1+pa.ψ)*ss.μ/θ - (1 + pa.ψ)/θ * (ss.λ*ss.e- ss.ω/pa.α *nlb^(1-pa.α))*h_e;
    end

    candidate_nn= my_find_zeros(eq6_aux, BigFloat(nlb), BigFloat(nub), 10);

    if A < 0.0
        append!(candidate_nn,nlb) # Case z=0
    end
    append!(candidate_nn,nub)

    (K,)=size(candidate_nn)

    if K==0
        return nn, zz, ll, pp;
    else
        objective(z,n,l,p) = ss.uw^pa.ϕ*(h_w + p*h_e) + ss.μ*pa.χ*l^(1.0+pa.ψ)/θ + ss.λ*( ss.e*n^pa.α*p*h_e - pa.β*z^(1.0+pa.σ)/(1.0+pa.σ)*p*h_e - ss.uw*(h_w + p*h_e) - pa.χ*l^(1+pa.ψ)*h_w/(1+pa. ψ) ) + ss.ω*( θ*l*h_w - n*p*h_e ) + ss.ϕ_e*p ;
        nn= candidate_nn[1];

        zz= z_aux(nn);

        if zz == 0.0
            ll = lz_extremum(auxz0);
        elseif abs(zz^pa.σ*pa.β - 1.0) < 1e-14
            println(1.0-pa.β*zz^pa.σ)
            ll = lz_extremum(auxz1);
        else
            ll = l(aux(zz, nn));
        end
        pp= p(ll, zz, nn);
        potential_value = objective(zz,nn,ll,pp)
        current_max = isnan(potential_value) ? -Inf : potential_value;
#        println("objetive = ", current_max)
#        println("profits= ", Float64(ss.λ*ss.e*nn^pa.α-ss.ω*nn), " n = ", Float64(nn), " z = ", Float64(zz), " l = ", Float64(ll), " p= ", Float64(pp) )

        for j=2:K
            potential_n=candidate_nn[j];
            potential_z= z_aux(potential_n);
            if potential_z==0.0
                potential_l = lz_extremum(auxz0);
            elseif abs((potential_z + 1e-10)^pa.σ*pa.β - 1.0) < 1e-14
                potential_l = lz_extremum(auxz1);
                println(1.0-pa.β*potential_z^pa.σ)
            else
                potential_l = l(aux(potential_z, potential_n));
            end
            potential_p= p(potential_l, potential_z, potential_n);

            potential_value= objective(potential_z,potential_n,potential_l,potential_p);
#            println("objetive = ", potential_value)
#            println("profits= ", Float64(ss.λ*ss.e*potential_n^pa.α-ss.ω*potential_n), " n = ", Float64(potential_n), " z = ", Float64(potential_z), " l = ", Float64(potential_l), " p= ", Float64(potential_p) )

            if isnan(potential_value)

            elseif potential_value > current_max
                current_max = potential_value;
                nn= potential_n;
                zz= potential_z;
                ll= potential_l;
                pp= potential_p
            end
        end
    end
    println(" A= ", A,  " he = ", h_e)
    println("n= ", Float64(nn), " z= ", Float64(zz), " l= ",Float64(ll), " p= ",Float64(pp), " denominator= ", aux(zz,nn))
    println("profits= ", Float64(ss.λ*ss.e*nn^pa.α-ss.ω*nn), " value = ", Float64(current_max))
    nn, zz, ll, pp;
end



function init_controls( uw0, e0, ϕ_e0, λ0, ω0, pa)
    θ      = exp(pa.μ_w-3*pa.σ2_w);
    μ0     = 0.0;
    y_agg0 = 0.0;
    l_agg0 = 0.0;

    ss = State(e0, uw0, ϕ_e0, μ0, λ0, ω0);
    find_controls( θ, ss, pa);
end


function  eq6(n, ss, pa, θ)
    if n == Inf
        eq6 = Inf;
    else
        h_e = pa.he(θ, ss.e);
        A= ((ss.uw^pa.ϕ-ss.λ*ss.uw)*h_e + ss.ϕ_e) / (ss.λ*h_e);
        eq6= ss.e*n^pa.α - z(n, ss, pa, θ)/pa.σ + pa.β/((1+pa.σ)*pa.σ)*z(n, ss, pa, θ)^(pa.σ+1) - ss.ω/ss.λ*n + A;
    end
    if eq6 <0.0 && eq6>-1e-15
        eq6 = 0.0
    end
    eq6
 end


function  z(n, ss, pa, θ)
    h_e =  pa.he(θ, ss.e);
    A= ((ss.uw^pa.ϕ-ss.λ*ss.uw)*h_e + ss.ϕ_e) / (ss.λ*h_e);

    base_z = n*ss.ω*(1.0-pa.α)/pa.α/ss.λ + A

    if base_z< 1e-15*n && base_z>-1e-15*n
        return z=0.0
    else
        z = ( (1.0+pa.σ)/pa.β* base_z )^(1.0/(1.0+pa.σ))
        if abs(z^(pa.σ)*pa.β - 1.0) < 1e-15*n
            return z =(1.0/pa.β)^(1.0/pa.σ) - 1e-10
        end
    end
    z
end
