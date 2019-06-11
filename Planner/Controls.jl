function new_find_controls( θ, ss, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e= pa.he(θ, ss.e);
    h_w= pa.hw(θ, ss.e);

    #omega/lambda is often used
    ww = ss.ω/ss.λ;


    #Define A
    A= ((ss.uw^pa.ϕ-ss.λ*ss.uw)*h_e + ss.ϕ_e) / (ss.λ*h_e);

    #Define bounds for z, and make sure there exist feasible values for n and z
    z_ub = (1.0/pa.β)^(1.0/pa.σ)
    pa.β*(z_ub^(1.0+pa.σ))/(1.0+pa.σ) < A &&  error("n<0,adjust A downward")

    if A<0.0
        z_lb = eps()
    else
        z_lb = ( (1.0+pa.σ)*A/pa.β )^(1/(1.0+pa.σ))
    end
    pa.β*(z_ub^(1.0+pa.σ))/(1.0+pa.σ) < A &&  error("no range for z, adjust A downward")


    #Define grid for z
    #grid size - defining the grid inside the function is inefficient, can be adjusted if need to speed up
    Nz = 500;
    zstep = (z_ub-z_lb)/(Nz-1);
    zgrid = z_lb:zstep:z_ub;

    #Define the hamiltonian as function of controls given states
    objective(z,n,l,p) = ss.uw^pa.ϕ*(h_w + p*h_e) + ss.μ*pa.χ*l^(1.0+pa.ψ)/θ + ss.λ*( ss.e*n^pa.α*p*h_e - pa.β*z^(1.0+pa.σ)/(1.0+pa.σ)*p*h_e - ss.uw*(h_w + p*h_e) - pa.χ*l^(1+pa.ψ)*h_w/(1+pa. ψ) ) + ss.ω*( θ*l*h_w - n*p*h_e ) + ss.ϕ_e*p ;
    #Define function n(z)
    n_fun(z) = pa.α /( ww* (1.0 - pa.α) ) * ( pa.β*z^(1.0+pa.σ)/(1.0+pa.σ) - A );
    #Define function l(z,n)
    denominator_of_l(z,n) =  ss.λ*pa.χ*h_w - (ss.μ + z*ss.λ*h_e / (pa.σ* n^pa.α) )*(1+pa.ψ)*pa.χ/θ;
    l_fun(denominator) = denominator > 0.0 ? (ss.ω*θ*h_w/denominator)^(1.0/pa.ψ) : Inf;
    #Define function l(z,n,l)
    p(l,z,n) = pa.χ*l^(1.0+pa.ψ) / (θ*n^pa.α*(1.0-pa.β*z^pa.σ));

    #Itialize values for objective and potencial maximizers
    current_max = -Inf
    nn = NaN;
    zz = NaN;
    ll = NaN;
    pp = NaN;

    #Define a variable counting the number of gridpoints where no solution if found
    nosol = 0;

    #Check corner solution at z=0
    (potential_n,potential_z,potential_l,potential_p) = z_zero(θ, A, ss, pa);
    potential_value= objective(potential_z,potential_n,potential_l,potential_p);

    if isnan(potential_value)
        nosol = nosol + 1;
    elseif potential_value > current_max
        current_max = potential_value;
        nn= potential_n;
        zz= potential_z;
        ll= potential_l;
        pp= potential_p
    end

    #loop over grid of z. Recover n, l and p, and compute hamiltonian value
    for j=1:Nz
        #println(j)
        potential_z= zgrid[j];
        potential_n= n_fun(potential_z);

        denom = denominator_of_l(potential_z,potential_n)
        if denom <= 0.0
            potential_value = NaN;
        else
            potential_l = l_fun(denom);
            potential_p= p(potential_l, potential_z, potential_n);
            potential_value= objective(potential_z,potential_n,potential_l,potential_p);
        end

        if isnan(potential_value)
            nosol = nosol + 1;
        elseif potential_value > current_max
            current_max = potential_value;
            nn= potential_n;
            zz= potential_z;
            ll= potential_l;
            pp= potential_p
        end
        #println("z = ", potential_z, " n = ", potential_n, " l = ", potential_l, " p = ", potential_p, " value = ", potential_value)
    end

    nosol == Nz + 1 && error("no solution to hamiltonian for any z, adjust mu downward")

    zz, nn, ll, pp
end


function recover_controls!(ctrlvec::Array{Float64}, θvec::Array{Float64}, solvec::Array{Float64})
    (Nspan,~)=size(solvec)

    for j=1:Nspan
      θ = θvec[j];
      uw    = solvec[j,1];
      μ     = solvec[j,2];
      e     = solvec[j,3];
      ϕ_e   = solvec[j,4];
      λ     = solvec[j,6];
      ω     = solvec[j,8];
      ss = State(e, uw, ϕ_e, μ, λ, ω);

      (zz, nn, ll, pp) = new_find_controls( θ, ss, pa)
      ctrlvec[j,1] = zz;
      ctrlvec[j,2] = nn;
      ctrlvec[j,3] = ll;
      ctrlvec[j,4] = pp;
    end
    nothing
end

function z_zero(θ, A, ss, pa)
    A > 0.0 && ("No corner solution when A >0. Decrease A.")

    h_e= pa.he(θ, ss.e);
    h_w= pa.hw(θ, ss.e);
    ww = ss.ω/ss.λ

    nn =  -pa.α/(1.0-pa.α)*A/ww;
    kkappa = ss.λ*h_e*( ss.e + (A/(1.0-pa.α)) * (-(1.0-pa.α)/pa.α*ww/A)^pa.α )

    kkappa > 0.0 && ("Complementary slackness does not hold for z=0. Increase wage.")

    denominator_l = ss.λ*pa.χ*h_w - (ss.μ + kkappa )*(1+pa.ψ)*pa.χ/θ

    denominator_l < 0.0 && ("Planner wants infinite labor. Decrease mu")

    ll = (ss.ω*θ*h_w/denominator_l)^(1.0/pa.ψ)

    pp= pa.χ*ll^(1.0+pa.ψ) / (θ*nn^pa.α);

    (nn,0.0,ll,pp)
end



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
