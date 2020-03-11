function new_find_controls( θ, ss, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions:
    h_e = pa.he(θ, ss.e);
    h_w = pa.hw(θ, ss.e);

    #Solutions depending on the bound of A:
    #We have 2 major cases:
    n_full_info = ((ss.λ*pa.α*ss.e)/ss.ω)^(1.0/(1.0-pa.α));
    z_max  = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    A_cons = (pa.indicator*ss.uw^pa.ϕ-ss.λ*ss.uw)+ss.ϕ_e/h_e;
    cond   = -(1.0-pa.α)/pa.α*ss.ω*n_full_info
    #println("A_cons = ", A_cons, "  n_full = ", n_full_info, "  condition", cond )

    #Defining the functions we are using:
    n_opt(z)   = 1.0/ss.ω*pa.α/(1.0-pa.α)*(ss.λ*pa.β/(1.0+pa.σ)*z^(1.0+pa.σ)-A_cons);
    fun_z(z,n) = (ss.λ*ss.e*n^pa.α-ss.ω*n-ss.λ/pa.σ*z*(1.0-pa.β/(1.0+pa.σ)*z^pa.σ)+A_cons);
    fun_nz(z)  = fun_z(z,n_opt(z));

    den_l(n,z)  = (ss.λ*pa.χ*h_w-pa.χ/θ*(1.0+pa.ψ)*(ss.μ+h_e/(n^pa.α*(1.0-pa.β*z^pa.σ))*(ss.λ*ss.e*n^pa.α-
                   ss.λ*pa.β/(1.0+pa.σ)*z^(1.0+pa.σ)-ss.ω*n + A_cons)));
    ll_opt(n,z) = ((ss.ω*θ*h_w)/den_l(n,z))^(1.0/pa.ψ);
    den_p(n,z)  = θ*n^pa.α*(1.0-pa.β*z^pa.σ);

        #The Hamiltonian:
    objective(z,n,l,p) = pa.indicator*ss.uw^pa.ϕ*(h_w+p*h_e) +ss.μ*pa.χ/θ*l^(1.0+pa.ψ) +
                         ss.λ*(ss.e*n^pa.α*p*h_e-pa.β/(1.0+pa.σ)*z^(1.0+pa.σ)*p*h_e-ss.uw*(h_w+p*h_e)-pa.χ/(1.0+pa.ψ)*l^(1.0+pa.ψ)*h_w)+
                         ss.ω*(θ*l*h_w-n*p*h_e)+ss.ϕ_e*p; #Hamiltonian:

    nn_z0 = -A_cons*1.0/ss.ω*pa.α/(1.0-pa.α);
    #Defining bounds (limits of z) for cases 1 and 3:
    if nn_z0 >= n_full_info

        println("1 case: A_cons implies n>n_full_info, z==0.")
        zz = 0.0;
        nn = nn_z0;
        den_l(nn,zz) <= 0.0 && error("Denominator of L is negative or zero.") #Stop when denominator is negative.
        ll = ll_opt(nn,zz);
        pp = (pa.χ*ll^(1.0+pa.ψ))/den_p(nn,zz);
        κnz = 0.0;

    else
        z_lwbar = 0.0; #The smallest possible z.
        z_1     = z_max;
        z_2     = ((1.0+pa.σ)/(ss.λ*pa.β)*(A_cons+(1.0-pa.α)/pa.α*ss.ω*n_full_info))^(1.0/(1.0+pa.σ));
        #Keep the lower number:
        z_upbar = min(z_1,z_2);

        if A_cons <= 0

            println("1 case: A_cons<0. A_cons > OtherValue.")
            #1. Evaluating the interior solution:
            #Using bisection method:
            if fun_nz(z_lwbar)*fun_nz(z_upbar)<0
                potential_z = find_zero(fun_nz, (z_lwbar,z_upbar), Bisection());
            else
                error("Bisection: signs equal --> Cannot solve.")
            end

            potential_n = n_opt(potential_z);
            den_l(potential_n,potential_z) <= 0.0 && error("Denominator of L is negative or zero.") #Stop when denominator is negative.
            potential_l = ll_opt(potential_n,potential_z);
            potential_p = (pa.χ*potential_l^(1.0+pa.ψ))/den_p(potential_n,potential_z);

            interior_hamiltonian = objective(potential_z,potential_n,potential_l,potential_p);

            #2. Evaluating the corner solution:
            corner_zz = 0.0;
            corner_nn = nn_z0;
            den_l(corner_nn,corner_zz) <= 0.0 && error("Denominator of L (corner) is negative or zero.") #Stop when denominator is negative.
            corner_ll = ll_opt(corner_nn,corner_zz);
            corner_pp = (pa.χ*corner_ll^(1.0+pa.ψ))/den_p(corner_nn,corner_zz);

            corner_hamiltonian = objective(corner_zz,corner_nn,corner_ll,corner_pp);

            #3. Defining the solution we are keeping:
            if corner_hamiltonian > interior_hamiltonian
                #zz = corner_zz;
                #nn = corner_nn;
                #ll = corner_ll;
                #pp = corner_pp;
                error("corner_hamiltonian > interior_hamiltonian.")
            else
                zz = potential_z;
                nn = potential_n;
                ll = potential_l;
                pp = potential_p;
                κnz = 0.0;
            end

        else

        println("3 case: A_cons > 0.")

                n_lwbar = 1.0e-10;
                n_upbar = n_full_info;

                fun_zagrzero(n)     = ss.e*n^pa.α - ss.ω/ss.λ*n;
                lkappa_nzagrzero(n) = pa.σ*pa.β*fun_zagrzero(n)^(pa.σ-1.0)/(1.0-pa.β*fun_zagrzero(n)^pa.σ)*(A_cons +
                                      ss.λ*ss.e*n^pa.α - ss.ω*n - ss.λ*fun_zagrzero(n)/pa.σ*(1.0-pa.β/(1.0+pa.σ)*fun_zagrzero(n)^pa.σ));
                funct_nagrzero(n)   = pa.α/n*(-(1.0-pa.α)/pa.α*ss.ω*n + ss.λ/(1.0+pa.σ)*pa.β*fun_zagrzero(n)^(1.0+pa.σ) - A_cons) +
                                      lkappa_nzagrzero(n)*(pa.α*ss.e*n^(pa.α-1.0) - ss.ω/ss.λ);

                #Solving for n using bisection:
                if funct_nagrzero(n_lwbar)*funct_nagrzero(n_upbar)<0
                    potential_n= find_zero(funct_nagrzero, (n_lwbar,n_upbar), Bisection());
                else
                    error("Bisection: signs equal --> Cannot solve.")
                end

                potential_z = fun_zagrzero(potential_n);
                potential_n,potential_z <= 0.0 && error("Z is negative or zero.") #Stop when denominator is negative.
                den_l(potential_n,potential_z) <= 0.0 && error("Denominator of L is negative or zero.") #Stop when denominator is negative.
                potential_l = ll_opt(potential_n,potential_z);
                potential_p = pa.χ*potential_l^(1.0+pa.ψ)/den_p(potential_n,potential_z);
                κnz = lkappa_nzagrzero(potential_n)/ss.λ;

                zz = potential_z;
                nn = potential_n;
                ll = potential_l;
                pp = potential_p;

        end
    end

    #Final Output:
    #println("zz = ", zz, "nn = ", nn, "ll = ", ll, "pp = ", pp)
    zz, nn, ll, pp, κnz;
end

function recover_controls!(ctrlvec::Array{Float64}, θvec::Array{Float64}, solvec::Array{Float64})
    (Nspan,~)=size(solvec)

    for j=Nspan:-1:1
      θ     = θvec[j];
      uw    = solvec[j,1];
      μ     = solvec[j,2];
      e     = solvec[j,3];
      ϕ_e   = solvec[j,4];
      λ     = solvec[j,6];
      ω     = solvec[j,8];
      ss    = State(e, uw, ϕ_e, μ, λ, ω);

      #(zz, nn, ll, pp) = new_find_controls( θ, ss, pa)
      (zz, nn, ll, pp) = new_find_controls( θ, ss, pa)
      ctrlvec[j,1] = zz;
      ctrlvec[j,2] = nn;
      ctrlvec[j,3] = ll;
      ctrlvec[j,4] = pp;
    end
    nothing
end
