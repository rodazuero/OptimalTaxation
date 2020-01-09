function new_find_controls( θ, ss, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e = pa.he(θ, ss.e);
    h_w = pa.hw(θ, ss.e);

    #1. Interior solution:
    #Definitions:
    n_max   = ((ss.λ*pa.α^2.0*ss.e)/ss.ω)^(1.0/(1.0-pa.α));
    z_max   = (1.0/((1.0+pa.σ)*pa.β))^(1.0/pa.σ);
    n_lwbar = n_max;
    n_upbar = ((ss.λ*pa.α*ss.e)/ss.ω)^(1.0/(1.0-pa.α));
    f_n(n)  = pa.α*ss.e*n^pa.α - ss.ω/ss.λ*n;
    g_z(z)  = pa.α/pa.σ*z*(1.0-pa.β*z^pa.σ);
    z_lwbar = 0.0000016; #Random Number --> Not zero, because 0 is the corner solution.

    #Hamiltonian:
    objective(z,n,l,p,κ) = pa.indicator*ss.uw^pa.ϕ*(h_w+p*h_e) +ss.μ*pa.χ/θ*l^(1.0+pa.ψ) +
                           ss.λ*(ss.e*n^pa.α*p*h_e-pa.β/(1.0+pa.σ)*z^(1.0+pa.σ)*p*h_e-ss.uw*(h_w+p*h_e)-pa.χ/(1.0+pa.ψ)*l^(1.0+pa.ψ)*h_w)+
                           ss.ω*(θ*l*h_w-n*p*h_e)+ss.ϕ_e*p+κ*z;

    #Defining the upper bar for z for the grid and the grid:
    if g_z(z_max) < f_n(n_max)

        z_upbar = z_max;

    else

        fun_n_max = f_n(n_max);
        fun_z(z)  = pa.α/pa.σ*z*(1.0-pa.β*z^pa.σ) - fun_n_max;
        z_use     = find_zero(fun_z, (0.0,z_max), Bisection());

        z_upbar = z_use;

    end

    Nz = 500; #Number of Zs
    #zstep = (z_upbar-z_lwbar)/(Nz-1);
    #zgrid = z_lwbar:zstep:z_upbar;
    zgrid = range(z_lwbar, stop = z_upbar, length = Nz);
    #println("z_grid = ", size(zgrid), "z_lwbar = ", z_lwbar, "z_upbar = ", z_upbar)

    #Defining globals for loops --> If not defined the loops do not run:
    global nosol = 0; #Variable counting the number of gridpoints where no solution is found
    global current_max  = NaN;
    global current_res  = NaN;
    global nn = NaN;
    global zz = NaN;
    global ll = NaN;
    global pp = NaN;
    global nn_res = NaN;
    global zz_res = NaN;
    global ll_res = NaN;
    global pp_res = NaN;

    #Finding the solutions for the equations:
    for j = 1:Nz

        #Solution for f_n:
        #println("j = ", j)
        #println("j = ", j, "z_grid = ", zgrid[j])
        potential_z = zgrid[j];
        fun_z_grid  = g_z(potential_z);

        fun_n(n) = pa.α*ss.e*n^pa.α - ss.ω/ss.λ*n - fun_z_grid;

        if fun_n(n_lwbar)*fun_n(n_upbar)>0
            potential_n=NaN;
        else
            potential_n = find_zero(fun_n, (n_lwbar,n_upbar), Bisection());
        end

        potential_n <0.0 && error("N Global is negative.")

        potentialresidual = ss.λ*(ss.e*potential_n^pa.α-ss.ω/ss.λ*potential_n-potential_z/pa.σ*(1.0-pa.β/(1.0+pa.σ)*potential_z^pa.σ))*h_e +
                            (pa.indicator*ss.uw^pa.ϕ-ss.λ*ss.uw)*h_e + ss.ϕ_e;

        #Solution for l and p:
        den_l = (ss.λ*pa.χ*h_w-pa.χ/θ*(1.0+pa.ψ)*(ss.μ+ss.λ/pa.σ*potential_z/(potential_n^pa.α)*h_e));

        if den_l <= 0.0

            potential_l     = NaN;
            potential_p     = NaN;
            potential_value = NaN;

        else

            potential_l = ((ss.ω*θ*h_w)/den_l)^(1.0/pa.ψ);
            den_p       = θ*potential_n^pa.α*(1.0-pa.β*potential_z^pa.σ);
            potential_p = (pa.χ*potential_l^(1.0-pa.ψ))/den_p;
            potential_value = objective(potential_z,potential_n,potential_l,potential_p,0.0);
        end

        #Evaluating the value of the Hamiltonian:
        if isnan(potential_value)

            global nosol = nosol + 1;

        elseif isnan(current_max)
            global current_max = potential_value;
            global nn          = potential_n;
            global zz          = potential_z;
            global ll          = potential_l;
            global pp          = potential_p;

        elseif potential_value > current_max
            global current_max = potential_value;
            global nn          = potential_n;
            global zz          = potential_z;
            global ll          = potential_l;
            global pp          = potential_p;
        end

        #Evaluating the value of the Residual:
        if isnan(current_res)
            global current_res = potentialresidual;
            global nn_res     = potential_n;
            global zz_res     = potential_z;
            global ll_res     = potential_l;
            global pp_res     = potential_p;

        elseif abs(potentialresidual) < abs(current_res)
            global current_res = potentialresidual;
            global nn_res     = potential_n;
            global zz_res     = potential_z;
            global ll_res     = potential_l;
            global pp_res     = potential_p;
        end

        #To print solution eliminate # in the following:
        #println("z = ", potential_z, " n = ", potential_n, " l = ", potential_l, " p = ", potential_p,
        #         " value = ", potential_value)
    end

    #Check solutions for Z, N, L and P are the same in Hamiltonian and residual:
    if zz_res != zz
        println("zz = ", zz, "zz_res = ", zz_res)
        #println("nn = ", nn, "nn_res = ", nn_res)
        #println("ll = ", ll, "ll_res = ", ll_res)
        #println("pp = ", pp, "pp_res = ", pp_res)
    end

    #zz_res != zz && error("Solutions for z do not match.")
    #nn_res != nn && error("Solutions for n do not match.")
    #ll_res != ll && error("Solutions for l do not match.")
    #pp_res != pp && error("Solutions for p do not match.")

    #2. Corner Solution (if σ<1):
    corner_nn = ((ss.λ*pa.α*ss.e)/ss.ω)^(1/(1-pa.α));
    condition = (pa.indicator*ss.uw^pa.ϕ-ss.λ*ss.uw)*h_e-ss.λ*ss.e*corner_nn^pa.α*h_e-ss.ω*corner_nn*h_e+ss.ϕ_e;

    if pa.σ<1.0 && condition==0.0

        cor_l_den = ss.λ*pa.χ*h_w-ss.μ*(1.0+pa.ψ)*pa.χ/θ;
        corner_ll = ((ss.ω*θ*h_w)/cor_l_den)^(1.0/pa.ψ);
        corner_pp = (pa.χ*corner_ll^(1+pa.ψ))/(θ*corner_nn^pa.α);

        corner_hamiltonian = objective(0.0,corner_nn,corner_ll,corner_pp,1.0);

        if corner_hamiltonian > current_max
            global current_max = corner_hamiltonian;
            global nn          = corner_nn;
            global zz          = 0.0;
            global ll          = corner_ll;
            global pp          = corner_pp;
        end

    else

        global nosol = nosol + 1;

    end

    #Finding if there is no solution at all:
    nosol == Nz + 1 && error("No solution to Hamiltonian for any z.")

    #Final Output:
    #println("zz = ", zz, "nn = ", nn, "ll = ", ll, "pp = ", pp)
    zz, nn, ll, pp

end

function recover_controls!(ctrlvec::Array{Float64}, θvec::Array{Float64}, solvec::Array{Float64})
    (Nspan,~)=size(solvec)

    for j=1:Nspan
      θ     = θvec[j];
      uw    = solvec[j,1];
      μ     = solvec[j,2];
      e     = solvec[j,3];
      ϕ_e   = solvec[j,4];
      λ     = solvec[j,6];
      ω     = solvec[j,8];
      ss    = State(e, uw, ϕ_e, μ, λ, ω);

      (zz, nn, ll, pp) = new_find_controls( θ, ss, pa)
      ctrlvec[j,1] = zz;
      ctrlvec[j,2] = nn;
      ctrlvec[j,3] = ll;
      ctrlvec[j,4] = pp;
    end
    nothing
end
