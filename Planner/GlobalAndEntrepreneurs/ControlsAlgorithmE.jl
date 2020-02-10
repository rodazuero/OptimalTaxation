function new_find_controlse( θ, θe, sse, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e= pa.he(θ, θe);

    #Definitions:
    n_max = ((sse.λe*pa.α^2.0*θe)/sse.ωe)^(1.0/(1.0-pa.α));
    z_max = (1.0/((1.0+pa.σ)*pa.β))^(1.0/pa.σ);
    n_lwbar = n_max;
    n_upbar = ((sse.λe*pa.α*θe)/sse.ωe)^(1.0/(1.0-pa.α));
    f_n(n) = pa.α*θe*n^pa.α - sse.ωe/sse.λe*n;
    g_z(z) = pa.α/pa.σ*z*(1.0-pa.β*z^pa.σ);
    global z_lwbar=0.0000016;

    #Hamiltonian:
    objective(ze,ne,κe) = pa.indicator*sse.ue^pa.ϕ*h_e + sse.λe*( θe*ne^pa.α - pa.β*ze^(1.0+pa.σ)/(1.0+pa.σ) - sse.ue)*h_e - sse.ωe*ne*h_e + sse.μe*ne^pa.α*(1.0-pa.β*ze^pa.σ) + κe*ze;

    #1. Intererior solution:
    #Defining the upper bar for z for the grid and the grid:
    if g_z(z_max) < f_n(n_max)

        global z_upbar = z_max;

    else

        fun_n_max = f_n(n_max);
        fun_z(z)  = pa.α/pa.σ*z*(1-pa.β*z^pa.σ) - fun_n_max;
        z_use     = find_zero(fun_z, (0,z_max), Bisection());

        z_upbar = z_use;

    end

    #println(z_upbar)

    Nz = 1000; #Number of Zs
    zstep = (z_upbar - z_lwbar)/(Nz-1);
    #zgrid = z_lwbar:zstep:z_upbar;
    zgrid = range(z_lwbar, stop = z_upbar, length = Nz);
    #println("z_grid = ", size(zgrid), "z_lwbar = ", z_lwbar, "z_upbar = ", z_upbar)
    #println("z_last = ", collect(zgrid)[end])

    #Define a variable counting the number of gridpoints where no solution if found
    global nosol = 0;
    global current_max  = NaN;
    global current_res  = NaN;
    global nne = NaN;
    global zze = NaN;
    global nne_res = NaN;
    global zze_res = NaN;
    #current_max  = NaN;
    #current_res  = NaN;
    #nne = NaN;
    #zze = NaN;
    #nne_res = NaN;
    #zze_res = NaN;

    #Finding the solutions for the equations:
    for j = 1:Nz

        #Solution for f_n:
        #println("j = ", j, "ze_grid = ", zgrid[j])
        potential_ze = zgrid[j];
        fun_z_grid   = g_z(potential_ze);

        fun_n(n) = pa.α*θe*n^pa.α - sse.ωe/sse.λe*n - fun_z_grid;

        if fun_n(n_lwbar)*fun_n(n_upbar)>0
            potential_ne = NaN;
            distance_f_n_g_ze = fun_n(n_lwbar);
            #println("distance_f_n_g_ze = ", distance_f_n_g_ze)
        else
            potential_ne = find_zero(fun_n, (n_lwbar,n_upbar), Bisection());
        end

        potential_ne <0 && error("N Entrepreneurs is negative.")

        potentialresidual = -sse.λe*pa.β*potential_ze^pa.σ*h_e-sse.μe*potential_ne^pa.α*pa.β*potential_ze^(pa.σ-1.0)*pa.σ;

        potential_value =  objective(potential_ze,potential_ne,0);

        #Evaluating the value of the Hamiltonian:
        if isnan(potential_value)

            global nosol = nosol + 1;

        elseif isnan(current_max)
            global current_max = potential_value;
            global nne         = potential_ne;
            global zze         = potential_ze;
            #current_max = potential_value;
            #nne         = potential_ne;
            #zze         = potential_ze;

        elseif potential_value > current_max
            global current_max = potential_value;
            global nne         = potential_ne;
            global zze         = potential_ze;
            #current_max = potential_value;
            #nne         = potential_ne;
            #zze         = potential_ze;

        end

        #Evaluating the value of the Residual:
        if isnan(current_res)
            global current_res = potentialresidual;
            global nne_res     = potential_ne;
            global zze_res     = potential_ze;
            #current_res = potentialresidual;
            #nne_res     = potential_ne;
            #zze_res     = potential_ze;

        elseif abs(potentialresidual) < abs(current_res)
            global current_res = potentialresidual;
            global nne_res     = potential_ne;
            global zze_res     = potential_ze;
            #current_res = potentialresidual;
            #nne_res     = potential_ne;
            #zze_res     = potential_ze;

        end

        #println("ze = ", potential_ze, " ne = ", potential_ne, " value = ", potential_value)
    end

    #Check solutions for N and Z are the same in Hamiltonian and residual:
    #if zze_res != zze
        #println("zze = ", zze, "zze_res = ", zze_res)
        #println("nne = ", nne, "nne_res = ", nne_res)
    #end
    #zze_res != zze && error("Solutions for z do not match.")
    #nne_res != nne && error("Solutions for n do not match.")

    #2. Corner Solution (if σ>1):
    if pa.σ>=1 && (sse.λe*θe*h_e*pa.α+sse.μe*pa.α)>0

        nn_num     = sse.λe*pa.α*θe*h_e + pa.α*sse.μe;
        corner_nne = (nn_num/(sse.ωe*h_e))^(1.0/(1.0-pa.α));

        corner_hamiltonian = objective(0.0,corner_nne,1.0);

        if corner_hamiltonian > current_max
            global current_max = corner_hamiltonian;
            global nne         = corner_nne;
            global zze         = 0.0;
        end

    else

        global nosol = nosol + 1;

    end

    #Finding if there is no solution at all:
    nosol == Nz + 1 && error("No solution to Hamiltonian for any z.")

    #Final Output:
    #println("zze = ", zze, "nne = ", nne)
    zze, nne

end

function recover_controlse!(ctrlvec::Array{Float64}, θw::Float64 ,θvec::Array{Float64}, solvec::Array{Float64})
    (Nspan,~)=size(solvec)

    for j=1:Nspan
      θe = θvec[j];
      ue    = solvec[j,1];
      μe     = solvec[j,2];
      λe     = solvec[j,4];
      ωe     = solvec[j,6];
      sse = StateE(ue, μe, λe, ωe);

      (zze, nne) = new_find_controlse( θw, θe, sse, pa)
      ctrlvec[j,1] = zze;
      ctrlvec[j,2] = nne;
    end
    nothing
end
