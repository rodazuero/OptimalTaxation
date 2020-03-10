function new_find_controlse( θ, θe, sse, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e= pa.he(θ, θe);

    #Defining bounds we use in various cases (limits of z):
    n_full_info = ((sse.λe*pa.α*θe)/sse.ωe)^(1.0/(1.0-pa.α));
    z_max   = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    z_lwbar = eps(); #The smallest possible z.
    z_1     = z_max;
    z_2     = -pa.σ*sse.μe*n_full_info^pa.α/(sse.λe*h_e);
    #Keep the lower number:
    z_upbar = min(z_1,z_2);

    Nz    = 1000; #Number of Zs
    zstep = (z_upbar - z_lwbar)/(Nz-1);
    zgrid = range(z_lwbar, stop = z_upbar, length = Nz);

    #Defining the functions we are using:
    n_opt(z)   = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
    fun_z(z,n) = (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ));
    fun_nz(z)  = fun_z(z,n_opt(z));

    #Hamiltonian:
    objective(ze,ne) = pa.indicator*sse.ue^pa.ϕ*h_e + sse.λe*( θe*ne^pa.α - pa.β*ze^(1.0+pa.σ)/(1.0+pa.σ) - sse.ue)*h_e - sse.ωe*ne*h_e + sse.μe*ne^pa.α*(1.0-pa.β*ze^pa.σ);

    #Defining globals for loop:
    global nosol = 0;
    global current_max  = NaN;
    global nne = NaN;
    global zze = NaN;

    #1. The interior solution is:
    #Finding the solutions for the equations:
    for j = 1:Nz

        #Solution for f_n:
        #println("j = ", j, "ze_grid = ", zgrid[j])
        potential_ze = zgrid[j];

        potential_ne    = n_opt(potential_ze);;
        potential_value = objective(potential_ze,potential_ne);

        #println("z =   ", potential_ze, "  n =", potential_ze, "  value =", potential_value )

        #Evaluating the value of the Hamiltonian:
        if isnan(potential_value)

            global nosol = nosol + 1;

        elseif isnan(current_max)
            global current_max = potential_value;
            global nne         = potential_ne;
            global zze         = potential_ze;

        elseif potential_value > current_max
            global current_max = potential_value;
            global nne         = potential_ne;
            global zze         = potential_ze;

        end
    end

    #2. The corner solution is:
    corner_zz = 0.0;
    corner_nn = 0.0;

    corner_hamiltonian = objective(corner_zz,corner_nn);

    if corner_hamiltonian > current_max
        global current_max = corner_hamiltonian;
        global zze = corner_zz;
        global nne = corner_nn;
    end

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
