
function find_states!(du,u,pa,θ,ini)
    μ     = u[1];
    e     = u[2];

    #Construct state object
    #ss = State(μ,e);
    ss = State(ini[1], ini[2]);

    #Find optimal controls
    (p) = new_find_controls( ini[3], ss, pa);

    h_e   = pa.he( θ, e);
    h_w   = pa.hw( θ, e);
    h_tot = h_w + p*h_e;

    λ = 1.0;
    du[1] = λ*h_tot;
    du[2] = p;

    nothing
end
