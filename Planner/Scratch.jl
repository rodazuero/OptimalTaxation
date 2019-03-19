Scratch
#Define nodes and weights for integration
nodes, weights = gausslegendre( 100 )

#Define the h functions
function hw_fun(θ::Real, e::Real , g::Function, node::Array, weights::Array, θe_low::Real)
    args= (nodes+1)(e - θe_low)/2 + θe_low;
    hw= (e - θe_low)/2 * dot( weights, g(θ, args) )
end

function he_fun(θ::Real, e::Real , g::Function, node::Array, weights::Array, θw_low::Real)
    args= (nodes+1)(θ - θw_low)/2 + θw_low;
    he= (θ - θ_low)/2 * dot( weights, g(args, e) )
end

function h_fun(p::Real, θ::Real, e::Real , g::Function, node::Array, weights::Array, θw_low::Real, θe_low::Real)
    h= hw_fun(θ, e, g, nodes, weights, θe_low) + p*he_fun(θ, e, g, nodes, weights, θw_low)
end




#Define function
DiffEq = @ode_def SystemOfDifferentialEquations begin
  de = p
  du = χ l^(1+ψ)/θ
  dμ = λ h_fun(e,g)
end χ ψ



function jacobian_find_z!(jac_fzero,z, θ, λ, ω, state, p)
    h_e= p.he( θ, state.e);
    h_w= p.hw(θ, state.e);

    numerator= p.β*z^(1+p.σ)/(1+p.σ)h_e-state.ϕ_e/(λ*h_e)-(1-p.α)*z*(1-p.β*z^p.σ)/p.σ;
    denominator = p.α*(1-p.α)*state.e;

    nα = numerator/denominator;

    term1 = (p.β*z^p.σ*(2+1/p.σ)-1)/denominator;
    term2 = 1/p.α *nα^(p.α/(1-p.α)) *term1;

    jac_fzero = λ*state.e*p.α*term1 - ω*term2 - λ*p.α*(1-(p.β*z^p.σ)*(1+p.σ))/p.σ;
end


function integrate_dg_de(θ,e,pa)
    x1(s) = log(s);
    x2 = log(e)
    ff(s) = 1/e*( -2*pa.σ_we*(x1(s)-pa.μ_w)/((pa.σ2_w*pa.σ2_e)^0.5) + 2*(x2 - pa.μ_e)/pa.σ2_e ) *
            (2.0*π*(pa.σ2_w*pa.σ2_e*(1-pa.σ_we^2.0))^0.5 )^(-1.0) * (-1/(2*(1-(1-pa.σ_we^2.0)))) *
            exp((2.0*(1-pa.σ_we^2.0))*(-(x2-pa.μ_e)^2.0/(pa.σ2_e))) *
            exp(pa.σ_we*(x1(s)-pa.μ_w)*(x2 - pa.μ_e)/((1-pa.σ_we^2.0)*(pa.σ2_w*pa.σ2_e)^0.5)) *
            (2.0*(1-pa.σ_we^2.0)*pa.σ_w)^0.5*s;

    u(s) = (x1(s)-pa.μ_w)/((2.0*(1-pa.σ_we^2.0)*pa.σ2_w)^0.5);

    a = (0.0 -pa.μ_w)/((2.0*(1-pa.σ_we^2.0)*pa.σ2_w)^0.5);
    b = u(θ);
    v(t) = (a+b)/2.0 + (b-a)/2.0 * tanh(t);


    (N,)=size(pa.nodes);
    integral=0;
    for i=1:N
        t= u(pa.nodes[i]);
        integral +=( (b-a)/2.0 * ff(v(t)) * sech(t)^2.0 )*pa.weights[i];
    end
    integral
end



function find_z( θ, λ, ω, ss, pa)
    h_e= pa.he( log(θ), log(ss.e));
    h_w= pa.hw( log(θ), log(ss.e));

    param1 =  ( pa.ϕ * ss.uw^(pa.ϕ - 1) - λ * ss.uw ) * h_e + ss.ϕ_e;
    param2 = pa.α * (1+pa.σ);

    N(z) = λ*z/pa.σ * (1-pa.α - (pa.β*z^pa.σ)*( 1/(1+pa.σ) - pa.α ) ) - param1/h_e; #cf eq 21
    nn(z) = (N(z)^(1/pa.α)) / ((1-pa.α)*λ*ss.e);
    L(z) = λ*ss.e*pa.α*nn(z)^pa.α - ω*nn(z); #cf. eq 22
    R(z) = λ*pa.α*z/pa.σ*(1-pa.β*z^pa.σ);
    eq22(z)= L(z) - R(z);

    if param1 < 0.0
        z_lb = 0.0;
        z_ub = (1/pa.β)^(1/pa.σ);
        if  param2 > 1.0  # Convex case
            z = find_zeros( eq22, z_lb, z_ub );
        else #concave case
            if N(z_ub) > 0
                z = find_zeros( eq22, z_lb, z_ub );
            else
                eror("line 72, to be filed later")
                #find max z such that n is positive
            end
        end
    else
        if  param2 > 1.0 #Convex case
            z_ub = (1/pa.β)^(1/pa.σ) ;
            if N(z_ub) > 0
                z_lb = find_zeros( N, 0.0, z_ub );
                z = find_zeros( eq22, z_lb[1], z_ub );
                #No solution
            else
                error("Case to be resolved. No possible solution")
            end
        else
            z_max = ( (1-pa.α) / (pa.β*(1-pa.α*(1+pa.σ))) )^(1/pa.σ);
            if N(z_max)>0
                z_lb = find_zeros( N, 0.0, z_max );
                z_ub = find_zeros( N, z_max, Inf );
                z = find_zeros( eq22, z_lb[1], z_ub[1] );
            else
                error("Case to be resolved. No possible solution")
            end
        end
    end
    z
end





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
