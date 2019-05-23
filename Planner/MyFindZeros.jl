function my_find_zeros_inf(f::Function, lb::Real, N::Int64 )
    #This funtion looks for roots in the interval (lb, infty)
    roots= [];
    for i=-5:N
        a = (10.0^(i-1) - 10.0^-6) +lb ;
        b = (10.0^(i) - 10.0^-6)+lb;
#      println("a= ", Float64(a), " b= ", Float64(b), " f(a)= ", Float64(f(a)), " f(b)= ",Float64(f(b)))
        if f(a)*f(b) == 0.0
            root = f(a) == 0.0 ? a : b;
            append!(roots,root);
        elseif f(a)*f(b) < 0.0
            root = find_zero(f, (a, b), Bisection(),  atol = 4.0 * eps(real(Float64)), rtol = 4.0 * eps(real(Float64)) );
            append!(roots,root);
        else
            try
                root = find_zero(f, (a + b)/2, Order0());
                append!(roots,root);
            catch

            end
        end
    end
    a = (10.0^(N) -1) + lb;
    b = Inf;
    # println("a= ", a, " b= ",b, " f(a)= ", f(a), " f(b)= ",f(b))

    if f(a)*f(b) <= 0.0
        global i = -1;
        global flag = true;
        while flag
            global i += 1
            a = (10.0^(N + i) -1) + lb;
            b = (10.0^(N + i + 1 ) -1) + lb;
            global flag = f(a)*f(b) > 0.0
            println("a= ", a, " b= ",b, " f(a)= ", f(a), " f(b)= ",f(b))
            println("i= ", i, " flag ", flag)
        end
        a = (10.0^(N + i ) -1) + lb;
        b = (10.0^(N + i +1) -1) + lb;
#        println("a= ", a, " b= ",b, " f(a)= ", f(a), " f(b)= ",f(b))
        root = find_zero(f, (a, b), Bisection());
        append!(roots,root);
    else
        try
            root = find_zero(f, (a + b)/2, Order0());
            append!(roots,root);
        catch

        end
    end
    println("roots = ", roots)
    roots
end


function my_find_zeros(f::Function, lb::Real, ub::Real, N::Int64 )
    #This funtion looks for roots in the interval (lb, ub)
    roots= [];
    for i=1:N
        a = (i-1)*(ub-lb)/N + lb;
        b = i*(ub-lb)/N + lb;
#         println("a= ", Float64(a), " b= ", Float64(b), " f(a)= ", Float64(f(a)), " f(b)= ",Float64(f(b)))
         if f(a)*f(b) == 0.0
             root = f(a) == 0.0 ? a : b;
             append!(roots,root);
         elseif f(a)*f(b) < 0.0
            root = find_zero(f, (a, b), Bisection());
            append!(roots,root);
        else
            try
                root = find_zero(f, (a + b)/2, Order0());
                append!(roots,root);
            catch

            end
        end
    end

    roots
end



function my_newton(f::Function, fprime::Function, x0::Float64)
    #This funtion looks for roots in the interval (lb, infty)
    x=BigFloat(x0);
    ff =BigFloat(f(x));
    ftol=10^-5
    while abs(ff) > ftol
        next_x= x - f(x)/fprime(x);

        x=next_x;
        println("x = ", x)
        ff=BigFloat(f(x))
        println("f(x) = ", ff)
    end
    x
end




function my_bisection(f::Function, lb::Real, ub::Real, tol:: Real )
    #This funtion looks for roots in the interval (lb, ub)
    f(lb)*f(ub) > 0.0 && error("No bracketing interval")

    midpoint = (lb+ub)/2.0

    while f(midpoint) > tol

        if f(lb)*f(midpoint) <=0.0
            ub = midpoint
        else
            lb = midpoint
        end
        midpoint = (lb+ub)/2.0
    end

    midpoint

end
