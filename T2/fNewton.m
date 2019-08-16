function x = fNewton(xi)
    
    for i = 1 : 10
        deltaX = -f(xi)./fdf(xi);
        x = xi .+ deltaX;
        xi = x;
    end
    
end