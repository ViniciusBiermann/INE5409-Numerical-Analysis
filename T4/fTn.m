function Tn = fTn(n, a, b)
    h = (b-a)/n;
    x = a : h : b;
    y = 2/(sqrt(pi))*exp(-x.^2);
    
    Tn = h*0.5 * (y(1) + 2*sum(y(2:n)) + y(n+1));
    
end