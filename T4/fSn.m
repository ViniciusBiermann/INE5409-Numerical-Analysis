function Sn = fSn(n, a, b)
    
    if(mod(n,2) != 0)
        n = n + 1
    end
    
    h = (b-a) / n;
    
    x = a : h : b;
    y  = 2/(sqrt(pi))*exp(-x.^2);
    
    Sn = (h/3) * (y(1) + 4*sum(y(2 : 2 : n)) + 2*sum(y(3 : 2 : n-1)) + y(n+1));
    
end