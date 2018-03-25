function x = fGaussTRD(n, t, r, d, b)
    for i = 2 : n
        Aux = t(i)/r(i-1);
        t(i) = 0;
        r(i) = r(i) - Aux*d(i-1);
        b(i) = b(i) - Aux*b(i-1);
    end
    
    x(n) = b(n)/r(n);
    for i = n-1 : -1 : 1
        x(i) = (b(i) - d(i)*x(i+1))/r(i);
    end
end
