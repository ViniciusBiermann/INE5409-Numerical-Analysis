function rmax = fresiduoTRD(n, t, r, d, b, x)

    r(1) = abs(r(1)*x(1) + d(1)*x(2) - b(1));
    r(n) = abs(t(n)*x(n-1) + r(n)*x(n) - b(n));
    
    for i = 2 : n-1
        r(i) = abs(t(i)*x(i-1) + r(i)*x(i) + d(i)*x(i+1) - b(i));
    end
    
    rmax = max(r);
end