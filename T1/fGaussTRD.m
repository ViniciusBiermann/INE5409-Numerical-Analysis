function [x, flt_op] = fGaussTRD(n, t, r, d, b)
    flt_op = 0;
    for i = 2 : n
        Aux = t(i)/r(i-1);
        t(i) = 0;
        r(i) = r(i) - Aux*d(i-1);
        b(i) = b(i) - Aux*b(i-1);
        flt_op = flt_op + 5;
    end
    
    x(n) = b(n)/r(n);
    flt_op = flt_op + 1;
    for i = n-1 : -1 : 1
        x(i) = (b(i) - d(i)*x(i+1))/r(i);
        flt_op = flt_op + 3;
    end
    
end
