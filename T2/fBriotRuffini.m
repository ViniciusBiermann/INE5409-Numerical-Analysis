function [n, b, r] = fBriotRuffini(n, a, xi)
        b(1) = a(1);

        for i = 2 : n + 1
                b(i) = a(i) + b(i-1)*xi;
        end
        
        r = b(n+1);
        b = b(1 : n);
        n = n - 1;
end