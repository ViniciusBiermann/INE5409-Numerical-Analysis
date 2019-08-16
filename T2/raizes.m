function x = raizes(a, m)
        n = length(a) - 1;
        
        a = a(2 : n + 1)./a(1);
        a = [1 a];
        
        i = 1; % a primeira raiz
        
        while(n > 0)
            xi = fLocalizacaoPn(n, a);
            
            [x(i), M(i)] = fNewtonPn(xi, n, a, m);
        
            % Reduzindo o grau de n para n-1
            for j = 1 : M(i)
                [n, a, r] = fBriotRuffini(n, a, x(i));
            end
            
            x(i+1 : i+M(i)-1) = x(i);
            
            %n
            i = i + M(i);
        end
        
        x = transpose(x);
end