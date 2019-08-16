function [x, M] = fNewtonPn(xi, n, a, m)

        dx = 1;
        k = 0;
        tolerancia = 1e-14;
        nOriginal = n;
        aOriginal = a;

        while(abs(dx) > tolerancia && k < 100)
            k = k + 1;
            % n, a são grau e coeficientes do polinomio quociente resultante da divisão
            % r(1) é o resto
            % r(2) é da derivada
            
            n = nOriginal;
            a = aOriginal;
            
            for i = 1 : nOriginal
                [n, a, r(i)] = fBriotRuffini(n, a, xi);
            end

            r(nOriginal + 1) = 1;
            
            %[r(1) r(2) r(3) r(4)]
            
            if(m == 1)
                M = 1;
                dx = - r(1) / r(2);
            else
                % Estimar restos menores que um limite (estimar a multiplicidade)
                M = fMultiplicidade(r);
                dx = - r(M) / (M * r(M + 1));
            end
            
            x = xi + dx;
            xi = x;
        end
end
