function yp = fPnGregoryNewton(n, x, y, xp)
    % Calcular coeficientes
    % Diferença dividida ascendente com linha i e coluna k
    
    % para k = 1
    for i = 1 : n
        dda(i,1) = (y(i+1) - y(i)) / (x(i+1) - x(i));
    end
    
    % para k > 1
    for k = 2 : n
        for i = 1 : n + 1 - k
            dda(i,k) = (dda(i+1, k-1) - dda(i, k-1)) / (x(i+k) - x(i));
        end
    end
    % dda
    
    for k = 1 : length(xp)
        yp(k) = fPnValorNumGregoryNewton(xp(k), x, y, n, dda);
    end
    
end