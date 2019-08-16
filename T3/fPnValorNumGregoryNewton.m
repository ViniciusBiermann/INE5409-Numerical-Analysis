function yp = fPnValorNumGregoryNewton(xp, x, y, n, dda)
    yp = y(1);
    
    soma = 0;
    for k = 1 : n
    
        produtorio = 1;
        for i = 1 : k
            produtorio *= (xp - x(i));
        end
        
        soma += dda(1, k)*produtorio;
    end
    yp += soma;
end