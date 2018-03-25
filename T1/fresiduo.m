function rmax = fresiduo(A, b, x)
    n = size(A,1);
    
    for i = 1 : n
        SumAux = 0;
        for j = 1 : n
            SumAux = SumAux + A(i,j)*x(j);
        end
        r(i) = abs(SumAux - b(i));
    end
    % imprimir lista de residuos
    %r
    rmax = max(r);
end
