function x = fGauss(A,b)
    n = size(A,1);
    A = [A b;];
    
    % -- Fase de escalonamento --
    for k = 1 : n - 1
        A = fpivotpar(A, k, n);
        for i = k + 1 : n
            aux = A(i,k)/A(k,k);
            for j = k + 1 : n + 1
                A(i,j) = A(i,j) - aux*A(k,j);
            end
            A(i,k) = 0;
        end
    end
    % A
    
    % -- Fase de Retrosubstituicao --
    
    if abs(A(n,n)) <= 1e-14
        if abs(A(n,n+1)) <= 1e-14
            x(n) = 0;   % Valor qualquer escolhido
            fprintf("\n - Sistema possível e indeterminado.\n");
        else
            x(n) = NaN;     % Not a number
            fprintf("\n - Sistema impossível.\n");
        end
    else
        x(n)=A(n,n+1)/A(n,n);
        fprintf("\n - Sistema possível e determinado.\n");
    end
    
    x(n) = A(n,n+1)/A(n,n);
    for i = n - 1 : -1 : 1
        SumAux = 0;
        for j = i + 1 : n
            SumAux = SumAux + (A(i,j)*x(j));
        end
        x(i) = (A(i,n+1) - SumAux)/A(i,i);
    end
end
