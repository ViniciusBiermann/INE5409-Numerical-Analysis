function [x, flt_op] = fGauss(A,b)
    n = size(A,1);
    A = [A b;];
    
    flt_op = 0;
    % -- Fase de escalonamento --
    for k = 1 : n - 1
        A = fpivotpar(A, k, n);
        for i = k + 1 : n
            aux = A(i,k)/A(k,k);
            flt_op = flt_op + 1;
            for j = k + 1 : n + 1
                A(i,j) = A(i,j) - aux*A(k,j);
                flt_op = flt_op + 2;
            end
            A(i,k) = 0;
        end
    end
    % A
    
    % -- Fase de Retrosubstituicao --
    
    if abs(A(n,n)) <= 1e-14
        if abs(A(n,n+1)) <= 1e-14
            x(n) = 0;   % Valor qualquer escolhido
            %fprintf("\n - Sistema possível e indeterminado.\n");
        else
            x(n) = NaN;     % Not a number
            %fprintf("\n - Sistema impossível.\n");
        end
    else
        x(n)=A(n,n+1)/A(n,n);
        flt_op = flt_op + 1;
        %fprintf("\n - Sistema possível e determinado.\n");
    end
    
    %x(n) = A(n,n+1)/A(n,n);
    %flt_op = flt_op + 1;
    for i = n - 1 : -1 : 1
        SumAux = 0;
        for j = i + 1 : n
            SumAux = SumAux + (A(i,j)*x(j));
            flt_op = flt_op + 2;
        end
        x(i) = (A(i,n+1) - SumAux)/A(i,i);
        flt_op = flt_op + 2;
    end
    
end
