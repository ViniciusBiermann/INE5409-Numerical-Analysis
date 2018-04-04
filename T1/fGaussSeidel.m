function [x, flt_op] = fGaussSeidel(n, A, b, tolerancia)

    x = zeros(n,1);
    
    flt_op = 0;
    
    crit = 1;
    k = 0;
    lambda = 1.1;
    lambda_aux = 1 - lambda;
    
    while (crit > tolerancia && k < 100)
        aux = x;
        k = k + 1;
    
        x(1) = (lambda_aux)*aux(1) + lambda*((b(1) - A(1,2)*x(2))/A(1,1));
        flt_op = flt_op + 6;
    
        for i = 2 : n/2
            x(i) = (lambda_aux)*aux(i) + lambda*((b(i) - A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1))/A(i,i));
            flt_op = flt_op + 8;
        end
    
        for i = (n/2) + 1 : n - 1
            x(i) = (lambda_aux)*aux(i) + lambda*((b(i) - A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1))/A(i,i));
            flt_op = flt_op + 8;
        end
        
        x(n) = (lambda_aux)*aux(n) + lambda*((b(n) - A(n,n-1)*x(n-1))/A(n,n));
        flt_op = flt_op + 6;
        
        
        crit = max(abs(x - aux));
    end
    
    %k
    %crit
    %lambda
    
end
