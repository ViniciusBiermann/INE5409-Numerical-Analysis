function x = fGaussSeidel(n, A, b)

    x = zeros(n,1);
    
    crit = 1;
    k = 0;
    lambda = 1.1;
    
    while (crit > 1e-6 && k < 100)
        aux = x;
        k = k + 1;
    
        x(1) = (1-lambda)*aux(1) + lambda*((b(1) - A(1,2)*x(2))/A(1,1));
    
        for i = 2 : n/2
            x(i) = (1-lambda)*aux(i) + lambda*((b(i) - A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1))/A(i,i));
        end
    
        for i = (n/2) + 1 : n - 1
            x(i) = (1-lambda)*aux(i) + lambda*((b(i) - A(i,i-1)*x(i-1) - A(i,i+1)*x(i+1))/A(i,i));
        end
        
        x(n) = (1-lambda)*aux(n) + lambda*((b(n) - A(n,n-1)*x(n-1))/A(n,n));
        
        
        crit = max(abs(x - aux));
    end
    
    k
    crit
    lambda
    
end
