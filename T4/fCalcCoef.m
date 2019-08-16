function a = fCalcCoef(m,T,V)
    
    A(1,1) = sum(T.*T);
    A(1,2) = sum(T .* cos(T));
    A(2,1) = A(1,2);
    A(2,2) = sum(cos(T).*cos(T));
    
    B(1) = sum(T .* V); 
    B(2) = sum(V .* cos(T)); 

    a = fCholesky(2, A, B);
end