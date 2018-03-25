function A = fpivotpar(A, k, n)
    maxVal = abs(A(k,k));
    maxLine = k;
    
    for i = k+1 : n
        Aux = abs(A(i,k));
        if Aux > maxVal
            maxVal = Aux;
            maxLine = i;
        end
    end
    
    for j = 1 : n+1
        Aux = A(k,j);
        A(k,j) = A(maxLine, j);
        A(maxLine, j) = Aux;
    end
end
