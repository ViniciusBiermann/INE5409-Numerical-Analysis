function M = fMultiplicidade(r)
        % contador de |restos| =< limite
        numResto = length(r);
        
        sumResto = abs(r(1) + r(2));
        
        M = 1;
        limite = 1e-1;
        
        while(sumResto < limite)
            M = M + 1;
            sumResto += abs(r(M+1));
        end
end