function GTm = fGTm(m, fT)
    soma = 0;
    
    for k = 1 : m
        soma = soma + fT(cos((2*k - 1)*pi/(2*m)));
    end
    
    GTm = pi / m*soma;
end