function c = fCoefMaclaurin(grau)
    c(1:grau+1) = 0;
    for k = 0 : grau/2
        c(2*k + 1) = (-4)^k/factorial(k);
    end
end