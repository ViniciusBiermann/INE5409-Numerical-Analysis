function xi = fLocalizacao(A, B)
    h = 0.01;
    a = A;
    b = a + h;
    i = 0;
    lim = 6;
    
    while a < B
        if(f(a)*f(b)) < 0 && abs(f(a)) < lim && abs(f(b)) < lim
            ++i;
            xi(i) = 0.5*(a+b);
        end
        a = b;
        b = a + h;
    end
    
end