function x = fLocalizacaoPn(n, a)
        r = 1 + max(abs(a(2 : n + 1)));
        x = 0.5*r*complex(cos(pi/4), sin(pi/4));
end