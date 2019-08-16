function y = fPn(x, n, a)
    % Valor na base canônica
    for k = 1 : length(x)
        y(k) = a(1);
        xa = 1;
        for i = 2 : n + 1
            % Não é eficiente por causa da exponenciação
            % y = y + a(i) * x(k)^(i-1)
        
            xa *= x(k);
            y(k) += a(i)*xa;
        end
    end
    
end