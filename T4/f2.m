function y = f2(a, t, v)
    y(1) = sum( (log( a(1) .+ a(2) .* t.*t ) .- v ) .* 1./(a(1) .+ a(2) .* t.*t ));
    y(2) = sum( (log( a(1) .+ a(2) .* t.*t ) .- v ) .* (t.*t)./(a(1) .+ a(2) .* t.*t ));
end