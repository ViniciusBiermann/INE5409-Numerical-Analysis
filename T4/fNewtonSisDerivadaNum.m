function x = fNewtonSisDerivadaNum(xi, z, w)
        
        k = 0;
        tolerancia = 1e-12;
        
        fun = f2(xi, z, w);
        n = length(fun);

        for l = 1 : n
            dx(l) = 1e-6;
        end
        
        while(max(abs(dx)) > tolerancia && k < 100)
            k = k + 1;
            
            xn = xi;
            
            for j = 1 : n
                xn(j) = xi(j) .+ dx(j);
                g = (f2(xn, z, w) .- fun)/dx(j);
                for i = 1 : n
                    A(i,j) = g(i);
                end
                xn = xi;
            end

            b = -(f2(xi, z, w));
            
            dx = A\transpose(b);
            
            x = xi + dx;
            xi = x;
        end
        
        valorK = k;
        criterioDeParada = max(abs(dx));
end	 	  	 	 			      	  		      	 	
