

function [fi] = evol(fi, F, deltaT, deltaX, deltaY)

    m = size(fi,1);
    n = size(fi,2);
    sigmaFi = zeros(m,n);
    for i = 2:m-1
        for j= 2:n-1
           if F(i,j) > 0
               sigmaFi(i,j) = deltaT*F(i,j)*sqrt( ...
                   max((fi(i,j) - fi(i-1,j))/deltaX,0)^2 + ...
                   min((fi(i+1,j) - fi(i,j))/deltaX,0)^2 + ...
                   max((fi(i,j) - fi(i,j-1))/deltaY,0)^2 + ...
                   min((fi(i,j+1) - fi(i,j))/deltaY,0)^2);
           elseif F(i,j) < 0
               sigmaFi(i,j) = deltaT*F(i,j)*sqrt( ...
                   min((fi(i,j) - fi(i-1,j))/deltaX,0)^2 + ...
                   max((fi(i+1,j) - fi(i,j))/deltaX,0)^2 + ...
                   min((fi(i,j) - fi(i,j-1))/deltaY,0)^2 + ...
                   max((fi(i,j+1) - fi(i,j))/deltaY,0)^2);
           end    
        end
    end        
    
    %renovamos los valores internos, en el borde hay ceros, no cambia nada
    fi = fi - sigmaFi;
    
    %ahora hacemos interpolación lineal para los valores del borde
    fi(1,:) = 2*fi(2,:) - fi(3,:);
    fi(m,:) = 2*fi(m-1,:) - fi(m-2,:);
    fi(:,1) = 2*fi(:,2) - fi(:,3);
    fi(:,n) = 2*fi(:,n-1) - fi(:,n-2);
    
    %esquinas
    fi(1,1) = fi(2,1) + fi(1,2) - fi(2,2);
    fi(1,n) = fi(2,n) + fi(1,n-1) - fi(2,n-1);
    fi(m,1) = fi(m,2) + fi(m-1,1) - fi(m-1,2);
    fi(1,n) = fi(m,n-1) + fi(m-1,n) - fi(m-1,n-1);


end