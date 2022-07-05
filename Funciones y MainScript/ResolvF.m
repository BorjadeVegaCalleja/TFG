function [result] = ResolvF(u,lambda,FI,pSuma,t,x_cart,y_cart,X,Y,h0)

    FI_fundament = [1 1/2 1/2;1/2 1 1/2;1/2 1/2 1]/12;
    E = [1 1/2;1/2 1]/3;
    fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);
    n = size(t,1);
    integral = 0.0; 
    for i = 1:n     
        
        x1 = pSuma(t(i,1),:);
        x2 = pSuma(t(i,2),:);
        x3 = pSuma(t(i,3),:);

        x = (x1+x2+x3)/3;    

        if fiInterp(x(1,1),x(1,2)) < 0 
            ro = 1;
        else
            ro = 2;
        end
        U = u(t(i,:))*u(t(i,:))'; 
        G = [x2-x1;x3-x1]; 
        det = G(1,1)*G(2,2)-G(1,2)*G(2,1);
        jac = abs(det); 
        integral = integral + (ro * jac * sum(sum(U.*FI_fundament)));
    end
    A = lambda*(2-1)/integral;

    integral = 0.0;
    len = 0.0;
    edges = find_edges(t,pSuma,FI,h0,X,Y); 
    s = size(edges,1);
    for i = 1:s
        l = sqrt(sum((pSuma(edges(i,1),:)-pSuma(edges(i,2),:)).^2)); 
        len = len + l; 
        integral = integral + l * (u(edges(i,1))^2 + u(edges(i,2))^2)/2;
    end
    nu = -integral/len;
    uinterp = scatteredInterpolant(pSuma(:,1),pSuma(:,2),u);
    discretizacion = crear_Matriz(uinterp,x_cart,y_cart);
    U2 = discretizacion.^2; 
    result = A*(U2 + nu);

end