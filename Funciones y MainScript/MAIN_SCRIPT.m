% Main script

% El caso que está es el de la elipse, para poner el de los dos
% rectángulos simplemente hay que cambiar fi por la siguiente:
% fi = @(x,y) dunion(drectanglefi(x,y,0.2,0.8,0.1,0.4),drectanglefi(x,y,0.2,0.8,1.1,1.4));


clear, close all

bbox=[0,0;1,1.5];

deltaX = 0.025;
deltaY = 0.025;
deltaT = (deltaX/10)*0.5;

x_cart = bbox(1,1):deltaX:bbox(2,1);  
y_cart = bbox(1,2):deltaY:bbox(2,2);  

[X,Y] = meshgrid(x_cart,y_cart);

h0 = 0.01;
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      
p=[x(:),y(:)];                                       
sizp = size(p,1);
j = 1;
for i = 1:sizp
    a = p(i,1);
    if (p(i,1) <= 1)
        newp(j,:) = p(i,:);
        j = j + 1;
    end
end
p = newp;
pfixCuadrado = [bbox(1,1),bbox(1,2);bbox(1,1),bbox(2,2);bbox(2,1),bbox(1,2);bbox(2,1),bbox(2,2)];
fh = @huniform;
fi = @(x,y)sqrt((x-0.5).^2*0.1^2/0.4^2+(y-0.75).^2) - 0.1;
FI = crear_Matriz(fi,x_cart,y_cart);
fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);
fdCuadrado = @(p) ddiff(drectangle(p,0,1,0,1.5),fiInterp(p(:,1),p(:,2)));
[t,pSuma,pCuadrado,pFi,pfixFi] = triangulationPersson(p,fdCuadrado,FI,h0,fh,pfixCuadrado,X,Y,bbox);
[u, lambda] = fem(t,pSuma,FI,X,Y);
lambda = lambda(1,1);
u = u(:,1);
lambda_list(1) = lambda; 
result = ResolvF(u,lambda,FI,pSuma,t,x_cart,y_cart,X,Y,h0);
areaS(1) = calculo_area(t,pSuma,FI,X,Y);
%%

for i = 1:1
 
    FI = evol(FI, result, deltaT, deltaX, deltaY);
    fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);
    fdCuadrado = @(pSuma) ddiff(drectangle(pSuma,0,1,0,1.5),fiInterp(pSuma(:,1),pSuma(:,2)));
    [t,pSuma,pCuadrado,pFi,pfixFi] = triangulationPersson(pSuma,fdCuadrado,FI,h0,fh,pfixCuadrado,X,Y,bbox);
   
    [u, lambda] = fem(t,pSuma,FI,X,Y);
    lambda = lambda(1,1);
    u = u(:,1);
    lambda_list(i+1) = lambda; 
    iteraciones = i+1

    result = ResolvF(u,lambda,FI,pSuma,t,x_cart,y_cart,X,Y,h0);
    areaS(i+1) = calculo_area(t,pSuma,FI,X,Y);

end
    figure(1)
    hold on
    triplot(t,pSuma(:,1),pSuma(:,2),'c');
    plot(pFi(:,1),pFi(:,2),'.black');
    plot(pfixFi(:,1),pfixFi(:,2),'g.')
    axis equal
    axis([0 1 0 1.5])
    hold off
    figure(2)
    surf(X,Y,FI)
    hold on
    surf(X,Y,zeros(length(y_cart),length(x_cart)))


