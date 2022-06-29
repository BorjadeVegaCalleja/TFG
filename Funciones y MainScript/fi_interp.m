function  [fiInterpolada] = fi_interp(fi,X,Y)

    %x = linspace(0,1,100); %genero un vector de -1 a 1 con 100 puntos
    %y = linspace(0,1.5,100);
    %x = linspace(-1,1,100); %genero un vector de -1 a 1 con 100 puntos
    %y = linspace(-1,1,100);
    %x = linspace(0,1.5,100); %genero un vector de -1 a 1 con 100 puntos
    %y = linspace(0,1.5,100);


    fiInterpolada =@(xp,yp) interp2(X,Y,fi,xp,yp);

end