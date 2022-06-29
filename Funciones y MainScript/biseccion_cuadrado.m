function [result] = biseccion_cuadrado(x,deltaX,fd,h0)

    l = sqrt(deltaX(1)^2 + deltaX(2)^2);
    a = x;
    b = x + deltaX;
    tol = h0/1000;

    while l > tol
        m = (a+b)/2;
        sig = feval(fd,m);
        if sig > 0
            b = m;
        else
            a = m;
        end
        l = l/2;
    end
    result = a - x;

end