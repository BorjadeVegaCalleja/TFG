function [a] = findfixpoints(p,fiInterp,h0)
    
    siz = size(p);
    rows = siz(1);
    j = 1;
    tol = h0/100;
    %fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);

    for i = 1:rows
        if (abs(fiInterp(p(i,1),p(i,2)))<= tol)
            a(j,:) = p(i,:);
            j = j+1;
        end    
    end  

end