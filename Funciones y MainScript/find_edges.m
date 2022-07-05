function [result] = find_edges(t,p,FI,h0,X,Y)

    
   bars = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         
   order_bars = unique(sort(bars,2),'rows');            
   fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);
   
   tol = h0/100;
   siz = size(order_bars);
   rows = siz(1);
   s = 1;
    
   for k = 1:rows
        i = order_bars(k,1);
        j = order_bars(k,2); 
        if (abs(fiInterp(p(i,1),p(i,2)))<= tol && abs(fiInterp(p(j,1),p(j,2)))<= tol)
            result(s,:) = [i,j];
            s = s + 1;
        end
   end    
end