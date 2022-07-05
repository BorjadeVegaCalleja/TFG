function [p] = back2boundingbox(p,bb,h0)
    
    tol = h0/100;
    p(p(:,1) < bb(1,1)+tol,1) = bb(1,1); 
    p(p(:,2) < bb(1,2)+tol,2) = bb(1,2);
    p(p(:,1) > bb(2,1)+tol,1) = bb(2,1);
    p(p(:,2) > bb(2,2)+tol,2) = bb(2,2);
    

end