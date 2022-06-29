function [areaS] = calculo_area(t,p,FI,X,Y)
    
    sizt = size(t,1);
    fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);
    areaS = 0;
    for i = 1:sizt 
        pmid=(p(t(i,1),:)+p(t(i,2),:)+p(t(i,3),:))/3; 
        nodes = t(i,:);
		v = p(nodes,:); 
        x1 = v(1,:);
        x2 = v(2,:);
        x3 = v(3,:);
        G = [x2-x1;x3-x1]; 
        det = G(1,1)*G(2,2)-G(1,2)*G(2,1);
        area = abs(det)/2;
        if fiInterp(pmid(:,1),pmid(:,2)) < 0 
             areaS = areaS + area;
        else
            
        end
    end

end