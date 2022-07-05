

function [t,pSuma,pCuadrado,pFi,pfixFi] = triangulationPersson(p,fdCuadrado,FI,h0,fh,pfixCuadrado,X,Y,bbox)

    fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);    

    [pCuadrado, ~] = distmesh2d(p,fdCuadrado,fh,h0,pfixCuadrado);
    pfixFi = findfixpoints(pCuadrado,fiInterp,h0);
    [pFi, ~] = distmesh2dFI(p,FI,fh,h0,pfixFi,X,Y);


    
    pFi = pFi(fiInterp(pFi(:,1),pFi(:,2))<-h0/10,:);  
    
    
    pSuma = [pCuadrado' pFi']';
    pSuma = back2boundingbox(pSuma,bbox,h0); %Lleva al borde los puntos muy cercanos al él
    pSuma = unique(pSuma,'rows');
    

    t = delaunayn(pSuma);
 

end    
