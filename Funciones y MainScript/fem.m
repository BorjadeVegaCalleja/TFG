

% el model se crea con createpde() y se va rellenando con lo que necesite
% para crear la malla utilizamos geometryFromMesh y le pasamos el model, 
% los puntos y los triángulos
% Con applyBoundaryCondition aplicamos las condiciones para definir el tipo
% de problema que queremos resolver, en este caso es un problema de
% autovalores
function [u,lambda] = fem(t,p,FI,X,Y)

    model = createpde();
    
    geometryFromMesh(model,p',t');
    
    applyBoundaryCondition(model,'dirichlet', ...
                                 'Edge',1:model.Geometry.NumEdges, ...
                                 'u',0);

    ro = @(l,s) RO(l,FI,X,Y);
    
    specifyCoefficients(model,'m',0,'d',ro,'c',1,'a',0,'f',0);
    
    r = [7,20];
    
    results = solvepdeeig(model,r);
    u = results.Eigenvectors;
    lambda = results.Eigenvalues;
     

end

function d = RO(l,FI,X,Y)
    fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp);
    n = length(l.x);
    d = zeros(1,n);
    v = fiInterp(l.x,l.y);
    d(v<0) = 1;
    d(v>=0) = 2;
    
end

