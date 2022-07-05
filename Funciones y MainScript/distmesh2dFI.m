function [p,t]=distmesh2dFI(p,FI,fh,h0,pfix,X,Y,varargin)

dptol=.01; ttol=.1; Fscale=1.2; deltat=.2; geps=.01*h0; deps=sqrt(eps)*h0;
densityctrlfreq=30;

fiInterp = @(xp,yp) interp2(X,Y,FI,xp,yp); % interpolamos FI
% 2. Eliminamos los puntos fuera de la región
p=p(fiInterp(p(:,1),p(:,2))<geps,:);                 
r0=1./feval(fh,p,varargin{:}).^2;                   
p=p(rand(size(p,1),1)<r0./max(r0),:);               
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end    
pfix=unique(pfix,'rows'); nfix=size(pfix,1);
p=[pfix; p];                                        
N=size(p,1);                                        

count=0;
pold=inf;                                           
while 1
  count=count+1;
  % 3. Triangulación con Delaunay
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol          
    pold=p;                                          
    t=delaunayn(p);                                  
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;  
    t=t(fiInterp(pmid(:,1),pmid(:,2))<-geps,:);        
   
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         
    bars=unique(sort(bars,2),'rows');               
  end

 
  barvec=p(bars(:,1),:)-p(bars(:,2),:);             
  L=sqrt(sum(barvec.^2,2));                          
  hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));    
  
  if mod(count,densityctrlfreq)==0 & any(L0>2*L)
      p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
      N=size(p,1); pold=inf;
      continue;
  end
  
  F=max(L0-L,0);                                    
  Fvec=F./L*[1,1].*barvec;                          
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(1:size(pfix,1),:)=0;                           
  
  paux=p+deltat*Ftot;                                 
  deltaP = paux - p;
  d=fiInterp(paux(:,1),paux(:,2)); ix_out=d>0;
  x = p(ix_out,:);
  deltaX = deltaP(ix_out,:);
  for i = 1:size(x,1)
        deltaX(i,:) = biseccion_fi(x(i,:),deltaX(i,:),fiInterp,h0);
  end
  deltaP(ix_out,:) = deltaX;
  p = p + deltaP;

  if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
end

