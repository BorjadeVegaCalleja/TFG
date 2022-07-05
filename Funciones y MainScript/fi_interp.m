function  [fiInterpolada] = fi_interp(fi,X,Y)

    fiInterpolada =@(xp,yp) interp2(X,Y,fi,xp,yp);

end