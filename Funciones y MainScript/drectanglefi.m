function d=drectanglefi(x,y,x1,x2,y1,y2)


d=-min(min(min(-y1+y,y2-y),-x1+x),x2-x);