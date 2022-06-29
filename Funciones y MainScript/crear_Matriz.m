
function Discretizacion = crear_Matriz(funcion,x,y)

Discretizacion = zeros(length(y),length(x));

for j = 1:length(x)
    for i = 1:length(y)
        Discretizacion(i,j) = funcion(x(j),y(i));
    end
end    
 
end
    

