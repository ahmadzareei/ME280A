function result = my_guass1d(f,a,b)
%compute integral of f(x,y) in the domain [a,b];
%where [a,b] are the only parameters in the domain

xX = [0.0,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717];
wx = [0.0,0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881];

xm = 0.5*(b+a);
xr = 0.5*(b-a);


result = 0;
for i=2:6
        result = result + xr *wx(i).*(f(xm+xr*xX(i)) + f(xm-xr*xX(i)) );
end