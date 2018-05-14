
function result = my_guass(f,a,b)
x = [0.0,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717];
w = [0.0,0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881];

xm = 0.5*(b+a);
xr = 0.5*(b-a);

result = xr *sum((w.*(f(xm+xr*x) + f(xm-xr*x))));
