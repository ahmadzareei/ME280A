function result = my_guass(f,a,b)

%x = [0.0, 0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717];
%w = [0.0, 0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881];

x=[0.0,0.538469310105683,0.906179845938664];
w=[0.568888888888889/2,0.478628670499366,0.236926885056189];


xm = 0.5*(b+a);
xr = 0.5*(b-a);
result = 0;
for i=2:3
%for i=2:5
    result = result + xr *((w(i).*(f(xm+xr*x(i)) + f(xm-xr*x(i)))));
end
result = result + xr*w(1)*f(xm+xr*x(i));
