function result = my_guass(f,a,b,c,d)
%compute integral of f(x,y) in the domain [a,b] * [c,d];
%where [a,b] are in the x direction
%  and [c,d] are in the y direction
xX = [0.0,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717];
yY = [0.0,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717];
wx = [0.0,0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881];
wy = [0.0,0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881];

xm = 0.5*(b+a);
xr = 0.5*(b-a);

ym = 0.5*(d+c);
yr = 0.5*(d-c);

%dummy = @(YY) xr *sum((wx.*(f(xm+xr*xX,YY) + f(xm-xr*xX,YY))));
%result = yr *sum((wy.*(dummy(ym+yr*yY) + dummy(ym-yr*yY))));

result = 0;
for i=2:6
    for j = 2:6
        result = result + xr *wx(i).*(f(xm+xr*xX(i),ym+yr*yY(j)) + f(xm-xr*xX(i),ym+yr*yY(j)) +...
                                      f(xm+xr*xX(i),ym-yr*yY(j)) + f(xm-xr*xX(i),ym-yr*yY(j)) )* yr *wy(j);
    end
end