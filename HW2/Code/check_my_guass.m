function check_my_guass(f,a,b)

x10 = [0.0,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717];
w10 = [0.0,0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881];

x2= [0.0,1/sqrt(3)];
w2= [0.0,1.0];

x3=[0,sqrt(3/5)];
w3=[0.888888888888889/2,0.555555555555556];

x4=[0.0,0.339981043584856,0.861136311594053];
w4=[0.0,0.652145154862546,0.347854845137454];

x5=[0.0,0.538469310105683,0.906179845938664];
w5=[0.568888888888889/2,0.478628670499366,0.236926885056189];

xm = 0.5*(b+a);
xr = 0.5*(b-a);

result10= xr*sum((w10.*(f(xm+xr*x10) + f(xm-xr*x10))));
result2 = xr*sum((w2.*(f(xm+xr*x2) + f(xm-xr*x2))));
result4 = xr*sum((w4.*(f(xm+xr*x4) + f(xm-xr*x4))));
result3 = xr*sum((w3.*(f(xm+xr*x3) + f(xm-xr*x3))));
result5 = xr*sum((w5.*(f(xm+xr*x5) + f(xm-xr*x5))));
result = quad(f,a,b,1e-16);

tol = 1e-15;

if(abs(result-result2)<tol)
    fprintf('Quad 2 passed\n');
else
        fprintf('Quad 2 failed\n');
end

if(abs(result-result3)<tol)
    fprintf('Quad 3 passed\n');
else
        fprintf('Quad 3 failed\n');
end

if(abs(result-result4)<tol)
    fprintf('Quad 4 passed\n');
else
        fprintf('Quad 4 failed\n');
end

if(abs(result-result5)<tol)
    fprintf('Quad 5 passed\n');
else
        fprintf('Quad 5 failed\n');
end

if(abs(result-result10)<tol)
    fprintf('Quad 10 passed\n');
else
        fprintf('Quad 10 failed\n');
end
