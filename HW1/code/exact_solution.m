A1 = 0.2;
L=1;
C = (1-L^3/3/A1)/L;
x = linspace(0,1,500);
k = 2.^[0,1,2,3,4,5];

figure
for kp = k
    u = -L^2/A1/pi^2*sin(pi*kp*x/L) + x.^3/3/A1 + C*x;
    plot(x,u); hold on;
end

    