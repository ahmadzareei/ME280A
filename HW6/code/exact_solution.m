function exact_solution
figure
ro = 0.25;
ri = 0.10;
Tb = 100;
qb = -25;
z = -500;
k = 0.04;
N = 200;


dr = (ro-ri)/N;
dth = 2*pi/ N;
C1 = qb/k * ro + z/(2*k)*ro^2;
C2 = Tb + z/(4*k) *ri^2 - C1*log(ri);

[r,th] = meshgrid(ri:dr:ro,0:dth:2*pi);

f = -z/(4*k)*r.^2 + C1 * log(r) + C2;

surf(r.*cos(th), r.*sin(th), f);

view(0,90);
shading interp
colorbar
set(gca,'FontSize',16);
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
axis([-0.25,0.25,-0.25,0.25]);

title('T(x,y) - Exact Solution','Fontsize',20);
print -depsc exactsolution.eps