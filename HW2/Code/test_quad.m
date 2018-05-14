function test_quad

Nmax = [1348,94,34];

%testing linear case
fprintf('Linear Case test\n');
N=Nmax(1);
dx = 1./N;
a=floor(N/2)*dx;
b=a+dx;
points=[a,b];
h=dx;
f = @(w) (256*sin(0.75*pi*(points(1) + (w+1)*h/2)).*cos(16*pi*(points(1) + (w+1)*h/2))).*(1-w)/2;
check_my_guass(f,a,b)

%testing Quadratic case
fprintf('Quadratic Case test\n');
N=Nmax(2);
dx = 1./N;
a=floor(N/2)*dx;
b=a+dx/2;
c=a+dx;
p=[a,b,c];
h=dx;
phi1 = @(x) 0.5*x.*(x-1);
phi2 = @(x) (1-x).*(x+1);
phi3 = @(x) 0.5*x.*(x+1);
phi1p = @(x) 0.5*(2*x-1);
phi2p = @(x) (-2*x);
phi3p = @(x) 0.5*(2*x+1);
f = @(w) (256*sin(0.75*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w))).* ...
                  cos(16*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w))));
check_my_guass(f,a,b)
%g=@(x) phi1p(x).*phi1p(x);     
%check_my_guass(g,-1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%testing Cubic case
fprintf('Cubic Case test\n');
N=Nmax(3);
dx = 1./N;
a=0;
b=a+dx/3;
c=a+2*dx/3;
d=a+1*dx;
p=[a,b,c,d];
h=dx;
phi1 = @(x) (x+1/3).*(x-1/3).*(x-1)./((-1+1/3).*(-1-1/3).*(-1-1));
phi2 = @(x) (x+1).*(x-1/3).*(x-1)./((-1/3+1).*(-1/3-1/3).*(-1/3-1));
phi3 = @(x) (x+1).*(x+1/3).*(x-1)./((1/3+1).*(1/3+1/3).*(1/3-1));
phi4 = @(x) (x+1/3).*(x-1/3).*(x+1)./((1+1/3).*(1-1/3).*(1+1));
phi1p = @(x) (9*x)/8 - (27*x.^2)/16 + 1/16;
phi2p = @(x) (81*x.^2)/16 - (9*x)/8 - 27/16;
phi3p = @(x) 27/16 - (81*x.^2)/16 - (9*x)/8;
phi4p = @(x) (27*x.^2)/16 + (9*x)/8 - 1/16;

f = @(w) (256*sin(0.75*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w)+ p(4)*phi4(w))).* ...
                  cos(16*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w)+ p(4)*phi4(w))));
%g = @(x) phi1p(x).*phi1p(x);
check_my_guass(f,a,b)
%check_my_guass(g,-1,1)

