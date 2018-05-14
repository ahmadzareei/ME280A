function fem_main1

Nmax = [1348,94,34];
A1= 0.2;
f =@(x) -128/(A1 *16.75^2 *pi^2) *sin (16.75 *pi*x) + 128/(A1 * 15.25^2 *pi^2) *sin(15.25*pi*x) + 1/A1 * (3072/(4087*pi*sqrt(2))+1)*x;
x=linspace(0,1,5000);

%Linear Element 
err = 0;
index = 1;
N=1;
error =1;
figure
uu = f(x);
plot(x,uu);
hold on;
while N < Nmax(1)
   N=2*N;
   [u,err(index)] = fem_linear(N);
   plot([0:N]/N,u);
   index = index+1;
end
N=Nmax(1);
[u,err(index)] = fem_linear(N);
plot([0:N]/N,u);    
hold off
M = [2.^(1:index-2),Nmax(1),2.^(index-1)];
err = [err(1:index-2),err(index), err(index-1)];
figure
loglog(M,err);
caption = sprintf('Linear Case');
title(caption);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quadratic Element 
zeta = [-1:0.1:1];
phi1 = @(x) 0.5*x.*(x-1);
phi2 = @(x) (1-x).*(x+1);
phi3 = @(x) 0.5*x.*(x+1);

err = 0;
index = 1;
N=2;
error =1;
figure
uu = f(x);
plot(x,uu);
hold on;

while N < Nmax(2)
    points = (1/(2*N))*[0:2*N]';
    e=[[1:2:2*N]',[2:2:2*N+1]',[3:2:2*N+1]'];
    boundary=[1,2*N+1];
    [u,err(index)] = fem_quadratic(N);
    xx = zeros(size(e,1)*length(zeta),1);
    yy = zeros(size(e,1)*length(zeta),1);
    for element = 1:size(e)
        p = points(e(element,:));
        a_i = u(e(element,:));
        xx((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*p(1) + phi2(zeta)*p(2) +phi3(zeta)*p(3);
        yy((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*a_i(1) + phi2(zeta)*a_i(2) +phi3(zeta)*a_i(3);
    end
    plot(xx,yy);
    N=2*N;
    index = index+1;
end
N=Nmax(2);
points = (1/(2*N))*[0:2*N]';
e=[[1:2:2*N]',[2:2:2*N+1]',[3:2:2*N+1]'];
boundary=[1,2*N+1];
[u,err(index)] = fem_quadratic(N);
    xx = zeros(size(e,1)*length(zeta),1);
    yy = zeros(size(e,1)*length(zeta),1);
    for element = 1:size(e)
        p = points(e(element,:));
        a_i = u(e(element,:));
        xx((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*p(1) + phi2(zeta)*p(2) +phi3(zeta)*p(3);
        yy((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*a_i(1) + phi2(zeta)*a_i(2) +phi3(zeta)*a_i(3);
    end
plot(xx,yy);
hold off
M = [2.^(1:index-2),Nmax(2),2.^(index-1)];
err = [err(1:index-2),err(index), err(index-1)];
figure
loglog(M,err);
caption = sprintf('Quadratic Element');
title(caption);
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cubic Element 
zeta = [-1:0.1:1];
phi1 = @(x) (x+1/3).*(x-1/3).*(x-1)./((-1+1/3).*(-1-1/3).*(-1-1));
phi2 = @(x) (x+1).*(x-1/3).*(x-1)./((-1/3+1).*(-1/3-1/3).*(-1/3-1));
phi3 = @(x) (x+1).*(x+1/3).*(x-1)./((1/3+1).*(1/3+1/3).*(1/3-1));
phi4 = @(x) (x+1/3).*(x-1/3).*(x+1)./((1+1/3).*(1-1/3).*(1+1));

err = 0;
index = 1;
N=2;
error =1;
figure
uu = f(x);
plot(x,uu);
hold on;

while N < Nmax(2)
    points = (1/(3*N))*[0:3*N]';
    e=[[1:3:3*N]',[2:3:3*N+1]',[3:3:3*N+1]',[4:3:3*N+1]'];
    [u,err(index)] = fem_cubic(N);
    xx = zeros(size(e,1)*length(zeta),1);
    yy = zeros(size(e,1)*length(zeta),1);
    for element = 1:size(e)
        p = points(e(element,:));
        a_i = u(e(element,:));
        xx((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*p(1) + phi2(zeta)*p(2) +phi3(zeta)*p(3)+phi4(zeta)*p(4);
        yy((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*a_i(1) + phi2(zeta)*a_i(2) +phi3(zeta)*a_i(3)+ phi4(zeta)*a_i(4);
    end
    plot(xx,yy);
    N=2*N;
    index = index+1;
end
N=Nmax(3);
points = (1/(3*N))*[0:3*N]';
e=[[1:3:3*N]',[2:3:3*N+1]',[3:3:3*N+1]',[4:3:3*N+1]'];
[u,err(index)] = fem_cubic(N);
xx = zeros(size(e,2)*length(zeta),1);
yy = zeros(size(e,2)*length(zeta),1);
for element = 1:size(e)
        p = points(e(element,:));
        a_i = u(e(element,:));
        xx((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*p(1) + phi2(zeta)*p(2) +phi3(zeta)*p(3)+phi4(zeta)*p(4);
        yy((element-1)*length(zeta)+1:element*length(zeta)) = phi1(zeta)*a_i(1) + phi2(zeta)*a_i(2) +phi3(zeta)*a_i(3)+ phi4(zeta)*a_i(4);
%        xx = phi1(zeta)*p(1) + phi2(zeta)*p(2) +phi3(zeta)*p(3)+phi4(zeta)*p(4);
%        yy = phi1(zeta)*a_i(1) + phi2(zeta)*a_i(2) +phi3(zeta)*a_i(3)+ phi4(zeta)*a_i(4);
end
plot(xx,yy);
hold off
M = [2.^(1:index-2),Nmax(3),2.^(index-1)];
err = [err(1:index-2),err(index), err(index-1)];
figure
loglog(M,err);
caption = sprintf('Quadratic Element');
title(caption);

keyboard

end
