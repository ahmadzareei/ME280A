function [points,u,error,EI] = fem_linear_EI_b2(N)
%clear all;
%close all;
N1 = floor(N/3);
N2 = N - N1;

l=1.0; %Length of the domain
A1 = @(x) (x<1/3).*0.2 + (x>=1/3)*2; ; %constant in the equation

K=sparse(N+1,N+1);
L=zeros(N+1,1);

%points = (1/N)*[0:N]';
points1 = linspace (0,1/3,N1+1);
points2 = linspace (1/3,1,N2+1);
points = [points1(1:end-1), points2]';

e=[[1:N]',[2:N+1]'];
boundary=[1,N+1];
EI = zeros(N,1);
rhs = @(x) 256 * sin(3/4*pi*x) .* cos(16*pi*x);
    
u_true = @(x) 512 ./ (4087*pi*A1(x)) .*(268*sin(61*pi*x/4)/61/pi - 244*sin(67*pi*x/4)/67/pi) + ...
              (x<1/3).*(5+ 7680*sqrt(2)/4087/pi).*x + ...
              (x>=1/3).*...
              (2304*sqrt(2) *(8210-768*sqrt(3) + ...
                4087*pi)/(16703569*pi^2) + 3/2 + (1/2 + 768*sqrt(2)/4087/pi).*x);
up_true = @(x) 512 ./ (4087*pi*A1(x)) .* (268*cos(61*pi*x/4)/4 - 244*cos(67*pi*x/4)/4) + ...
               (x<1/3).*(5+ 7680*sqrt(2)/4087/pi) + ...
              (x>=1/3).*(1/2 + 768*sqrt(2)/4087/pi);
                
phi1 = @(x) 0.5*(1-x);
phi2 = @(x) 0.5*(1+x);

phi1p = @(x) -0.5;
phi2p = @(x) +0.5;     
          
for element = 1:size(e,1)
    p = points(e(element,:));
    h=p(2)-p(1);

    K11 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi1p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
    K12 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
    K22 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi2p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
    
    L_stamp = [-my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi1(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1);...
               -my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi2(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1)];
    K_stamp = [K11,K12;K12,K22];

    K(e(element,:),e(element,:)) = K(e(element,:),e(element,:)) + K_stamp;
    L(e(element,:)) = L(e(element,:))+L_stamp;      
end

K(1,:)= 0; 
K(1,1) = 1; 
L(1) = 0;  
L(N+1) = L(N+1) + 1;

u = K\L;

error = 0;
error_denom = 0;
for element = 1:size(e,1)
    p = points(e(element,:));
    h=p(2)-p(1);    
    ff =@(w) up_true(p(1) + (w+1)*h/2);
    gg = @(w) -u(element)/h + u(element+1)/h;
    EI(element) = quad(@(w) A1(p(1)+(w+1)*h/2).*(ff(w)-gg(w)).^2,-1,1)*h/2;
    error = error + EI(element);
    error_denom = error_denom + quad(@(w) A1(p(1)+(w+1)*h/2).*ff(w).^2,-1,1)*h/2;
    EI(element) = EI(element)/h;
end
error = sqrt(error)/sqrt(error_denom);
EI = EI*l/error_denom;
