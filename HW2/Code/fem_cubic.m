function [u,error] = fem_cubic(N)
%clear all;
%close all;

l=1.0; %Length of the domain
A1=0.2; %constant in the equation

K=sparse(3*N+1,3*N+1);
L=zeros(3*N+1,1);
points = (1/(3*N))*[0:3*N]';
e=[[1:3:3*N]',[2:3:3*N+1]',[3:3:3*N+1]',[4:3:3*N+1]'];
boundary=[1,3*N+1];
phi1 = @(x) (x+1/3).*(x-1/3).*(x-1)./((-1+1/3).*(-1-1/3).*(-1-1));
phi2 = @(x) (x+1).*(x-1/3).*(x-1)./((-1/3+1).*(-1/3-1/3).*(-1/3-1));
phi3 = @(x) (x+1).*(x+1/3).*(x-1)./((1/3+1).*(1/3+1/3).*(1/3-1));
phi4 = @(x) (x+1/3).*(x-1/3).*(x+1)./((1+1/3).*(1-1/3).*(1+1));
phi1p = @(x) (9*x)/8 - (27*x.^2)/16 + 1/16;
phi2p = @(x) (81*x.^2)/16 - (9*x)/8 - 27/16;
phi3p = @(x) 27/16 - (81*x.^2)/16 - (9*x)/8;
phi4p = @(x) (27*x.^2)/16 + (9*x)/8 - 1/16;
for element = 1:size(e,1)
    p = points(e(element,:));
    h=p(4)-p(1);

    K11 = my_guass(@(x) phi1p(x).*phi1p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K12 = my_guass(@(x) phi1p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K13 = my_guass(@(x) phi1p(x).*phi3p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K14 = my_guass(@(x) phi1p(x).*phi4p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);  
    K22 = my_guass(@(x) phi2p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K23 = my_guass(@(x) phi2p(x).*phi3p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K24 = my_guass(@(x) phi2p(x).*phi4p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K33 = my_guass(@(x) phi3p(x).*phi3p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K34 = my_guass(@(x) phi3p(x).*phi4p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);
    K44 = my_guass(@(x) phi4p(x).*phi4p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x) + p(3)*phi3p(x) + p(4)*phi4p(x)),-1,1);    
    K_stamp = A1*[K11,K12,K13,K14;K12,K22,K23,K24;K13,K23,K33,K34;K14,K24,K34,K44];

    f = @(w) (256*sin(0.75*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w)+ p(4)*phi4(w))).* ...
                  cos(16*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w)+ p(4)*phi4(w))));
    L_stamp = [-my_guass(@(w) f(w).*phi1(w).*(p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w) + p(4)*phi4p(w)),-1,1);...
               -my_guass(@(w) f(w).*phi2(w).*(p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w) + p(4)*phi4p(w)),-1,1);...
               -my_guass(@(w) f(w).*phi3(w).*(p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w) + p(4)*phi4p(w)),-1,1);...
               -my_guass(@(w) f(w).*phi4(w).*(p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w) + p(4)*phi4p(w)),-1,1)];       
    K(e(element,:),e(element,:)) = K(e(element,:),e(element,:)) + K_stamp;
    L(e(element,:)) = L(e(element,:))+L_stamp;      
end

K(1,:)= 0; 
K(1,1) = 1; 
L(1) = 0;  
L(3*N+1) = L(3*N+1) + 1;
u = K\L;


%exact solution repetition
A1= 0.2;
f =@(x) -128/(A1 *16.75^2 *pi^2) *sin (16.75 *pi*x) + ...
    128/(A1 * 15.25^2 *pi^2) *sin(15.25*pi*x) + 1/A1 * (3072/(4087*pi*sqrt(2))+1)*x;
%keyboard
%plot(p,u,'-*'); hold on; plot(p,uu,'.- b');
% 
%computing the error

error = 0;
error_denom = 0;
for element = 1:size(e,1)
    p = points(e(element,:));
    h=p(4)-p(1);    
    ff =@(w) -128/(A1 *16.75 *pi) *cos (16.75 *pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w)+p(4)*phi4(w))) + ...
              128/(A1 * 15.25 *pi) *cos(15.25*pi*(p(1)*phi1(w) + p(2)*phi2(w) + p(3)*phi3(w)++p(4)*phi4(w))) + ...
              1/A1 * (3072/(4087*pi*sqrt(2))+1);
    %ff = @(w) (-l*k/(A1 * pi)*cos(pi*k*(points(1) + (w+1)*h/2)/l) + (points(1) + (w+1)*h/2).^2/A1 + C);
    gg = @(w) (u(e(element,1))*phi1p(w) + u(e(element,2))*phi2p(w) + u(e(element,3))*phi3p(w)+ u(e(element,4))*phi4p(w))./ ...
                (p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w)+ p(4)*phi4p(w));
    error = error + my_guass(@(w) A1*(ff(w)-gg(w)).^2 .* (p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w)+p(4)*phi4p(w)),-1,1);
    error_denom = error_denom + my_guass(@(w) A1*ff(w).^2 .* (p(1)*phi1p(w) + p(2)*phi2p(w) + p(3)*phi3p(w)+ p(4)*phi4p(w)) ,-1,1);
end
error = sqrt(error)/sqrt(error_denom);