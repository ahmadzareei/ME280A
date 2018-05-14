function [p,u,error,EI] = fem_linear_EI_a(N)
%clear all;
%close all;

l=1.0; %Length of the domain
A1=1.0; %constant in the equation

K=sparse(N+1,N+1);
L=zeros(N+1,1);
p = (1/N)*[0:N]';
e=[[1:N]',[2:N+1]'];
boundary=[1,N+1];
rhs = @(x) -90*pi^2*sin(3*pi*x).*sin(36*pi*x.^3) + ...
        (10*sin(3*pi*x)+5).*(216*pi.*x.*cos(36*pi*x.^3)-11664*pi.^2.*x.^4.*sin(36*pi*x.^3)) ...
        + 6480 * pi^2 *x.^2 .* cos(3*pi*x).*cos(36*pi*x.^3);
    
u_true = @(x) (10*sin(3*pi*x) + 5).*sin(36*pi*x.^3) ;
up_true = @(x) 10*(3*pi*cos(3*pi*x).*sin(36*pi*x.^3) + 3*36*pi*x.^2 .* cos(36*pi*x.^3).*sin(3*pi*x)) + 5*3*36*pi*x.^2.*cos(36*pi*x.^3);
EI = zeros(N,1);

for element = 1:size(e,1)
    points = p(e(element,:));
    h=points(2)-points(1);
    K_stamp = A1/h*[1,-1;-1,1];
    f = @(w) rhs(points(1) + (w+1)*h/2).*(1-w)/2;
    g = @(w) rhs(points(1) + (w+1)*h/2).*(1+w)/2;
    L_stamp = [-my_guass(f,-1,1);-my_guass(g,-1,1)]*h/2;       
    K(e(element,:),e(element,:)) = K(e(element,:),e(element,:)) + K_stamp;
    L(e(element,:)) = L(e(element,:))+L_stamp;      
end

K(1,:)= 0; 
K(1,1) = 1; 
L(1) = 0;  
%L(N+1) = L(N+1) + 1;
K(N+1,:)= 0; 
K(N+1,N+1) = 1; 
L(N+1) = 0;  

u = K\L;
%u_exact = u_true(p);

error = 0;
error_denom = 0;
for element = 1:size(e,1)
    points = p(e(element,:));
    h=points(2)-points(1);    
    ff = @(w) up_true(points(1) + (w+1)*h/2);
    gg = @(w) -u(element)/h + u(element+1)/h;
    EI(element) = quad(@(w) A1*(ff(w)-gg(w)).^2,-1,1)*h/2;
    error = error + EI(element);
    error_denom = error_denom + quad(@(w) A1*ff(w).^2,-1,1)*h/2;
    EI(element) = EI(element)/h;
end
error = sqrt(error)/sqrt(error_denom);
EI = EI *l/error_denom;
