function [u,error] = fem_linear(N)
%clear all;
%close all;

l=1.0; %Length of the domain
A1=0.2; %constant in the equation

K=sparse(N+1,N+1);
L=zeros(N+1,1);
p = (1/N)*[0:N]';
e=[[1:N]',[2:N+1]'];
boundary=[1,N+1];

for element = 1:size(e,1)
    points = p(e(element,:));
    h=points(2)-points(1);
    K_stamp = A1/h*[1,-1;-1,1];
    f = @(w) (256*sin(0.75*pi*(points(1) + (w+1)*h/2)).*cos(16*pi*(points(1) + (w+1)*h/2))).*(1-w)/2;
    g = @(w) (256*sin(0.75*pi*(points(1) + (w+1)*h/2)).*cos(16*pi*(points(1) + (w+1)*h/2))).*(1+w)/2; 
    L_stamp = [-my_guass(f,-1,1);-my_guass(g,-1,1)]*h/2;       
    K(e(element,:),e(element,:)) = K(e(element,:),e(element,:)) + K_stamp;
    L(e(element,:)) = L(e(element,:))+L_stamp;      
end

K(1,:)= 0; 
K(1,1) = 1; 
L(1) = 0;  
L(N+1) = L(N+1) + 1;
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
    points = p(e(element,:));
    h=points(2)-points(1);    
    ff =@(w) -128/(A1 *16.75 *pi) *cos (16.75 *pi*(points(1) + (w+1)*h/2)) + ...
            128/(A1 * 15.25 *pi) *cos(15.25*pi*(points(1) + (w+1)*h/2)) + ...
            1/A1 * (3072/(4087*pi*sqrt(2))+1);
    %ff = @(w) (-l*k/(A1 * pi)*cos(pi*k*(points(1) + (w+1)*h/2)/l) + (points(1) + (w+1)*h/2).^2/A1 + C);
    gg = @(w) -u(element)/h + u(element+1)/h;
    error = error + quad(@(w) A1*(ff(w)-gg(w)).^2,-1,1)*h/2;
    error_denom = error_denom + quad(@(w) A1*ff(w).^2,-1,1)*h/2;
end
error = sqrt(error)/sqrt(error_denom);
