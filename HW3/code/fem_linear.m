function [u,u_cg,error,J,iter] = fem_linear(N,tole)
%clear all;
%close all;

l=1.0; %Length of the domain
% A1=0.2; %constant in the equation
A1 = [2.00,2.5,1.25,.25,4.00,1.75,0.5,0.75,3.25,1.00];

K=sparse(N+1,N+1);
L=zeros(N+1,1);
p = (1/N)*[0:N]';
e=[[1:N]',[2:N+1]'];
boundary=[1,N+1];
Ke_table = zeros(N,3);
Le = zeros(N+1,1);

for element = 1:size(e,1)
    points = p(e(element,:));
    mid = (points(1) + points(2))/2;
    h=points(2)-points(1);
    K_stamp = A1(floor(mid*10)+1)/h*[1,-1;-1,1];
    f = @(w) (256*sin(0.75*pi*(points(1) + (w+1)*h/2)).*cos(16*pi*(points(1) + (w+1)*h/2))).*(1-w)/2;
    g = @(w) (256*sin(0.75*pi*(points(1) + (w+1)*h/2)).*cos(16*pi*(points(1) + (w+1)*h/2))).*(1+w)/2; 
    L_stamp = [-my_guass(f,-1,1);-my_guass(g,-1,1)]*h/2;       
    K(e(element,:),e(element,:)) = K(e(element,:),e(element,:)) + K_stamp;
    Ke_table(element,:) =Ke_table(element,:)+[K_stamp(1,1),K_stamp(1,2),K_stamp(2,2)];
    L(e(element,:)) = L(e(element,:))+L_stamp;      
end

K(1,:)= 0; 
K(1,1) = 1;

Le = L;
Ke_table (1,1) =1;
% Since we set both Ke_12 and Ke_21 to zero we need to update Le(2)
% However since u(0) is zero nothing happens
Le(2) = Le(2) - 0*Ke_table(1,2); % Le(2) - u0 * K_21;
Ke_table(1,2) = 0;

L(1) = 0;  
Le(1) = 0;

L(N+1) = L(N+1) + 1;
Le(N+1) = Le(N+1) + 1;

u = K\L;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %exact solution repetition
% [x,y,yp] = exact_solution(N);
% plot(x,y,p,u)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % testing my_multiply function
% kkk = [1:N+1]';
% Kp = sparse(N+1,N+1);
% Kp = K;
% Kp(2,1) = 0;
% %These two should be same
% Kp*kkk;
% my_multiply(Ke_table,kkk,e);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% finding the preconditioner for the Matrix K
T = zeros(N+1,1);
T(1,1) = 1./sqrt(Ke_table(1,1));
T(N+1,1) = 1./sqrt(Ke_table(N,3));
for i = 2:N
    T(i) = 1./sqrt( Ke_table(i,1) + Ke_table(i-1,3));
end

% Doing the preconditioning
Ke_table_pre = precondition(Ke_table,T,e);
Le_pre = T.*Le;
%finding the solution using the conjugate gradient that we wrote
[u_cg,iter] = my_conjugate_gradient(Ke_table_pre,Le_pre,zeros(N+1,1),tole,@my_multiply,e);
u_cg = T.*u_cg;


% % % Computing the error
error = 0;
error_denom = 0;
J=0;
for element = 1:size(e,1)
    points = p(e(element,:));
    h=points(2)-points(1);    
    index = floor(10*(points(1)+ h/2))+1;
    if (index > 10)
        index = 10;
    end
    ff =@(w) u_Np((points(1) + (w+1)*h/2));
    hh =@(w) u_N((points(1) + (w+1)*h/2));

    gg = @(w) -u_cg(element)/h + u_cg(element+1)/h;
    error = error + A1(index) * my_guass(@(w) (ff(w)-gg(w)).^2,-1,1)*h/2;
    error_denom = error_denom + A1(index)* my_guass(@(w) ff(w).^2,-1,1)*h/2;
    J = J + A1(index)* my_guass(@(w) hh(w).^2,-1,1)*h/2;
end
error = sqrt(error)/sqrt(error_denom);