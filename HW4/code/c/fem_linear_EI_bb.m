function [points,u,error,EI] = fem_linear_EI_bb(N)
%clear all;
%close all;


l=1.0; %Length of the domain
A1 = @(x) (x<1/3).*0.2 + (x>=1/3)*2; ; %constant in the equation

K=zeros(N+1,N+1);
L=zeros(N+1,1);

%points = (1/N)*[0:N]';
points = (1/N)*[0:N]';

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




KK=sparse(K);
KK(1,:)= 0; 
KK(1,1) = 1;
LL=L;
LL(1) = 0;  
LL(N+1) = LL(N+1) + 1;

u = KK\LL;

error = 0;
error_denom = 0;
for element = 1:size(e,1)
    p = points(e(element,:));
    h=p(2)-p(1);    
    ff =@(w) up_true(p(1) + (w+1)*h/2);
    gg = @(w) -u(element)/h + u(element+1)/h;
    EI(element) = quad(@(w) A1(p(1)+(w+1)*h/2).*(ff(w)-gg(w)).^2,-1,1)*h/2;
    error = error + EI(element);
    error_denom = error_denom + my_guass(@(w) A1(p(1)+(w+1)*h/2).*ff(w).^2,-1,1)*h/2;
    EI(element) = EI(element)/h;
end
error = sqrt(error)/sqrt(error_denom);
EI = EI*l/error_denom;


tol_EI = 0.04;
%tol_EI = 0.05^2;
max_EI = max(EI);    

while max_EI > tol_EI 

    N_element = size(e,1);
    for element = 1:N_element
        if(EI(element) > tol_EI)            
            elem = e(element,:);
            p = points(elem);
            h=p(2)-p(1);
            %zero out the entries for the matrix that we had
            K11 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi1p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
            K12 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
            K22 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi2p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
    
            L_stamp = [-my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi1(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1);...
                       -my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi2(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1)];
            K_stamp = [K11,K12;K12,K22];
            K(elem,elem)= K(elem,elem) - K_stamp;%[0,0;0,0]; ; %
            L(elem) = L(elem)- L_stamp; %[0;0]; %
            %making the new node
            points(end+1) = (p(2)+p(1))/2;
            p_size = length(points);
            e(element,:) = [ elem(1), p_size];
            e(end+1,:) = [p_size, elem(2)];
            %setting the matrix for the first element
            elem = e(element,:);
            p = points(elem);
            h=p(2)-p(1);
            K11 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi1p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
            K12 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
            K22 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi2p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
    
            L_stamp = [-my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi1(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1);...
                       -my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi2(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1)];
            K_stamp = [K11,K12;K12,K22];
            K = [K,zeros(size(K,1),1);zeros(1,size(K,1)+1)];
            K(elem,elem) = K(elem,elem) + K_stamp;
            L = [L;0];
            L(elem) = L(elem)+L_stamp;
            %setting the matrix for the second element
            elem = e(end,:);
            p = points(elem);
            h=p(2)-p(1);
            K11 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi1p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
            K12 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi1p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
            K22 = my_guass(@(x) A1(p(1)*phi1(x) + p(2)*phi2(x)).*phi2p(x).*phi2p(x)./(p(1)*phi1p(x) + p(2)*phi2p(x)),-1,1);
    
            L_stamp = [-my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi1(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1);...
                       -my_guass(@(w) rhs(p(1)*phi1(w) + p(2)*phi2(w)).*phi2(w).*(p(1)*phi1p(w) + p(2)*phi2p(w)) ,-1,1)];
            K_stamp = [K11,K12;K12,K22];
            K(elem,elem) = K(elem,elem) + K_stamp;
            L(elem) = L(elem)+L_stamp;
        end
    end


    KK=sparse(K);
    LL=L;
    KK(1,:)= 0; 
    KK(1,1) = 1; 
    LL(1) = 0;  
    LL(N+1) = LL(N+1) + 1;
%    KK(N+1,:)= 0;         
%    KK(N+1,N+1) = 1; 
    u = KK\LL;
    
    error = 0;
    EI = zeros(size(e,1),1);

    for element = 1:size(e,1)
        p = points(e(element,:));
        h=p(2)-p(1);    
        ff =@(w) up_true(p(1) + (w+1)*h/2);
        gg = @(w) -u(e(element,1))/h + u(e(element,2))/h;
        EI(element) = my_guass(@(w) A1(p(1)+(w+1)*h/2).*(ff(w)-gg(w)).^2,-1,1)*h/2;
        error = error + EI(element);
        error_denom = error_denom + my_guass(@(w) A1(p(1)+(w+1)*h/2).*ff(w).^2,-1,1)*h/2;
        EI(element) = EI(element)/h;
    end
    error = sqrt(error)/sqrt(error_denom);
    EI = EI *l/error_denom;
    max_EI = max(EI);
    
end

pp=[points,u];
pp = sortrows(pp,1);
points = pp(:,1);
u = pp(:,2);
ppp=[points(e(:,1)),EI];
ppp = sortrows(ppp,1);
EI = ppp(:,2);
