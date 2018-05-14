function [p,u,error,EI] = fem_linear_EI_a(N)
%clear all;
%close all;

l=1.0; %Length of the domain
A1=1.0; %constant in the equation

K=zeros(N+1,N+1);

L=zeros(N+1,1);
p = (1/N)*[0:N]';
e=[[1:N]',[2:N+1]'];
boundary=[1,N+1];
rhs = @(x) -90*pi^2*sin(3*pi*x).*sin(36*pi*x.^3) + ...
        (10*sin(3*pi*x)+5).*(216*pi.*x.*cos(36*pi*x.^3)-11664*pi.^2.*x.^4.*sin(36*pi*x.^3)) ...
        + 6480 * pi^2 *x.^2 .* cos(3*pi*x).*cos(36*pi*x.^3);
    
u_true = @(x) (10*sin(3*pi*x) + 5).*sin(36*pi*x.^3) ;
up_true = @(x) 10*(3*pi*cos(3*pi*x).*sin(36*pi*x.^3) + 3*36*pi*x.^2 .* cos(36*pi*x.^3).*sin(3*pi*x)) + 5*3*36*pi*x.^2.*cos(36*pi*x.^3);
EI = ones(N,1);

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

L(1) = 0;  
%L(N+1) = L(N+1) + 1;
L(N+1) = 0;
KK = sparse(K);
KK(N+1,:)= 0; 
KK(N+1,N+1) = 1; 
KK(1,:)= 0; 
KK(1,1) = 1; 

u = KK\L;
%u_exact = u_true(p);

error = 0;
error_denom = 0;
for element = 1:size(e,1)
    points = p(e(element,:));
    h=points(2)-points(1);    
    ff = @(w) up_true(points(1) + (w+1)*h/2);
    gg = @(w) -u(element)/h + u(element+1)/h;
    EI(element) = my_guass(@(w) A1*(ff(w)-gg(w)).^2,-1,1)*h/2;
    error = error + EI(element);
    error_denom = error_denom + my_guass(@(w) A1*ff(w).^2,-1,1)*h/2;
    EI(element) = EI(element)/h;
end
error = sqrt(error)/sqrt(error_denom);
EI = EI *l/error_denom;


%tol_EI = 0.0148;
tol_EI = 0.05^2;

max_EI = max(EI);    

while max_EI > tol_EI 

    N_element = size(e,1);
    for element = 1:N_element
        if(EI(element) > tol_EI)            
            elem = e(element,:);
            points = p(elem);
            h=points(2)-points(1);
            %zero out the entries for the matrix that we had
            K_stamp = A1/h*[1,-1;-1,1];
            f = @(w) rhs(points(1) + (w+1)*h/2).*(1-w)/2;
            g = @(w) rhs(points(1) + (w+1)*h/2).*(1+w)/2;
            L_stamp = [-my_guass(f,-1,1);-my_guass(g,-1,1)]*h/2;
            K(elem,elem)= K(elem,elem) - K_stamp;%[0,0;0,0]; ; %
            L(elem) = L(elem)- L_stamp; %[0;0]; %
            %making the new node
            p(end+1) = (points(2)+points(1))/2;
            p_size = length(p);
            e(element,:) = [ elem(1), p_size];
            e(end+1,:) = [p_size, elem(2)];
            %setting the matrix for the first element
            elem = e(element,:);
            points = p(elem);
            h=points(2)-points(1);
            K_stamp = A1/h*[1,-1;-1,1];
            f = @(w) rhs(points(1) + (w+1)*h/2).*(1-w)/2;
            g = @(w) rhs(points(1) + (w+1)*h/2).*(1+w)/2;
            L_stamp = [-my_guass(f,-1,1);-my_guass(g,-1,1)]*h/2;       
            K = [K,zeros(size(K,1),1);zeros(1,size(K,1)+1)];
            K(elem,elem) = K(elem,elem) + K_stamp;
            L = [L;0];
            L(elem) = L(elem)+L_stamp;
            %setting the matrix for the second element
            elem = e(end,:);
            points = p(elem);
            h=points(2)-points(1);
            K_stamp = A1/h*[1,-1;-1,1];
            f = @(w) rhs(points(1) + (w+1)*h/2).*(1-w)/2;
            g = @(w) rhs(points(1) + (w+1)*h/2).*(1+w)/2;
            L_stamp = [-my_guass(f,-1,1);-my_guass(g,-1,1)]*h/2;       
            K(elem,elem) = K(elem,elem) + K_stamp;
            L(elem) = L(elem)+L_stamp;
        end
    end

    L(N+1) = 0;
    KK=sparse(K);
    KK(1,:)= 0; 
    KK(1,1) = 1; 
    L(1) = 0;  
    %L(N+1) = L(N+1) + 1;
    KK(N+1,:)= 0;         
    KK(N+1,N+1) = 1; 
    L(N+1) = 0;
    u = KK\L;
    
    error = 0;
    EI = zeros(size(e,1),1);

    for element = 1:size(e,1)
        points = p(e(element,:));
        h=points(2)-points(1);    
        ff = @(w) up_true(points(1) + (w+1)*h/2);
        gg = @(w) -u(e(element,1))/h + u(e(element,2))/h;
        EI(element) = my_guass(@(w) A1*(ff(w)-gg(w)).^2,-1,1)*h/2;
        error = error + EI(element);
        error_denom = error_denom + my_guass(@(w) A1*ff(w).^2,-1,1)*h/2;
        EI(element) = EI(element)/h;
    end
    error = sqrt(error)/sqrt(error_denom);
    EI = EI *l/error_denom;
    max_EI = max(EI);


end

pp=[p,u];
ppp=[p(e(:,1)),EI];
pp = sortrows(pp,1);
ppp = sortrows(ppp,1);
p = pp(:,1);
u = pp(:,2);
EI= ppp(:,2);
