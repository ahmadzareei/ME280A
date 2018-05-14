function time_dependent_lumped_penalty

ri = 0.1;
ro = 0.25;

NEr = 10;
NEth = 20;

rho = 2;
cp = 2;
dt = 0.01;

[p,e] = mesh (ri, ro, NEr, NEth);

N = size(p,1);
K = sparse(N,N);
L = sparse(N,N);
Linv = sparse(N,N);
R = zeros(N,1);
k = @(x,y) 0.04;
qn  = -25;
z = @(X,Y) -500;
Tbar = 100;
T = 100;

%x and y are in master elements
%shape functions
phi1 = @(x,y) 1/4 *(x-1).*(y-1);
phi2 = @(x,y) -1/4 *(x+1).*(y-1);
phi3 = @(x,y) 1/4 *(x+1).*(y+1);
phi4 = @(x,y) -1/4 *(x-1).*(y+1);
%shape functions x-derivatives
phi1x = @(x,y) 1/4 *(y-1);
phi2x = @(x,y) -1/4 *(y-1);
phi3x = @(x,y) 1/4 *(y+1);
phi4x = @(x,y) -1/4 *(y+1);
%shape functions y-derivatives
phi1y = @(x,y) 1/4 *(x-1);
phi2y = @(x,y) -1/4 *(x+1);
phi3y = @(x,y) 1/4 *(x+1);
phi4y = @(x,y) -1/4 *(x-1);


for element =1:size(e,1)
    nodes = p(e(element,:),:);
    % xx and yy are elements in real space     
    xx = @(x,y)  nodes(1,1).*phi1(x,y) + nodes(2,1).*phi2(x,y) + nodes(3,1).*phi3(x,y) + nodes(4,1).*phi4(x,y);
    yy = @(x,y)  nodes(1,2).*phi1(x,y) + nodes(2,2).*phi2(x,y) + nodes(3,2).*phi3(x,y) + nodes(4,2).*phi4(x,y);
    F11 = @(x,y) nodes(1,1)*phi1x(x,y) + nodes(2,1)*phi2x(x,y) + nodes(3,1)*phi3x(x,y) + nodes(4,1)*phi4x(x,y); 
    F12 = @(x,y) nodes(1,1)*phi1y(x,y) + nodes(2,1)*phi2y(x,y) + nodes(3,1)*phi3y(x,y) + nodes(4,1)*phi4y(x,y);
    F21 = @(x,y) nodes(1,2)*phi1x(x,y) + nodes(2,2)*phi2x(x,y) + nodes(3,2)*phi3x(x,y) + nodes(4,2)*phi4x(x,y); 
    F22 = @(x,y) nodes(1,2)*phi1y(x,y) + nodes(2,2)*phi2y(x,y) + nodes(3,2)*phi3y(x,y) + nodes(4,2)*phi4y(x,y);
    F = @(x,y) [F11(x,y), F12(x,y);F21(x,y),F22(x,y)];
    detF = @(x,y) F11(x,y).*F22(x,y) - F12(x,y) .* F21(x,y);
    Finv = @(x,y) inv(F(x,y));    
    M = @(x,y) Finv(x,y)*transpose(Finv(x,y));
    Mp= @(x,y) Finv(x,y)*transpose(Finv(x,y));
    K11 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi1x(x,y);phi1y(x,y)]) * [phi1x(x,y);phi1y(x,y)]) ;
    K12 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi1x(x,y);phi1y(x,y)]) * [phi2x(x,y);phi2y(x,y)]) ;
    K13 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi1x(x,y);phi1y(x,y)]) * [phi3x(x,y);phi3y(x,y)]) ;
    K14 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi1x(x,y);phi1y(x,y)]) * [phi4x(x,y);phi4y(x,y)]) ;
    K22 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi2x(x,y);phi2y(x,y)]) * [phi2x(x,y);phi2y(x,y)]) ;
    K23 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi2x(x,y);phi2y(x,y)]) * [phi3x(x,y);phi3y(x,y)]) ;
    K24 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi2x(x,y);phi2y(x,y)]) * [phi4x(x,y);phi4y(x,y)]) ;
    K33 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi3x(x,y);phi3y(x,y)]) * [phi3x(x,y);phi3y(x,y)]) ;
    K34 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi3x(x,y);phi3y(x,y)]) * [phi4x(x,y);phi4y(x,y)]) ;
    K44 = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y) .* ( transpose(M(x,y)*[phi4x(x,y);phi4y(x,y)]) * [phi4x(x,y);phi4y(x,y)]) ;
    K_stamp = [ my_guass(@(x,y) K11(x,y),-1,1,-1,1) , my_guass(@(x,y) K12(x,y),-1,1,-1,1), my_guass(@(x,y) K13(x,y),-1,1,-1,1), my_guass(@(x,y) K14(x,y),-1,1,-1,1) ...
                ;0,my_guass(@(x,y) K22(x,y),-1,1,-1,1),my_guass(@(x,y) K23(x,y),-1,1,-1,1),my_guass(@(x,y) K24(x,y),-1,1,-1,1)...
                ;0,0,my_guass(@(x,y) K33(x,y),-1,1,-1,1),my_guass(@(x,y) K34(x,y),-1,1,-1,1)...
                ;0,0,0,my_guass(@(x,y) K44(x,y),-1,1,-1,1)];
    K_stamp(2,1) = K_stamp(1,2);
    K_stamp(3,1) = K_stamp(1,3);
    K_stamp(4,1) = K_stamp(1,4);
    K_stamp(3,2) = K_stamp(2,3);
    K_stamp(4,2) = K_stamp(2,4);    
    K_stamp(4,3) = K_stamp(3,4);
    K(e(element,:),e(element,:)) = K(e(element,:),e(element,:)) + K_stamp;
    R_stamp = [my_guass(@(x,y) z(xx(x,y),yy(x,y)).*phi1(x,y).*detF(x,y),-1,1,-1,1);...
               my_guass(@(x,y) z(xx(x,y),yy(x,y)).*phi2(x,y).*detF(x,y),-1,1,-1,1);...
               my_guass(@(x,y) z(xx(x,y),yy(x,y)).*phi3(x,y).*detF(x,y),-1,1,-1,1);...
               my_guass(@(x,y) z(xx(x,y),yy(x,y)).*phi4(x,y).*detF(x,y),-1,1,-1,1)];
    if( abs(nodes(2,1).^2 + nodes(2,2).^2 - ro.^2) < 0.00001)
        R_stamp(2) = R_stamp(2) + my_guass1d(@(y) sqrt(([1,0]*(Mp(1,y)*[1;0]))).*qn .* phi2(1,y).*detF(1,y),-1,1);
        R_stamp(3) = R_stamp(3) + my_guass1d(@(y) sqrt(([1,0]*(Mp(1,y)*[1;0]))).*qn .* phi3(1,y).*detF(1,y),-1,1);       
    end
    R(e(element,:)) = R(e(element,:)) + R_stamp;
    
    
    L11 = @(x,y) rho*cp/dt * detF(x,y) .* phi1(x,y) * phi1(x,y);
    L12 = @(x,y) rho*cp/dt * detF(x,y) .* phi1(x,y) * phi2(x,y);
    L13 = @(x,y) rho*cp/dt * detF(x,y) .* phi1(x,y) * phi3(x,y);
    L14 = @(x,y) rho*cp/dt * detF(x,y) .* phi1(x,y) * phi4(x,y);
    L22 = @(x,y) rho*cp/dt * detF(x,y) .* phi2(x,y) * phi2(x,y);
    L23 = @(x,y) rho*cp/dt * detF(x,y) .* phi2(x,y) * phi3(x,y);
    L24 = @(x,y) rho*cp/dt * detF(x,y) .* phi2(x,y) * phi4(x,y);
    L33 = @(x,y) rho*cp/dt * detF(x,y) .* phi3(x,y) * phi3(x,y);
    L34 = @(x,y) rho*cp/dt * detF(x,y) .* phi3(x,y) * phi4(x,y);
    L44 = @(x,y) rho*cp/dt * detF(x,y) .* phi4(x,y) * phi4(x,y);
    
    L_stamp = [ my_guass(@(x,y) L11(x,y),-1,1,-1,1) , my_guass(@(x,y) L12(x,y),-1,1,-1,1), my_guass(@(x,y) L13(x,y),-1,1,-1,1), my_guass(@(x,y) L14(x,y),-1,1,-1,1) ...
                ;0,my_guass(@(x,y) L22(x,y),-1,1,-1,1),my_guass(@(x,y) L23(x,y),-1,1,-1,1),my_guass(@(x,y) L24(x,y),-1,1,-1,1)...
                ;0,0,my_guass(@(x,y) L33(x,y),-1,1,-1,1),my_guass(@(x,y) L34(x,y),-1,1,-1,1)...
                ;0,0,0,my_guass(@(x,y) L44(x,y),-1,1,-1,1)];
    L_stamp(2,1) = L_stamp(1,2);
    L_stamp(3,1) = L_stamp(1,3);
    L_stamp(4,1) = L_stamp(1,4);
    L_stamp(3,2) = L_stamp(2,3);
    L_stamp(4,2) = L_stamp(2,4);    
    L_stamp(4,3) = L_stamp(3,4);
    L(e(element,:),e(element,:)) = L(e(element,:),e(element,:)) + L_stamp;
    
end

%Direct Dirichelet Implementation
for element =1:size(e,1)
    nodes = p(e(element,:),:);
    if( abs(nodes(1,1).^2 + nodes(1,2).^2 - ri.^2) < 0.00001)
        K(e(element,1),:) = 0;
        K(e(element,4),:) = 0;
        K(e(element,1),e(element,1)) = 1;
        K(e(element,4),e(element,4)) = 1;
        R(e(element,1)) = Tbar;
        R(e(element,4)) = Tbar;
    end   
end

%Sol = K\R;
Sol = 100*ones(size(R));
%Lumped Matrix Assumption
Linv = diag(sum(L'));
Linv = inv(Linv);
t=0;
clf;
patch('Faces',e,'Vertices',[p,Sol],'FaceVertexCData',Sol,'FaceCol','interp');colorbar; axis square
set(gca,'FontSize',16);
xlabel('x','FontSize',16);
ylabel('y','FontSize',16);
axis([-0.25,0.25,-0.25,0.25]);
title(sprintf('T(x,y) - Numerical Solution at time %5.2f',t),'Fontsize',16);
nam1  = sprintf('penalty_sol_start_%d.eps',floor(t/dt));
print('-depsc',nam1);

for t=0:dt:100
    Sol = Linv*L*Sol - Linv *K*Sol + Linv*R ;
%     for element =1:size(e,1)
%         nodes = p(e(element,:),:);
%         if( abs(nodes(1,1).^2 + nodes(1,2).^2 - ri.^2) < 0.00001)
%             Sol(e(element,1)) = Tbar;
%             Sol(e(element,4)) = Tbar;
%         end   
%     end

    if(mod(floor(t/dt),50) == 0)
        clf;
        patch('Faces',e,'Vertices',[p,Sol],'FaceVertexCData',Sol,'FaceCol','interp');colorbar; axis square
        set(gca,'FontSize',16);
        xlabel('x','FontSize',16);
        ylabel('y','FontSize',16);
        axis([-0.25,0.25,-0.25,0.25]);
        title(sprintf('T(x,y) - Numerical Solution at time %5.2f',t),'Fontsize',16);
        nam1  = sprintf('penalty_sol_%d.eps',floor(t/dt));
        print('-depsc',nam1);
    end
end

keyboard

