function [Sol, p, e,error] = fem_penalty_method(ri,ro,NEr, NEth,PenaltyConstant)

[p,e] = mesh (ri, ro, NEr, NEth);

N = size(p,1);
K = sparse(N,N);
R = zeros(N,1);
K_penalty = sparse(N,N);
R_penalty = zeros(N,1);

k = @(x,y) 0.04;
qn  = -25;
z = @(X,Y) -500;
Tbar = 100;

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
    %M = @(x,y) transpose(Finv(x,y))*Finv(x,y);
    M = @(x,y) Finv(x,y)*transpose(Finv(x,y));
    %Mp = @(x,y) transpose(Finv(x,y))*Finv(x,y);
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
        
    
   %Penalty part of the Stiffness matrix and Loading Vector
    if( abs(nodes(1,1).^2 + nodes(1,2).^2 - ri.^2) < 0.00001)
        Kp11 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi1(-1,y) .* phi1(-1,y).*detF(-1,y),1,-1);
        Kp12 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi1(-1,y) .* phi2(-1,y).*detF(-1,y),1,-1);
        Kp13 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi1(-1,y) .* phi3(-1,y).*detF(-1,y),1,-1);
        Kp14 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi1(-1,y) .* phi4(-1,y).*detF(-1,y),1,-1);
        Kp22 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi2(-1,y) .* phi2(-1,y).*detF(-1,y),1,-1);
        Kp23 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi2(-1,y) .* phi3(-1,y).*detF(-1,y),1,-1);
        Kp24 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi2(-1,y) .* phi4(-1,y).*detF(-1,y),1,-1);
        Kp33 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi3(-1,y) .* phi3(-1,y).*detF(-1,y),1,-1);
        Kp34 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi3(-1,y) .* phi4(-1,y).*detF(-1,y),1,-1);
        Kp44 = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* phi4(-1,y) .* phi4(-1,y).*detF(-1,y),1,-1);
        
        Kp_stamp = [Kp11, Kp12, Kp13, Kp14;...
                    0,    Kp22, Kp23, Kp24; ...
                    0,    0,    Kp33, Kp34; ...
                    0,    0,     0,   Kp44];  
        Kp_stamp(2,1) = Kp_stamp(1,2);
        Kp_stamp(3,1) = Kp_stamp(1,3);
        Kp_stamp(4,1) = Kp_stamp(1,4);
        Kp_stamp(3,2) = Kp_stamp(2,3);
        Kp_stamp(4,2) = Kp_stamp(2,4);    
        Kp_stamp(4,3) = Kp_stamp(3,4);
        K_penalty(e(element,:),e(element,:)) = K_penalty(e(element,:),e(element,:)) + Kp_stamp;
        Rp_stamp = [0;0;0;0];
        Rp_stamp(1) = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* Tbar .* phi1(-1,y).*detF(-1,y),1,-1);
        Rp_stamp(4) = my_guass1d(@(y) sqrt(([-1,0]*(Mp(-1,y)*[-1;0]))).* Tbar .* phi4(-1,y).*detF(-1,y),1,-1);       

        R_penalty(e(element,:)) = R_penalty(e(element,:)) + Rp_stamp;
    end    
end


Pstar = PenaltyConstant * max(diag(K));

K = K + Pstar * K_penalty;
R = R + Pstar * R_penalty;

% for element =1:size(e,1)
%     nodes = p(e(element,:),:);
%     if( abs(nodes(1,1).^2 + nodes(1,2).^2 - ri.^2) < 0.00001)
%         K(e(element,1),:) = 0;
%         K(e(element,4),:) = 0;
%         K(e(element,1),e(element,1)) = 1;
%         K(e(element,4),e(element,4)) = 1;
%         R(e(element,1)) = Tbar;
%         R(e(element,4)) = Tbar;
%     end   
% end

Sol = K\R;

%Computing the error 
C1 = @(x,y) qn/k(x,y) * ro + z(x,y)/(2*k(x,y))*ro^2;
C2 = @(x,y) Tbar + z(x,y)/(4*k(x,y)) *ri^2 - C1(x,y)*log(ri);
f = @(x,y) -z(x,y)/(4*k(x,y))* (x.^2 + y.^2) + C1(x,y) * log(sqrt(x.^2 + y.^2)) + C2(x,y);

fpx = @(x,y) (-z(x,y)/(2*k(x,y))* sqrt(x.^2 + y.^2) + C1(x,y)./(sqrt(x.^2 + y.^2))).*x./sqrt(x.^2 + y.^2);
fpy = @(x,y) (-z(x,y)/(2*k(x,y))* sqrt(x.^2 + y.^2) + C1(x,y)./(sqrt(x.^2 + y.^2))).*y./sqrt(x.^2 + y.^2);

num = 0;
denom = 0;

for element =1:size(e,1)
    elem = e(element,:);
    nodes = p(e(element,:),:);

    xx = @(x,y)  nodes(1,1).*phi1(x,y) + nodes(2,1).*phi2(x,y) + nodes(3,1).*phi3(x,y) + nodes(4,1).*phi4(x,y);
    yy = @(x,y)  nodes(1,2).*phi1(x,y) + nodes(2,2).*phi2(x,y) + nodes(3,2).*phi3(x,y) + nodes(4,2).*phi4(x,y);
    F11 = @(x,y) nodes(1,1)*phi1x(x,y) + nodes(2,1)*phi2x(x,y) + nodes(3,1)*phi3x(x,y) + nodes(4,1)*phi4x(x,y); 
    F12 = @(x,y) nodes(1,1)*phi1y(x,y) + nodes(2,1)*phi2y(x,y) + nodes(3,1)*phi3y(x,y) + nodes(4,1)*phi4y(x,y);
    F21 = @(x,y) nodes(1,2)*phi1x(x,y) + nodes(2,2)*phi2x(x,y) + nodes(3,2)*phi3x(x,y) + nodes(4,2)*phi4x(x,y); 
    F22 = @(x,y) nodes(1,2)*phi1y(x,y) + nodes(2,2)*phi2y(x,y) + nodes(3,2)*phi3y(x,y) + nodes(4,2)*phi4y(x,y);
    F = @(x,y) [F11(x,y), F12(x,y);F21(x,y),F22(x,y)];
    detF = @(x,y) F11(x,y).*F22(x,y) - F12(x,y) .* F21(x,y);
    Finv = @(x,y) inv(F(x,y));    
    %M = @(x,y) transpose(Finv(x,y))*Finv(x,y);
    M = @(x,y) Finv(x,y)*transpose(Finv(x,y));
    
    TNErNEthpx = @(x,y) Sol(elem(1))*phi1x(x,y) + Sol(elem(2))*phi2x(x,y) + Sol(elem(3))*phi3x(x,y)+ Sol(elem(4))*phi4x(x,y);
    TNErNEthpy = @(x,y) Sol(elem(1))*phi1y(x,y) + Sol(elem(2))*phi2y(x,y) + Sol(elem(3))*phi3y(x,y)+ Sol(elem(4))*phi4y(x,y);
    numerator = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y).*transpose(M(x,y)*(F(x,y)*[fpx(xx(x,y),yy(x,y));fpy(xx(x,y),yy(x,y))]-[TNErNEthpx(x,y);TNErNEthpy(x,y)]))...
                       *(F(x,y)*[fpx(xx(x,y),yy(x,y));fpy(xx(x,y),yy(x,y))]-[TNErNEthpx(x,y);TNErNEthpy(x,y)]);
    denominator = @(x,y) k(xx(x,y),yy(x,y)).*detF(x,y)*transpose(M(x,y)*[fpx(xx(x,y),yy(x,y));fpy(xx(x,y),yy(x,y))])*[fpx(xx(x,y),yy(x,y));fpy(xx(x,y),yy(x,y))];
    num = num + my_guass(@(x,y) numerator(x,y),-1,1,-1,1);
    denom = denom + my_guass(@(x,y) denominator(x,y),-1,1,-1,1);
end

error = sqrt(num)/sqrt(denom);








