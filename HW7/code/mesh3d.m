

Nr = 5;
Nth = 20;
Nc = 20;

Rci = 1;
t = 0.5;
Rs = 4;
Rco = Rci + t;

[theta1, theta2]  = meshgrid(linspace(0,2*pi,Nth+1),linspace(0,2*pi,Nc+1));

x=  (Rs + Rco * cos(theta2)).*cos(theta1);
y= (Rs + Rco * cos(theta2)).*sin(theta1);
z = Rco *sin(theta2);

% mesh(x,y,z);
% axis equal
% 
% figure

% [r,th] = meshgrid(linspace(Rci, Rco, Nr+1),linspace(0,2*pi, Nc+1));
% mesh(Rs+r.*cos(th), r.*sin(th),1.*ones(size(r)))


[r,th] = meshgrid(linspace(Rci, Rco, 2),linspace(0,2*pi, Nc+1));
mesh(Rs+r.*cos(th), r.*sin(th),1.*ones(size(r)));



X = Rs+r.*cos(th);
Z = r.*sin(th);

% X1 = X(:,1); Z1 = Y(:,1); Y1 = zeros(size(X(:,1)));
% X2 = X(:,2); Z2 = Y(:,2); Y2 = zeros(size(X(:,2)));
%polyarea(X2,Z2)-polyarea(X1,Z1)
answ = polyarea(X,Z);
answ(2)-answ(1)
pi*(1.5^2-1)

[r,th] = meshgrid(linspace(Rs-Rco,Rs+Rco, 2),linspace(0,2*pi, Nc+1));
X = r.*cos(th);
Z = r.*sin(th);

%X1 = X(:,1); Z1 = Y(:,1); Y1 = zeros(size(X(:,1)));
%X2 = X(:,2); Z2 = Y(:,2); Y2 = zeros(size(X(:,2)));

%polyarea(X2,Z2)-polyarea(X1,Z1)
answ = polyarea(X,Z);
answ(2)-answ(1)
pi*(5.5^2-2.5^2)
keyboard


% X = (Rs+r.*cos(th))*cos(2*pi/20);
% Y = (Rs+r.*cos(th))*sin(2*pi/20);
% Z = r.*sin(th);
% 
% X3 = X(:,1); Z3 = Y(:,1); Y3 = zeros(size(X(:,1)));
% X4 = X(:,2); Z4 = Y(:,2); Y4 = zeros(size(X(:,2)));
% 
keyboard
