function [x,y,yp]=exact_solution(N)
% This function gives the exact equations 
% N: number of nodes 
% x: N+1 nodes at [0,1] 
% y: solution at points x
% yp: derivative of solution at points x

% A1 constant that are given 
A1= [2.00,2.5,1.25,.25,4.00,1.75,0.5,0.75,3.25,1.00];
X=[0.1:0.1:0.9];
dx = 0.1;
f =@(x,A) -128/(A *16.75^2 *pi^2) *sin (16.75 *pi*x) + 128/(A * 15.25^2 *pi^2) *sin(15.25*pi*x);
g =@(x,A) -128/(A *16.75 *pi) *cos (16.75 *pi*x) + 128/(A*15.25 *pi) *cos(15.25*pi*x);


C = zeros(10,1);
D=zeros(10,1);

C(10)  = 1/A1(10) - g(1,A1(10));

for i =9:-1:1
    C(i) = A1(i+1)/A1(i)*(C(i+1) + g(i*dx,A1(i+1))) - g(i*dx,A1(i));
end

D(1) = -f(0,A1(1));
for i=1:9
    D(i+1) = (D(i) + C(i) * dx*i + f(dx*(i),A1(i))) - ( C(i+1)*dx*i+ f(dx*(i),A1(i+1)));
end

x=0:1/N:1;
y=zeros(N+1,1);
yp=zeros(N+1,1);
for i =1:N+1
    e = floor(x(i)*10)+1;
    if(i == N+1)
        e = 10;
    end
    y(i) = f(x(i),A1(e)) + C(e) *x(i) + D(e);
    yp(i) = A1(e)*(g(x(i),A1(e)) + C(e));       
end
keyboard