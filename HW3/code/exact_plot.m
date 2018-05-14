function exact_plot

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

x=[];
y=[];
yp=[];
for i=1:10
    xx = linspace((i-1)*dx,i*dx,20);
    yy = f(xx,A1(i)) + C(i) *xx + D(i);
    yyp = A1(i)*(g(xx,A1(i)) + C(i));
    x=[x,xx];
    y=[y,yy];
    yp=[yp,yyp];
end

subplot(2,1,1)
plot(x,y)
title('Continuity of exact solution');
xlabel('x');ylabel('u(x)');
subplot(2,1,2)
plot(x,yp);
title('Continuity of A1 du/dx ');
xlabel('x');ylabel('A_1 du(x)/dx');
keyboard
