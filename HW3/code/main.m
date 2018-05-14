function main 

close all;
clear all;
tole = 1e-12;
error= zeros(3,1);
J= zeros(3,1);
iter = zeros(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;
[u1,u_cg1,error(1),J(1),iter(1)] = fem_linear(N,tole);
x1= 1/N*[0:1:N];

plot(x1,u_cg1,'ko--','LineWidth',2, ...
     'MarkerSize',6);
hold on;
plot(x1,u1,'b*:','LineWidth',2, ...
     'MarkerSize',6);

title('conjugate gradient u(x) and u(x) from backslash','FontSize',16);
xlabel('x');
ylabel('u(x)');
legend('u conjugate gradient','u backslash','Location','northwest');
grid on
set(gca,'FontSize',16); 
print -depsc N100.eps
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
N = 1000;
[u2,u_cg2,error(2),J(2),iter(2)] = fem_linear(N,tole);
x2= 1/N*[0:1:N];

plot(x2,u_cg2,'ko--','LineWidth',2, ...
     'MarkerSize',6);
hold on;
plot(x2,u2,'b*:','LineWidth',2, ...
     'MarkerSize',6);

title('conjugate gradient u(x) and u(x) from backslash','FontSize',16);
xlabel('x');
ylabel('u(x)');
legend('u conjugate gradient','u backslash','Location','northwest');
grid on
set(gca,'FontSize',16); 
print -depsc N1000.eps
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
N = 10000;
[u3,u_cg3,error(3),J(3),iter(3)] = fem_linear(N,tole);
x3= 1/N*[0:1:N];

plot(x3,u_cg3,'ko--','LineWidth',1, ...
     'MarkerSize',3);
hold on;
plot(x3,u3,'b*:','LineWidth',1, ...
     'MarkerSize',3);

title('conjugate gradient u(x) and u(x) from backslash','FontSize',16);
xlabel('x');
ylabel('u(x)');
legend('u conjugate gradient','u backslash','Location','northwest');
grid on
set(gca,'FontSize',16); 
print -depsc N10000.eps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard