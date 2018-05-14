function main
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%case a
x=linspace(0,1,3000);
u_true = @(x) (10*sin(3*pi*x) + 5).*sin(36*pi*x.^3) ;
plot(x,u_true(x),'k-','LineWidth',2);

hold on 
[p,u,error,EI] = fem_linear_EI_a(16);
plot(p,u,'o--','LineWidth',2);
set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('u(x)','FontSize',16)
legend(sprintf('Adaptive with N = %d',size(EI,1)),'true solution','Location', 'SouthWest');
grid on
print -depsc 'adaptive-case-a.eps'

hold off
figure 
plot(p(1:end-1),EI,'o-','LineWidth',2);
set(gca,'FontSize',16);
xlabel('x_I','FontSize',16)
ylabel('EI(x)','FontSize',16)
legend(sprintf('EI with N = %d',size(EI,1)),'Location', 'NorthWest');
grid on
print -depsc 'EI_case_a.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case (b)
figure
x=linspace(0,1,3000);
A1 = @(x) (x<1/3).*0.2 + (x>=1/3)*2; ; %constant in the equation
u_true = @(x) 512 ./ (4087*pi*A1(x)) .*(268*sin(61*pi*x/4)/61/pi - 244*sin(67*pi*x/4)/67/pi) + ...
              (x<1/3).*(5+ 7680*sqrt(2)/4087/pi).*x + ...
              (x>=1/3).*...
              (2304*sqrt(2) *(8210-768*sqrt(3) + ...
                4087*pi)/(16703569*pi^2) + 3/2 + (1/2 + 768*sqrt(2)/4087/pi).*x);
plot(x,u_true(x),'k-','LineWidth',2);

hold on 
[p1,u1,error1,EI1] = fem_linear_EI_b(16);
plot(p1,u1,'o--','LineWidth',2);
set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('u(x)','FontSize',16)
legend('true solution',sprintf('Adaptive with N = %d',size(EI1,1)),'Location', 'SouthEast');
grid on
print -depsc 'adaptive-case-b.eps'
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%case (b) - 2
figure
x=linspace(0,1,3000);
A1 = @(x) (x<1/3).*0.2 + (x>=1/3)*2; ; %constant in the equation
u_true = @(x) 512 ./ (4087*pi*A1(x)) .*(268*sin(61*pi*x/4)/61/pi - 244*sin(67*pi*x/4)/67/pi) + ...
              (x<1/3).*(5+ 7680*sqrt(2)/4087/pi).*x + ...
              (x>=1/3).*...
              (2304*sqrt(2) *(8210-768*sqrt(3) + ...
                4087*pi)/(16703569*pi^2) + 3/2 + (1/2 + 768*sqrt(2)/4087/pi).*x);
plot(x,u_true(x),'k-','LineWidth',2);

hold on 
[p2,u2,error2,EI2] = fem_linear_EI_b2(16);
plot(p2,u2,'o--','LineWidth',2);
set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('u(x)','FontSize',16)
legend('true solution',sprintf('Adaptive with N = %d',size(EI2,1)),'Location', 'SouthEast');
grid on
print -depsc 'adaptive-case-b2.eps'

hold off

figure 
plot(p1(1:end-2),EI1(1:end-1),'o-',p2(1:end-2),EI2(1:end-1),'k *--','LineWidth',2);
set(gca,'FontSize',16);
xlabel('x_I','FontSize',16)
ylabel('EI(x)','FontSize',16)
legend(sprintf('EI with N = %d',size(EI1,1)),sprintf('EI with N = %d, Fixed Node',size(EI2,1)),'Location', 'NorthWest');
grid on
print -depsc 'EI_case_b.eps'
keyboard
