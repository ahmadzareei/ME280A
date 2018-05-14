function part5

clear all;
close all;

ri = 0.1;
ro = 0.25;

Const = 1000;
NEr = 5;
NEth = 20;

[Sol, p, e,error] = fem_penalty_method(ri,ro,NEr, NEth,Const);
error


% nam1 = sprintf('penalty_sol_%d.eps',Const);
% nam2 = sprintf('penalty_sol_%d_r.eps',Const);
% 
% [Sol, p, e,error] = fem_penalty_method(ri,ro,NEr, NEth,Const);
% 
% patch('Faces',e,'Vertices',[p,Sol],'FaceVertexCData',Sol,'FaceCol','interp');colorbar; axis square
% set(gca,'FontSize',16);
% xlabel('x','FontSize',16);
% ylabel('y','FontSize',16);
% axis([-0.25,0.25,-0.25,0.25]);
% title(sprintf('T(x,y) - Numerical Solution, penalty method\n  with NEr=%d, NEth=%d and P* = %d Max Kii',NEr,NEth,Const),'Fontsize',16);
% 
% print('-depsc',nam1);
% 
% 
% %exact solution 
% r = linspace(ri,ro, 300);
% Tb = 100;
% qb = -25;
% z = -500;
% k = 0.04;
% C1 = qb/k * ro + z/(2*k)*ro^2;
% C2 = Tb + z/(4*k) *ri^2 - C1*log(ri);
% f = -z/(4*k)*r.^2 + C1 * log(r) + C2;
% figure
% plot([ri:(ro-ri)/NEr:ro]', Sol(1:NEr+1),'r.-','LineWidth',2);hold on; plot(r,f,'b--','LineWidth',2);
% xlabel('r','FontSize',16);
% ylabel('T','FontSize',16);
% set(gca,'FontSize',16);
% title(sprintf('T(r) Numerical Solution, penalty method \n with NEr=%d, NEth=%d with and P* = %d Max Kii ',NEr,NEth,Const),'Fontsize',16);
% legend('Numerical solution with Penalty Method','Exact Solution');
% print ('-depsc',nam2);
