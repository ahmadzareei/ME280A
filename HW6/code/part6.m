function part6


clear all;
close all;

ri = 0.1;
ro = 0.25;

NEr = 20;
NEth = 20;

[Sol, p, e,J] = two_phase(ri,ro,NEr, NEth);

J

% nam1 = sprintf('penalty_sol_2_%d_%d.eps',NEr,NEth);
% nam2 = sprintf('penalty_sol_2_%d_%d_2.eps',NEr,NEth);
% 
% [Sol, p, e,J] = two_phase(ri,ro,NEr, NEth);
% 
% patch('Faces',e,'Vertices',[p,Sol],'FaceVertexCData',Sol,'FaceCol','interp');colorbar; axis square
% set(gca,'FontSize',16);
% xlabel('x','FontSize',16);
% ylabel('y','FontSize',16);
% axis([-0.25,0.25,-0.25,0.25]);
% title(sprintf('T(x,y) - Numerical Solution \n penalty method with NEr=%d, NEth=%d',NEr,NEth),'Fontsize',16);
% 
% print('-depsc',nam1);
% 
% colorbar; 
% axis square
% view(170,24);
% set(gca,'FontSize',16);
% xlabel('x','FontSize',16);
% ylabel('y','FontSize',16);
% zlabel('z','FontSize',16);
% title(sprintf('T(x,y) - Numerical Solution \n penalty method with NEr=%d, NEth=%d',NEr,NEth),'Fontsize',16);
% print('-depsc',nam2);