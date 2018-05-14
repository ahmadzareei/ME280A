function main_a

% Case a
N1 = 1456;
N2 = 369;
N3 = 252;

[p,u,error,EI] = fem_linear_EI_a(N1);
plot(p(1:end-1),EI,'o-','LineWidth',2);
set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('EI','FontSize',16)
legend(sprintf('case(a) - N = %d',N1))
figure

hold on ;

[p,u,error,EI] = fem_linear_EI_b(N2);
plot(p(1:end-1),EI,'--*','LineWidth',2);

[p,u,error,EI] = fem_linear_EI_b2(N3);
plot(p(1:end-1),EI,'b+-','LineWidth',2);

set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('EI','FontSize',16)
legend(sprintf('case(b) - N = %d',N2),sprintf('case (b) - fixed node- N = %d',N3),'Location', 'SouthWest');
keyboard
