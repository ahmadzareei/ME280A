function part7


J =[ 5,10, -1.2316e+04;
5,20, -1.3214e+04;
10,10,-1.2332e+04;
10,20,-1.3233e+04;
20,10, -1.2336e+04;
20,20, -1.3239e+04];


plot(J(1:2:end,1),J(1:2:end,3),'-*','LineWidth',2)
set(gca,'FontSize',16);
xlabel('Nr','FontSize',16);
ylabel('Potential Function J','FontSize',16);
grid on
title(sprintf('Potential function J vs NE_r where NEth=%d',J(1,2)),'Fontsize',16);
print('-depsc','J_Nth10');

plot(J(2:2:end,1),J(2:2:end,3),'-*','LineWidth',2)
set(gca,'FontSize',16);
xlabel('Nr','FontSize',16);
ylabel('Potential Function J','FontSize',16);
grid on
title(sprintf('Potential function J vs NE_r where NEth=%d',J(2,2)),'Fontsize',16);
print('-depsc','J_Nth20');

plot(J(1:2:end,1),J(1:2:end,3),'-*','LineWidth',2)
set(gca,'FontSize',16);
xlabel('Nr','FontSize',16);
ylabel('Potential Function J','FontSize',16);
grid on
title(sprintf('Potential function J vs NE_r where NEth=%d',J(1,2)),'Fontsize',16);
ylim([-14000,-12000]);
print('-depsc','J_Nth10_2');

plot(J(2:2:end,1),J(2:2:end,3),'-*','LineWidth',2)
set(gca,'FontSize',16);
xlabel('Nr','FontSize',16);
ylabel('Potential Function J','FontSize',16);
grid on
title(sprintf('Potential function J vs NE_r where NEth=%d',J(2,2)),'Fontsize',16);
ylim([-14000,-12000]);
print('-depsc','J_Nth20_2');



