function main_a

% Case a
N=8;
error =1;
while error > 0.05
     N=2*N;
     [u,error] = fem_linear_a(N);
end
M = N/2;

while abs(N-M)>1
    mid = floor((N+M)/2);
    [u,error] = fem_linear_a(mid);
    if(error>0.05)
        M=mid;
    else 
         N=mid;
    end
end
[u,error] = fem_linear_a(N);
fprintf('Linear, N = %d, error = %5.3f\n',N,error);


figure
x = 1/N*[0:N];
plot(x,u,'b+-','LineWidth',2);

hold on

N1 = floor(N/4);
N1 = floor(N1/100) * 100;
x = 1/N1*[0:N1];
[u,error] = fem_linear_a(N1);
plot(x,u,'o-','LineWidth',2);

N2 = N1*2;
N2 = floor(N2/100) * 100;
x = 1/N2*[0:N2];
[u,error] = fem_linear_a(N2);
plot(x,u,'--*','LineWidth',2);

x=linspace(0,1,3000);
u_true = @(x) (10*sin(3*pi*x) + 5).*sin(36*pi*x.^3) ;
plot(x,u_true(x),'k-','LineWidth',2);

set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('u(x)','FontSize',16)
legend(sprintf('N = %d',N),sprintf('N = %d',N1),sprintf('N = %d',N2),'true solution','Location', 'SouthWest');
grid on 
