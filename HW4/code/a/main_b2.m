function main_b2

% Case a
N=8;
error =1;
while error > 0.05
     N=2*N;
     [x,u,error] = fem_linear_b2(N);
end
M = N/2;

while abs(N-M)>1
    mid = floor((N+M)/2);
    [x,u,error] = fem_linear_b2(mid);
    if(error>0.05)
        M=mid;
    else 
         N=mid;
    end
end
[x,u,error] = fem_linear_b2(N);
fprintf('Linear, N = %d, error = %5.3f\n',N,error);


figure
plot(x,u,'b+-','LineWidth',2);

hold on

N1 = floor(N/4);
[x,u,error] = fem_linear_b2(N1);
plot(x,u,'o-','LineWidth',2);

N2 = N1*2;
[x,u,error] = fem_linear_b2(N2);
plot(x,u,'--*','LineWidth',2);

x=linspace(0,1,3000);
A1 = @(x) (x<1/3).*0.2 + (x>=1/3)*2; ; %constant in the equation
u_true = @(x) 512 ./ (4087*pi*A1(x)) .*(268*sin(61*pi*x/4)/61/pi - 244*sin(67*pi*x/4)/67/pi) + ...
              (x<1/3).*(5+ 7680*sqrt(2)/4087/pi).*x + ...
              (x>=1/3).*...
              (2304*sqrt(2) *(8210-768*sqrt(3) + ...
                4087*pi)/(16703569*pi^2) + 3/2 + (1/2 + 768*sqrt(2)/4087/pi).*x);
plot(x,u_true(x),'k-','LineWidth',2);

set(gca,'FontSize',16);
xlabel('x','FontSize',16)
ylabel('u(x)','FontSize',16)
legend(sprintf('N = %d',N),sprintf('N = %d',N1),sprintf('N = %d',N2),'true solution','Location', 'SouthWest');
keyboard