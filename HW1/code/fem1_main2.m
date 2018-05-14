function fem1_main2

Nmax = [92,146,338,712,1444,2899];
l=1.0; %Length of the domain
A1=0.2; %constant in the equation
C = (1-l^3/3/A1)/l;
x=linspace(0,1,5000);
err = 0;


for i=0:5
%for i=0
    index = 1;
    k=2^i;
    N=1;
    error =1;
    figure
    uu = -l^2/A1/pi^2*sin(pi*k*x/l) + x.^3/3/A1 + C*x;
    plot(x,uu);
    hold on;
    while N < Nmax(i+1)
        N=2*N;
        [err(index),u] = fem1_v3(N,k);
        plot([0:N]/N,u);
        index = index+1;
    end
    N=Nmax(i+1);
    [err(index),u] = fem1_v3(N,k);
    plot([0:N]/N,u);    
    hold off
    M = [2.^(1:index-2),Nmax(i+1),2.^(index-1)];
    err = [err(1:index-2),err(index), err(index-1)];
    figure
    loglog(M,err);
    caption = sprintf('K= %d',k);
    title(caption);
end