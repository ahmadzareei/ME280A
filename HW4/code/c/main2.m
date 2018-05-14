function main2
N= 16;
[p1,u,error,EI] = fem_linear_EI_a(N);
[p2,u,error,EI] = fem_linear_EI_b(N);

x = 1/N*[0:N];
number1 = zeros(N,1);
number2 = zeros(N,1);
number3 = zeros(N,1);
for i = 1:N
    number1(i) = sum ( p1>=x(i) & p1 <x(i+1));
    number2(i) = sum ( p2>=x(i) & p2 <x(i+1));
end


[p3,u,error,EI] = fem_linear_EI_b2(N);
N1 = 5;
N2 = 11;
points1 = linspace (0,1/3,N1+1);
points2 = linspace (1/3,1,N2+1);
x = [points1(1:end-1), points2]';
for i = 1:N
    number3(i) = sum ( p3>=x(i) & p3 <x(i+1));
end

keyboard