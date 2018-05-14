function test1
%this function will test u_N(x) which gives the exact solution at point
% x vs the function [x,y,yp]=exact_solution(N)
% which gives the exact solution by solving the exact equations 
% to give the solution at N+1 nodes 
% solution will y, its derivative is yp
%points are x
N=100;
[x,y,yp]=exact_solution(N);

for i = 1:N+1
    yy(i) = u_N(x(i));
end

max(abs(yy'-y))