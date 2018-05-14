function sol = u_Np(w)
% This function gives the derivative of the exact equations at point w 
% exact solution is u_N (x) = f(x) + C x + D
% its derivative is : u_N '(x) = f'(x) + C = g(x) + C
% where f(x) is the general solution and C x + D is the perticular solution
% this function read values of C and D

A1= [2.00,2.5,1.25,.025,4.00,1.75,0.5,0.75,3.25,1.00];
g =@(x,A) -128./(A .*16.75 .*pi) .*cos (16.75 .*pi.*x) + 128./(A.*15.25 .*pi) .*cos(15.25*pi*x);
coeff = load('C.mat');
C = coeff.C;


index = floor(10*w)+1;
if (w ==1)
    index = 10;
end


sol = g(w,A1(index)) + C(index);