function [a,iter] = my_conjugate_gradient(K,R,a1,tol,mult,el)

iter =1;
r1 = R - mult(K,a1,el);

z1 = r1 ;

l1 = sum(z1.*r1)/sum(z1.*(mult(K,z1,el)));

a2 = a1 + l1*z1;

error = 1;

while (error > tol)
    iter = iter+1;
    r2 = R - mult(K,a2,el);
    th2= - sum(r2.* mult(K,z1,el))/sum(z1.*mult(K,z1,el));
    z2 = r2 + th2*z1;
    l2 = sum (z2.*r2)/sum(z2.*mult(K,z2,el));
    a3 = a2 + l2*z2;
    error = l2 * sum(z2.*mult(K,z2,el))/sum(a2.*mult(K,a2,el));
    a2 = a3;
    z1=z2;
    
end

a = a2;

