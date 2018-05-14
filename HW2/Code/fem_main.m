function fem_main
% Linear Case
N=1;
error =1;
while error >0.01
     N=2*N;
     [u,error] = fem_linear(N);
end
M = N/2;
while abs(N-M)>1
    mid = floor((N+M)/2);
    [u,error] = fem_linear(mid);
    if(error>0.01)
        M=mid;
    else 
         N=mid;
    end
end
[u,error] = fem_linear(N);
fprintf('Linear, N = %d, error = %5.3f\n',N,error);

% Quadratic Case
N=1;
error =1;
while error >0.01
     N=2*N;
     [u,error] = fem_quadratic(N);
end
M = N/2;
while abs(N-M)>1
    mid = floor((N+M)/2);
    [u,error] = fem_quadratic(mid);
    if(error>0.01)
        M=mid;
    else 
         N=mid;
    end
end
[u,error] = fem_quadratic(N);
fprintf('Quadratic, N = %d, error = %5.3f\n',N,error);

%Cubic Case
N=1;
error =1;
while error >0.01
     N=2*N;
     [u,error] = fem_cubic(N);
end
M = N/2;
while abs(N-M)>1
    mid = floor((N+M)/2);
    [u,error] = fem_cubic(mid);
    if(error>0.01)
        M=mid;
    else 
         N=mid;
    end
end
[u,error] = fem_cubic(N);
fprintf('Cubic, N = %d, error = %5.3f\n',N,error);

