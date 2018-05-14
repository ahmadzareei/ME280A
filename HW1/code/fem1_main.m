function fem1_main

k = 2.^[0,1,2,3,4,5];


for i=0:5
    k=2^i;
    N=1;
    error =1;
    while error >0.01
        N=2*N;
        error = fem1_v2(N,k);
    end
    M = N/2;
    while abs(N-M)>1
        mid = floor((N+M)/2);
        error = fem1_v2(mid,k);
        if(error>0.01)
            M=mid;
        else 
            N=mid;
        end
    end
    %mid = floor((N+M)/2);
    error = fem1_v2(N,k);
    fprintf('k = %d, N = %d, error = %5.3f\n',k,N,error);
    
end
