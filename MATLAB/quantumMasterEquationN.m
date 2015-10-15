clear all

%eps = linspace(0,1e4,1e3);
eps = logspace(-2,4,1000);
for q=1:length(eps)
    N = 100; %number of levels
    
    Delta = 1e-8; %Level Spacing
    beta=1;
    
    aa = 5; %transition rates
    bb = 6;
    a = aa^2;
    b = bb^2;
    lambda = aa*bb;
    gamma = eps(q); %bath coupling
    gammaUp = gamma*exp(-beta*Delta/2);
    gammaDown = gamma*exp(beta*Delta/2);
    
    A = zeros(3*N-6);
    
    V_1 = zeros(N,N-2);
    for m=1:N
        for n=1:N-2
            if (m==n || m==n+2)
                V_1(m,n)=-1;
            end
            if(m==n+1)
                V_1(m,n)=2;
            end
            
        end
    end
    V_1 = V_1*lambda;
    V_2 = zeros(N-2,N-4);
    for m=1:N-2
        for n=1:N-4
            if (m==n || m==n+2)
                V_2(m,n)=-1;
            end
            if(m==n+1)
                V_2(m,n)=2;
            end
            
        end
    end
    V_2 = V_2*lambda;
    
    A(N+1:2*N-2,1:N) = V_1.';
    A(2*N-1:end,N+1:2*N-2) = V_2.';
    A = A + A';
    
    
    A_1 = -(a+b+gammaUp+gammaDown)* diag(ones(N,1),0) + ...
        (a+gammaDown)* diag(mod(1:N-1,2),1)+(b+gammaDown)*diag(mod(0:N-2,2),1) +...
        (a+gammaUp)* diag(mod(1:N-1,2),-1) + (b+gammaUp)*diag(mod(0:N-2,2),-1);
    
    A_2 = -2*(a+b+gamma + i*Delta)*diag(ones(N-2,1),0) + ...
        2*a* diag(mod(1:N-3,2),1)  + 2*b*diag(mod(0:N-4,2),1)+...
        2*a* diag(mod(1:N-3,2),-1) + 2*b*diag(mod(0:N-4,2),-1);
    
    
    if (N>3)
        A_3 = -(a+b+gamma+i*2*Delta)*diag(ones(N-4,1),0) + ...
            a* diag(mod(1:N-5,2),1)  + b*diag(mod(0:N-6,2),1) +...
            a* diag(mod(1:N-5,2),-1)/2 + b*diag(mod(0:N-6,2),-1);
        
        A = A + blkdiag(A_1,A_2,A_3);
        A(2*N-2,2*N-2) =A(2*N-2,2*N-2)+ a;
        A = [A;[ones(1,N),zeros(1,2*N-6)]];
        B=[zeros(3*N-6,1);1];
    else
        A = A + blkdiag(A_1,A_2);
        A = [A;[ones(1,N),zeros(1,2*N-5)]];
        B=[zeros(3*N-5,1);1];
    end
    A(1,1) = A(1,1)+b+gammaDown;
    A(N,N) = A(N,N)+a+gammaUp;
    A(N+1,N+1) = A(N+1,N+1) + b;
    
    %Normalization
    
    
    
    %Solution to steady state equation
    [X,R] = linsolve(A,B);
    P = X(1:N).';
    %Heating rate
    Wdot(q) = sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)));
end
figure;plot(eps,-Wdot);
%eps = logspace(-2,2,100); %driving intensity
