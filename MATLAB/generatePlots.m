% clear all

%eps = linspace(01e4,1e3);
eps = logspace(-3,2,100);
% ControlParameter = 'driving'
ControlParameter = 'bath'
for q=1:length(eps)
    N = 3; %number of levels
    
    Delta = 1/N; %Level Spacing
    beta=1;
    
    if(strcmp(ControlParameter,'driving'))
        a = eps(1); %transition rates a,b == w_a,w_b of the notes
        b = eps(q);
        gamma = 1; %bath coupling
    elseif (strcmp(ControlParameter,'bath'))
        a = 5;
        b = 10;
        gamma = eps(q);
    end
    lambda = sqrt(a*b);
    
    gammaUp = gamma*exp(-beta*Delta/2);
    gammaDown = gamma*exp(beta*Delta/2);
    
    A = zeros(3*N-6);
    V_1 = zeros(N,N-2);
    V_2 = zeros(N-2,N-4);
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
        A_3 = -(a+b+2*gamma+i*4*Delta)*diag(ones(N-4,1),0) + ...
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
    if(mod(N,2)~=0)
        A(1,1) = A(1,1)+b+gammaDown;
        A(N,N) = A(N,N)+a+gammaUp;
        A(N+1,N+1) = A(N+1,N+1) + b;
    else
        A(1,1) = A(1,1)+b+gammaDown;
        A(N,N) = A(N,N)+b+gammaUp;
        A(N+1,N+1) = A(N+1,N+1) + b;
        
    end
    
    
    %Solution to steady state equation
%      [reX,R] = linsolve((A+conj(A))/2,B);
%      [imX,R] = linsolve((A-conj(A))/2/i,B);
%      X = reX+i*imX;
 
 [X,R] = linsolve(A,B);
  P = real(X(1:N).');
    %Heating rate
%     Wdot(q) = sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)))/Delta/gamma;
    Wdot(q) = sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)));
end

figure;plot(eps,-Wdot,'LineWidth',2);
grid;%axes('FontSize',15);
% axis([0,eps(end),0,1.2])
xlabel('\epsilon','FontSize',27);
ylabel('1/D_B * dW/dt','FontSize',27);
title(['N = ',num2str(N)]);
% print(gcf, '-depsc2', '/Users/daniel/Desktop/Thesis/code/slrt/figures/dotWvsGammaQuantum.eps')

