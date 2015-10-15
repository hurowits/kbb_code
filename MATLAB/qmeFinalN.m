% clear all

%% init model parameters
eps = logspace(-2,2,100);
ControlParameter = 'driving'
% ControlParameter = 'bath'
Wdot = zeros(length(eps),1);
for q=1:length(eps)
    N = 25; %number of levels
    
    Delta = 1/N; %Level Spacing
    beta=1;
    
    if(strcmp(ControlParameter,'driving'))
        a = eps(1); %transition rates a,b == w_a,w_b of the notes
        b = eps(q);
        gamma = 1; %bath coupling
        gammaPerp = 2;
    elseif (strcmp(ControlParameter,'bath'))
        a = 1;
        b = 5;
        gamma = eps(q);
    end
    lambda = sqrt(a*b);
    
    gammaUp = gamma*exp(-beta*Delta/2);
    gammaDown = gamma*exp(beta*Delta/2);
    
    %% construct superoperator matrix
    W = zeros(ceil(N/2)*(1+N)/2);
    
    W_cl =  -(a+b+gammaUp+gammaDown)* diag(ones(N,1),0) + ...
        (a+gammaDown)* diag(mod(1:N-1,2),1)+(b+gammaDown)*diag(mod(0:N-2,2),1) +...
        (a+gammaUp)* diag(mod(1:N-1,2),-1) + (b+gammaUp)*diag(mod(0:N-2,2),-1);
    W_cl(1,1) = W_cl(1,1)+b+gammaDown;
    W_cl(N,N) = W_cl(N,N)+a+gammaUp;
    
    A_x = 2*eye(N+1) - diag(ones(N,1),1)-diag(ones(N,1),-1);
    A_x = A_x*lambda/2;
    

    W(N+1:2*N-2,1:N) = A_x(1:N,2:N-1).';
    
    
    idx = N+1;
    for k = N-2:-2:1
        
        A = -2*(a+b+gammaPerp + 1i*(N-k+2)/2*Delta)*diag(ones(k,1),0) + ...
            2*a* diag(mod(1:k-1,2),1)  + 2*b*diag(mod(0:k-2,2),1)+...
            2*a* diag(mod(1:k-1,2),-1) + 2*b*diag(mod(0:k-2,2),-1);
        A(1,1) = A(1,1) + b/2;
        A(end) = A(end) + a/2;
        
        W(idx:idx+k-1,idx:idx+k-1) = A;
        A_x_k = A_x(1:k,2:k-1);
        W(idx:idx+k-1, idx+k:idx+k+size(A_x_k,2)-1) = A_x_k;
        W(idx+k:idx+k+size(A_x_k,2)-1, idx:idx+k-1) = A_x_k';
        idx = idx+k;
        
    end
    
    WW = [W 1i*W];
    WW(1:N,1:N) = W_cl;
%     WW(1:N,N+1:2*N-2) = 2*A_x(1:N,2:N-1);
    NormCondition = zeros(1,size(WW,2));
    NormCondition(1,1:N)=ones(1,N);
    %% Solve steady state equation
    
    B = zeros(1,size(WW,1)+1);
    B(end)=1;
    
    [X,R] = linsolve([WW;NormCondition],B');
    %X == [P, Re{\rho}, 0, Im{\rho}]
    
    %% Calculate heating rate
    P = X(1:N).';
    %Heating rate
    Wdot(q) = 2*sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)))/Delta/gamma;
%         Wdot(q) = sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)));
end
%% Generate plots
figure;plot(eps,-Wdot,eps,abs(lambda./(a+eps+gammaPerp+i*Delta)),'LineWidth',2);
grid;%axes('FontSize',15);
axis([0,eps(end),0,1.2])
xlabel('\epsilon','FontSize',27);
ylabel('1/D_B * dW/dt','FontSize',27);
title(['N = ',num2str(N)]);
legend(['dW/dt';'alpha'],'Location','NorthEast')
print(gcf, '-depsc2', 'qme.eps')


% print(gcf, '-depsc2', '/Users/daniel/Desktop/Thesis/code/slrt/figures/dotWvsGammaQuantum.eps')

