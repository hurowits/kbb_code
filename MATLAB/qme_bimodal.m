function [Wdot, rho] = qme_bimodal(w_a,w_b,N,eps,gamma,Delta)


%% init model parameters
% % eps = logspace(-4,4,100);
% N = 15; %number of levels
ControlParameter = 'driving'

SupMatSize = ceil(N/2)*(1+N)/2;
Wdot = zeros(length(eps),1);
rho = zeros(SupMatSize);
alpha = zeros(length(eps),1);


beta=1;


gammaPerp = gamma;
gammaUp = gamma*exp(-beta*Delta/2);
gammaDown = gamma*exp(beta*Delta/2);
for q=1:length(eps)
    
    
    
    a = w_a*eps(q); %transition rates a,b == w_a,w_b of the notes
    b = w_b*eps(q);
    
    lambda = sqrt(a*b);
    alpha(q) = lambda/((a+b)/2+gammaPerp+1i*Delta);
    
    
    %% construct superoperator matrix
    
    W = zeros(SupMatSize);
    
    %Classical Master equation matrix
    W_cl =  -(a+b+gammaUp+gammaDown)* diag(ones(N,1),0) + ...
        (a+gammaDown)* diag(mod(1:N-1,2),1)+(b+gammaDown)*diag(mod(0:N-2,2),1) +...
        (a+gammaUp)* diag(mod(1:N-1,2),-1) + (b+gammaUp)*diag(mod(0:N-2,2),-1);
    W_cl(1,1) = W_cl(1,1)+b+gammaDown;
    W_cl(N,N) = W_cl(N,N)+a+gammaUp;
    
    % \Lambda: Coupling matrix
    A_x = 2*eye(N+1) - diag(ones(N,1),1)-diag(ones(N,1),-1);
    A_x = A_x*lambda/2;
    
    
    %W_perp
    
    W(N+1:2*N-2,1:N) = A_x(1:N,2:N-1).';
    %     W(1:N,N+1:2*N-2) = A_x(1:N,2:N-1);
    
    idx = N+1;
    for k = N-2:-2:1
        
        A = -2*(a+b+gammaPerp + 1i*(N-k+2)/2*Delta)*diag(ones(k,1),0) + ...
            2*a* diag(mod(1:k-1,2),1)  + 2*b*diag(mod(0:k-2,2),1)+...
            2*a* diag(mod(1:k-1,2),-1) + 2*b*diag(mod(0:k-2,2),-1);
        A(1,1) = A(1,1) + b;
        A(end) = A(end) + a;
        
        W(idx:idx+k-1,idx:idx+k-1) = A;
        A_x_k = A_x(1:k,2:k-1);
        
        W(idx:idx+k-1, idx+k:idx+k+size(A_x_k,2)-1) = A_x_k;
        W(idx+k:idx+k+size(A_x_k,2)-1, idx:idx+k-1) = A_x_k';
        idx = idx+k;
        
    end
    
    %Set up the matrix, so as to take only the real part of rho when needed
    
    W(1:N,1:N) = W_cl;
    WW = [W 1i*W];
    WW(1:N,N+1:2*N-2) = 2*A_x(1:N,2:N-1);
    WW(:,SupMatSize+1:SupMatSize+N)=[];

    
    W_3 = [-a-gammaUp, a+gammaDown, 0, -lambda;...
            a+gammaUp, -(a+b+gammaUp+gammaDown), b+gammaDown, 2*lambda;...
            0, b+gammaUp, -b-gammaDown, -lambda;...
            -lambda/2, lambda, -lambda/2, (((a+b)/2+gamma)^2+Delta^2)/((a+b)/2+gamma);...
            1, 1, 1, 0];
    
    %% Solve steady state equation
    NormCondition = zeros(1,size(WW,2));
    NormCondition(1,1:N)=ones(1,N);
    B = zeros(size(WW,1)+1,1);
    B(end)=1;
    
    % X -> {p,re[rho],im[rho]}
%     [X] = linsolve([WW;NormCondition],B);
    [X] = linsolve(W_3,B);

%     rho(:,q)=X;
    
    %% Calculate heating rate
    P = X(1:N).';

    % Bath normalized heating rate (as in the classical graph):
    Wdot(q) = -sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)))/(Delta*gamma);
    
end


