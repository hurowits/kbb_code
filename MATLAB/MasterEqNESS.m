
% runs vs eps  (eps/gamma)
% runs for various Delta (i.e. change gamma/Delta)
% look for circumstances in which the "off diagonals" are large.


clear all
N =25; %odd
Delta_vec = 1;
gamma_vec = 1;%logspace(-2,5,100);
beta=4;
wa_vec = 1;
wb_vec = logspace(-2,5,100);
% Estrength_vec = [1:500];
Estrength_vec =1;% logspace(-2,2,100);




% Mpoints=length(Estrength_vec);
Mpoints=max([length(Estrength_vec),length(wa_vec),length(wb_vec),length(Delta_vec),length(gamma_vec)]);
alpha = zeros(Mpoints,1);
SupMatSize = ceil(N/2)*(1+N)-N;
Wdot = zeros(Mpoints,1);
rho = zeros(SupMatSize,Mpoints);
complex_rho = zeros((SupMatSize-N)/2,Mpoints);
W = zeros(SupMatSize);
qWdot = zeros(1,Mpoints);
cWdot = zeros(1,Mpoints);
rho_cl = zeros(N,Mpoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP

for q=1:Mpoints
    
    %% init model parameters
    %% INPUT(w_a,w_b,N,Estrength,gamma,Delta)
    waitbar(q/Mpoints);
    
    Delta = Delta_vec;
    gamma = gamma_vec;
    w_a = wa_vec;
    w_b = wb_vec(q);
%     Estrength = Estrength_vec(q);
    Estrength = Estrength_vec;
    a = w_a*Estrength;             %transition rates a,b == w_a,w_b of the notes
    b = w_b*Estrength;
    lambda = sqrt(a*b);
    gammaPerp = gamma;
    gammaUp = gamma*exp(-beta*Delta/2);
    gammaDown = gamma*exp(beta*Delta/2);
    
    
    alpha(q) = lambda/((a+b)/2+gammaPerp+1i*Delta);
    S(q) = lambda*((a+b)/2+gammaPerp)/(((a+b)/2+gammaPerp)^2+Delta^2);
    
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
    W(1:N,1:N) = W_cl;
    
    W(1:N,N+1:2*N-2) = 2*A_x(1:N,2:N-1);
    W(N+1:2*N-2,1:N) = A_x(1:N,2:N-1).';
    
    
    idx = N+1;
    for k = N-2:-2:1
        
        
        A = -((a+b)+gammaPerp + 1i*(N-k)/2*Delta)*diag(ones(k,1),0) + ...
            1*a* diag(mod(1:k-1,2),1)  + 1*b*diag(mod(0:k-2,2),1)+...
            1*a* diag(mod(1:k-1,2),-1) + 1*b*diag(mod(0:k-2,2),-1);
        A(1,1) = A(1,1) + b/2;
        A(end) = A(end) + a/2;
        
        W(idx:idx+k-1,idx:idx+k-1) = real(A);
        W(idx:idx+k-1, idx+k:idx+2*k-1) = -imag(A);
        
        W(idx+k:idx+2*k-1,idx:idx+k-1) = imag(A);
        W(idx+k:idx+2*k-1,idx+k:idx+2*k-1) = real(A);
        
        
        A_x_k = A_x(1:k,2:k-1);
        if(~isempty(A_x_k))
            
            W(idx:idx+k-1, idx+2*k:idx+3*k-3) = A_x_k;
            W(idx+k:idx+2*k-1, idx+3*k-2:idx+4*k-5) = A_x_k;
            
            W( idx+2*k:idx+3*k-3,idx:idx+k-1) = A_x_k';
            W(idx+3*k-2:idx+4*k-5, idx+k:idx+2*k-1) = A_x_k';
            
        end
        idx = idx+2*k;
        
    end
    
    %
    %     W_3 = [-a-gammaUp, a+gammaDown, 0, -lambda;...
    %             a+gammaUp, -(a+b+gammaUp+gammaDown), b+gammaDown, 2*lambda;...
    %             0, b+gammaUp, -b-gammaDown, -lambda;...
    %             -lambda/2, lambda, -lambda/2, (((a+b)/2+gamma)^2+Delta^2)/((a+b)/2+gamma);...
    %             1, 1, 1, 0];
    
    %% Solve steady state equation
    NormCondition = zeros(1,size(W,2));
    NormCondition(1,1:N)=ones(1,N);
    B = zeros(size(W,1)+1,1);
    B(end)=1;
    
    % X -> {p,re[rho],im[rho]}
    %     [X] = linsolve([WW;NormCondition],B);
    [X] = linsolve([W;NormCondition],B);
    
    rho(:,q)=X;
    
    l=N+1;
    for k=N-2:-2:1
     complex_rho(l:l+k-1,q)=rho(l:l+k-1,q) +1i*rho(l+k:l+2*k-1,q);
     l=l+2*k;    

    end
%     complex_rho(1:N,q) = X(1:N);
    %% Calculate heating rate
    P = X(1:N).';
    
    % Bath normalized heating rate (as in the classical graph):
    qWdot(q) = -sum(diff(P).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)))/(Delta*gamma);
    
    %stochastic FGR
    Pcl = linsolve([W_cl;ones(1,N)],[zeros(N,1);1]).';
    rho_cl(:,q)=Pcl;
    cWdot(q) = -sum(diff(Pcl).*(a*mod(1:N-1,2)+b*mod(0:N-2,2)))/(Delta*gamma);
   
% Classical Analytical solution for bimodal case
    
    AA = (a+gamma*exp(-beta*Delta/2))/(a+gamma*exp(beta*Delta/2));
    BB = (b+gamma*exp(-beta*Delta/2))/(b+gamma*exp(beta*Delta/2));
    p0 = sum((AA.^ceil((0:N-1)/2)).*(BB.^floor((0:N-1)/2)))^(-1);
    P_cl_an(:,q) = (AA.^ceil((0:N-1)/2)).*(BB.^floor((0:N-1)/2))*p0;
    caWdot(q) = -sum(diff(P_cl_an(:,q).').*(a*mod(1:N-1,2)+b*mod(0:N-2,2)))/(Delta*gamma);
    
end

% END OF MAIN LOOP
% OUTPUT(qWdot,cWdot,rho,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FIGURES

%Fig1 - $rho_cl$ and $rho_qm$ profiles for the last run

%Fig2a - heating rate vs run
%fig2b - error = -log(pN/p0)  vs run
%fig2c - max(rho offdiagonal)/p0  vs run


%%%%%%%%%%%%%
figure(1);
% plot(rho(1:N,1),'-or')
hold on
% 
%  plot(abs(complex_rho(N+1:end,end)),'-.r') %quantum numeric
plot(rho(1:N,end),'r') %quantum numeric
plot(P_cl_an(:,end),'b'); %classical analytic
xlabel('n','FontSize',24);
ylabel('p_n','FontSize',24);
hold off
% title('Quantum(red) vs. Stochastic (blue) Probabilities')

%Reshape rho as matrix
rho_mat = zeros(N);
k=0;
ind = 1;

for q=N:-2:1
    
rho_mat = rho_mat + diag(complex_rho(ind:ind+q-1,end),k);
k = k + 2;
ind = ind+q;
end

figure(2); imagesc(abs(rho_mat)+abs(rho_mat).'+diag(P));axis image; 
title('Density Matrix');

%%%%%%%%%%

% figure;
% axes('FontSize',14)
% 
% loglog(wa_vec/wb_vec,cWdot,'b','LineWidth',1);
% hold on
% loglog(wa_vec/wb_vec,qWdot,'r','LineWidth',1);
% % loglog(Estrength_vec,abs(alpha),'g');
% 
% plot(wa_vec/wb_vec,log((rho(1,:)./(rho(N,:)))),'g','LineWidth',1);
% loglog(wa_vec/wb_vec,max(abs(complex_rho),[],1)./rho(1,:),'m','LineWidth',1);
% % plot(Estrength_vec,-log(rho_cl(N,:)./rho_cl(1,:)),'k','LineWidth',1);
% hold off;
% 
% grid;
% %axis([0,Estrength(end),0,2])
% 
% % title(['N = ',num2str(N)]);
% % xlabel('\epsilon / \gamma','FontSize',24);
% xlabel('w_a / w_b','FontSize',24);
% % ylabel('dW/dt','FontSize',24);
% legend(['classical  ';'quantum    ';'-logP_N/P_0';'max coh/p_0'],'Location','NorthEast');
% % print(gcf, '-depsc2', 'fig2_wa_N99.eps');
% 



figure;
axes('FontSize',14)

% loglog(wb_vec/Delta,cWdot,'r','LineWidth',1);

loglog(wb_vec/Delta,qWdot,'b','LineWidth',1);
hold on
loglog(wb_vec/Delta,caWdot,'r','LineWidth',1);
% loglog(Estrength_vec,abs(alpha),'g');
% 
% plot(wb_vec/Delta,log10((rho(1,:)./(rho(N,:)))),'c','LineWidth',1);
 if(wa_vec & wb_vec)
 plot(wb_vec/Delta,max(abs(complex_rho),[],1)./rho(1,:),'k','LineWidth',1);
 end
plot(wb_vec/Delta,abs(alpha),'m','LineWidth',1);
%plot(Estrength_vec,-log(rho_cl(N,:)./rho_cl(1,:)),'k','LineWidth',1);
hold off;

grid;
%axis([0,Estrength(end),0,2])

% title(['N = ',num2str(N)]);
% xlabel('\epsilon / \gamma','FontSize',24);
xlabel('w_b','FontSize',24);
ylabel('dW/dt','FontSize',24);
legend(['stochastic ';...
    'quantum    ';...
%     'analytical ';...
%     '-logP_N/P_0';...
    'max coh/p_0';...
    '|alpha|    '],'Location','NorthEast');
% print(gcf, '-depsc2', 'run1.eps');

c=cosh(beta*Delta/2);
s=sinh(-beta*Delta/2);
w=a
mu = ((w+gamma)^2+Delta^2)/(w*(w+gamma));
s^2/(mu*(w/gamma+c^2-0.5*(w/gamma+exp(-beta*Delta/2))*(w/gamma+exp(beta*Delta/2))-3*w/2/gamma*(c+w/gamma)));
