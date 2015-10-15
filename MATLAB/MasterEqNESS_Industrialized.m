%% Initialization
clear all
flagSave = 0;
N =25; %odd
Delta_vec = 1;
gamma_vec = 0.1;
beta_vec=0.1;

load randInd;
% s_vec = 10.^(-linspace(6,0.001,100));
% s_vec = linspace(eps,1-eps,100);
% 
% s_vec(end)=10^(-eps);%(logspace(-2,-eps,100);%[1-eps,0.5,0.1,0.01,0.001,0.0001,1e-7,1e-9]
% s_vec = 1e-15;
% s_vec = linspace(eps,1,100); %here we use sigma instead of s, but leave variable name unchanged
% % s_vec = 100;
% % s_vec = logspace(1,10)
% s_vec=linspace(1e-1,1-1e-7,100);
% s_vec = 10.^(-linspace(10,1,50));
% s_vec = [s_vec linspace(0.1,1-1e-8,50)];

% s_vec=linspace(1e-4,0.01,100);
% s_vec=logspace(-4,-2,100)
% s_vec=[s_vec,10.^(linspace(-2,2,300))];
s_vec=0.2
randX = genRand(s_vec,randInd,N,'sigma');

% 
% s_vec=logspace(-6,-eps,100);
% 
% % s_vec=0.5;
% s_vec = 0.99;
% randX = genRand(s_vec,randInd,N,'s');
Estrength_vec = logspace(-2,4,100);
Estrength_vec = linspace(0,1e4,1000);
Estrength_vec = 0.1;
% 
% Estrength_vec = 1e4;

Mpoints=max([length(Estrength_vec),length(Delta_vec),length(beta_vec),length(s_vec)]);
alpha = zeros(Mpoints,1);
SupMatSize = ceil(N/2)*(1+N)-N;
Wdot = zeros(Mpoints,1);
rho = zeros(SupMatSize,Mpoints);
complex_rho = zeros((SupMatSize-N)/2,Mpoints);
rho_mat = zeros(N,N,Mpoints);
W = zeros(SupMatSize);
qWdot = zeros(length(s_vec),length(Estrength_vec));
cWdot = zeros(length(s_vec),length(Estrength_vec));
cQdot = zeros(length(s_vec),length(Estrength_vec));
qQdot = zeros(length(s_vec),length(Estrength_vec));
cPrtNum = zeros(length(s_vec),length(Estrength_vec));
qPrtNum = zeros(length(s_vec),length(Estrength_vec));
rho_cl = zeros(N,Mpoints);
flag=zeros(Mpoints,1);
T_eff_q = zeros(length(s_vec),length(Estrength_vec));
T_eff_c = zeros(length(s_vec),length(Estrength_vec));
T_eff_ss= zeros(length(s_vec),length(Estrength_vec));
T_mix= zeros(length(s_vec),length(Estrength_vec));
T_mix_rmse=zeros(length(s_vec),length(Estrength_vec));
P_cl = zeros(N,length(s_vec),length(Estrength_vec));
P_r = zeros(N,length(s_vec),length(Estrength_vec));
E_r = zeros(N,length(s_vec),length(Estrength_vec));
dE_r = zeros(N,length(s_vec),length(Estrength_vec));
dispersion1 = zeros(length(s_vec),length(Estrength_vec));
dispersion2 = zeros(length(s_vec),length(Estrength_vec));
prtNum_r = zeros(length(s_vec),length(Estrength_vec));
prtNum_W = zeros(length(s_vec),length(Estrength_vec));
quantumCorrection=zeros(N,length(s_vec),length(Estrength_vec));
P_n = zeros(N,length(Estrength_vec));
P_n_q=zeros(N,length(s_vec),length(Estrength_vec));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
for iGamma = 1:length(gamma_vec)
    for iDelta = 1:length(Delta_vec)
        for iBeta = 1:length(beta_vec)
            for iS=1:length(s_vec)
                wb=waitbar(iS/100);
                x = randX(iS,:);
                x = x/mean(x); %normalize to unitary mean
                %                  a = 1;
                %                  b = s_vec(iS);%1e-5;
                %                  x = a*mod(1:N-1,2)+b*mod(0:N-2,2);
                %         x=ones(1,24);
                beta=beta_vec(iBeta);
                
                for q=1:length(Estrength_vec)
                    
                    Delta = Delta_vec(iDelta);
                    gamma = gamma_vec(iGamma);
                    
                    gammaPerp = gamma;
                    gammaUp = gamma*exp(-beta*Delta/2)/cosh(-beta*Delta/2);
                    gammaDown = gamma*exp(beta*Delta/2)/cosh(-beta*Delta/2);
                    
                    E_strength = Estrength_vec(q);
                    
                    w = E_strength*sqrt(x);%matrix element. transition rate = matrixelement^2
                    
                    T_eff_ss(iS,q) = (mean(gamma./(w.^2+gamma)))^-1 / beta;
                    
                    % w =sqrt(E_strength*( 1*mod(1:N,2)+10*mod(0:N-1,2) ));
                    
                    %Classical Master equation matrix
                    
                    W_cl = -diag([w(1:N-1).^2+gammaUp,0],0) - diag([0,w(1:N-1).^2+gammaDown],0) +...
                        diag(w(1:N-1).^2+gammaDown,1) + diag(w(1:N-1).^2+gammaUp,-1);
                    
                    
                    %W_perp
                    W(1:N,1:N) = W_cl;
                    [ii,jj]=meshgrid(1:N-2,1:N);
                    Lambda = zeros(N,N-2);
                    Lambda(ii==jj) = -0.5 *w(1:N-2).*w(2:N-1);
                    Lambda(ii==jj-1) = w(1:N-2).*w(2:N-1);
                    Lambda(ii==jj-2) = -0.5 *w(1:N-2).*w(2:N-1);
                    W(1:N,N+1:2*N-2) = 2*Lambda;
                    W(N+1:2*N-2,1:N) = Lambda';
                    
%                     Lambda0(:,:,q) = Lambda;
                    
                    idx = N+1;
                    for k = N-2:-2:1
                        
                        diagIdx = N-k;
                        
                        A = diag(w(1:k-1).*w(diagIdx+1:N-1),-1)+diag(w(1:k-1).*w(diagIdx+1:N-1),1)-...
                            0.5*(diag([0,w(1:k-1).^2],0)+diag(w(1:k).^2,0)+diag(w(diagIdx:N-1).^2,0)+...
                            diag([w(diagIdx+1:N-1).^2,0],0))-...
                            diag(ones(1,k),0)*(gammaPerp+1i*(N-k)/2*Delta);
                        
                        
                        W(idx:idx+k-1,idx:idx+k-1) = real(A);
                        W(idx:idx+k-1, idx+k:idx+2*k-1) = -imag(A);
                        
                        W(idx+k:idx+2*k-1,idx:idx+k-1) = imag(A);
                        W(idx+k:idx+2*k-1,idx+k:idx+2*k-1) = real(A);
                        
                        [ii,jj]=meshgrid(1:k-2,1:k);
                        Lambda = zeros(k,k-2);
                        Lambda(ii==jj) = -0.5 *w(diagIdx+1:N-2).*w(diagIdx+2:N-1);
                        Lambda(ii==jj-1) = w(1:k-2).*w(diagIdx+2:N-1);
                        Lambda(ii==jj-2) = -0.5*w(1:k-2).*w(2:k-1);
                        
                        
                        if(~isempty(Lambda))
                            
                            W(idx:idx+k-1, idx+2*k:idx+3*k-3) = Lambda;
                            W(idx+k:idx+2*k-1, idx+3*k-2:idx+4*k-5) = Lambda;
                            
                            W( idx+2*k:idx+3*k-3,idx:idx+k-1) = Lambda';
                            W(idx+3*k-2:idx+4*k-5, idx+k:idx+2*k-1) = Lambda';
                            
                        end
                        idx = idx+2*k;
                        %
                    end
                    %
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
                    [X,R] = linsolve([W;NormCondition],B);
                    
                    rho(:,q)=X;
                    
%                     quantumCorrection(:,q,iS) = Lambda0(:,:,q)*rho(26:48,q);
                    
                    
                    fst_l=N+1;
                    sl_l = N+1;
                    for k=N-2:-2:1
                        complex_rho(sl_l:sl_l+k-1,q)=rho(fst_l:fst_l+k-1,q) +1i*rho(fst_l+k:fst_l+2*k-1,q);
                        fst_l=fst_l+2*k;
                        sl_l = sl_l + k;
                        
                    end
                    P = X(1:N).';
                    P_n(:,q) = P; %P in energy basis
                    P_n_q(:,iS,q)=P;
                    %Reshape rho as matrix
                    rho_mat(:,:,q) = rho_vec2mat (complex_rho,P,N);
                    
                    %P in diagonal basis of rho
                    %                     [v(:,:,q),D(:,:,q)] = eig(rho_mat(:,:,q));
                    %                     P_d(:,q) = diag(D(:,:,q));
                    %                     E_d(:,q) = abs(v(:,:,q)').^2*(1:N)';
                    
                    
                    %P in basis that diagonalizes W
                    W_driving = diag(w,1)+diag(w,-1);
                    [r,D_r] = eig(W_driving);
                    P_r(:,iS,q) = diag(inv(r)*rho_mat(:,:,q)*r);
                    E_r(:,iS,q) = ((1:N)*abs(r).^2)';
                    dE_r(:,iS,q) = sqrt(abs(r').^2*((1:N)'-E_r(:,iS,q)).^2);
                    dispersion1(iS,q) = std(E_r(:,iS,q));
                    dispersion2(iS,q) = max(E_r(:,iS,q))-min(E_r(:,iS,q));
                    prtNum_r(iS,q) = sum(P_r(:,iS,q).^2)^-1;
                    prtNum_W(iS,q) = mean(sum(abs(r').^4,1).^-1);

                    
                    % P in position basis
                    %                     z = exp(1i*2*pi*(0:N-1)/N);
                    %                     for iF = 1:N
                    %                         F(iF,:) = z.^(iF-1);
                    %                     end
                    %                     F = F/sqrt(N);
                    %                     rho_xx = inv(F)*rho_mat(:,:,q)*F;
                    %
                    %
                    %                     P_x(:,q) = diag(rho_xx);
                    %
                    %Effective Temperature
                    T_n = 1./log(P(1:N-1)./P(2:N)) ;
                    T_eff_q(iS,q) = 1./mean(1./T_n);
                   % dT_eff_q(q) = 1./std(1./T_n);
                    
                    [fittedP,gof] = fit(real(E_r(:,iS,q)),real(P_r(:,iS,q)),'exp1');
                    
                    C = coeffValues(fittedP);
                    T_mix(iS,q) = -1/C(2);
                    T_mix_rmse(iS,q) = gof.rmse;
%                    T_n_r(:,q) = 1./log(P_r(1:N-1,iS,q)./P_r(2:N,iS,q)).* (E_r(2:N,iS,q)-E_r(1:N-1,iS,q));
%                    T_mix2(iS,q) = 1./mean(1./T_n_r(:,q));

                    %% Calculate heating rate
                    %                     figure(1);plot(abs(v(:,25,100)));
                    qPrtNum(iS,q) = sum(P.^2)^-1;
                    %                     S_tr(q,iS) = sqrt(sum(abs(complex_rho(:)).^2));
                    %S_tr(iS,q) = sqrt(sum(X(N+1:end).^2));
                    % Bath normalized heating rate (as in the classical graph):
                    
                    qWdot(iS,q) = -sum(diff(P).*w(1:N-1).^2)/(Delta*gamma);
                    qQdot(iS,q) = (gammaDown*( 1-P(1) ) - gammaUp*( 1 - P(N) ))/(Delta*gamma);
                    
                    if (sum(P<0))
                        flag(q)=1;
                    end
                    
                    %%stochastic FGR
                    Pcl = linsolve([W_cl;ones(1,N)],[zeros(N,1);1]).';
                    cPrtNum(iS,q) = sum(Pcl.^2)^-1;
                    rho_cl(:,q)=Pcl;
                    %         cWdot(q,qq) =-sum(diff(Pcl).*w(1:N-1).^2)/(Delta*gamma);
                    cWdot(iS,q) =-sum(diff(Pcl).*w(1:N-1).^2)/(Delta*gamma);
                    
                    cQdot(iS,q) = gammaDown*( 1-Pcl(1) ) - gammaUp*( 1 - Pcl(N) );
                    P_cl(:,iS,q) = Pcl.';
                    %Effective Temperature
                    T_n =1./ log(Pcl(1:N-1)./Pcl(2:N));
                    T_eff_c(iS,q) = 1./mean(1./T_n);
                    %dT_eff_c(q) = 1./std(1./T_n);
                    %                 a=w(1)^2 ;
                    %                 AA = ((a+gamma*exp(-beta*Delta/2))/(a+gamma*exp(beta*Delta/2))).^(0:N-1);
                    %                 p0 = sum(AA)^(-1);
                    %                 P_cl_an(:,q) = AA*p0;
                    %                 caWdot(q) = -sum(diff(P_cl_an(:,q).')*a)/(Delta^2*gamma);
                    %
                    %
                end

            end
        end
    end
end

%% Save and reload complete data sets
if (flagSave)
    save Temperature5 T_eff_q T_eff_c T_eff_ss T_mix dispersion1 dispersion2 s_vec Estrength_vec randX gamma Delta beta N
    save Probability P_cl P_n_q P_r
    save Wdot5        cWdot   qWdot            s_vec Estrength_vec randX gamma Delta beta N
    save PN5          cPrtNum qPrtNum          s_vec Estrength_vec randX gamma Delta beta N
    
    load Temperature5
    load Wdot5
    load PN5
end
%% Figs. Probabilities
%                 d_sorted = sortrows([E_d,P_d]);
iS = 1
iE = 1

r_sorted = sortrows([squeeze(E_r(:,iS,iE)),squeeze(P_r(:,iS,iE))],1);
figure;
axes('FontSize',14)
plot(1:N,P_cl(:,iS,iE),'-b','LineWidth',3);
hold on
plot(1:N,P_n_q(:,iS,iE),'-dr','MarkerSize',5,'LineWidth',2);
plot(r_sorted(:,1),r_sorted(:,2),'-sg','MarkerSize',5,'LineWidth',2);
% plot(E_r(:,iS,iE),P_r(:,iS,iE),'.-g','MarkerSize',5,'LineWidth',2);
plot(1:N, exp(-beta*(1:N))./sum(exp(-beta*(1:N))),'--k','LineWidth',3);
legend(['Stochastic     ';...
        'E Basis        ';...
        'V Basis        ';...
        'exp(-\beta E_n)' ]);
% axis([1,N,0.01,0.08]);
xlabel('E_n','FontSize',24);
ylabel('p','FontSize',24);

%% Plots of Heating Rates
figure;iS=1;
x = randX(iS,:);
x = x./mean(x);
axes('FontSize',14)
E=Estrength_vec;%*sqrt(mean(x));
loglog(E.^2,cWdot(iS,:),'b','LineWidth',3);

hold on
%     loglog(E,caWdot,'b','LineWidth',1);
loglog(E.^2,qWdot(iS,:),'--r','LineWidth',3);


loglog(E.^2,E.^2*mean(x)*beta/gamma,'g--','LineWidth',3);    %LRT
plot(E.^2,ones(length(E),1)*beta,'g--','LineWidth',3); %bath limited saturation

plot(ones(length(E),1)*((gamma)/mean(x)),linspace(1e-8,100,length(E)),'--k','LineWidth',3); %cross over to SLRT
plot(ones(length(E),1)*((gamma)/min(x)),linspace(1e-8,100,length(E)),'--k','LineWidth',3); %crossover to saturation

% plot(ones(length(E),1)*(mean(sqrt(x))-sqrt(mean(sqrt(x))^2*(1+4*gamma^2)-4*geomean(sqrt(x))^2*(gamma^2+Delta^2)))/2/(geomean(sqrt(x))^2-mean(sqrt(x))^2),linspace(1e-8,100,length(E)),'g--','LineWidth',1); %crossover to saturation


text(((gamma)/mean(x))/25,1e-2,'LRT','FontSize',14,'FontName','Times New Roman','FontWeight','bold' );

% text(((gamma)/min(x))/400,1e-2,'SLRT','FontSize',14,'FontName','Times New Roman','FontWeight','bold' );
   
text(((gamma)/median(x))*10^2,1e-2,'Saturation','FontSize',14,'FontName','Times New Roman','FontWeight','bold');

hold off;
%         grid;
xlabel('\epsilon^2','FontSize',24);
ylabel('Rate of heating','FontSize',24,'FontName','Times New Roman');

L=legend(['Stochastic  ';...
    'Quantum     '],'Location','SouthEast');

axis([Estrength_vec(1).^2, Estrength_vec(end)^2,1e-3,1]);
%% Images of Heating Rates
figure;
E=Estrength_vec;
ax=axes('FontSize',14);
h = imagesc(cWdot);%
set(h,'XData',[log10(E(1)^2) log10(E(end)^2)])
set(h,'YData',[log10(s_vec(1)) log10(s_vec(end))]);
set(ax,'YDir','normal')
axis tight; 
% colorbar;
xlabel('log \epsilon^2','FontSize',24);
ylabel('log s','FontSize',24);
% set(gca,'YTick',[1 35 69])
% set(gca,'XTick',[34 67 length(E)])
% set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
% set(gca,'YTick',[1 35 69])
% set(gca,'XTick',[34 67 length(E)])
% set(gca,'YTickLabel',[num2str(floor(log10(s_vec(9)))),'|',num2str(ceil(log10(s_vec(35)))),'|',num2str(ceil(log10(s_vec(69))))]);
% set(gca,'XTickLabel',[num2str(floor(log10(E(34)^2))),'|',num2str(floor(log10(E(67)^2))),'|',num2str(floor(log10(E(end)^2)))]);


figure;
ax=axes('FontSize',14);
h = imagesc(qWdot/(Delta));
%colorbar
set(h,'XData',[log10(E(1)^2) log10(E(end)^2)])
set(h,'YData',[log10(s_vec(1)) log10(s_vec(end))]);
set(ax,'YDir','normal')
axis tight;
xlabel('log \epsilon^2','FontSize',24);
ylabel('log s','FontSize',24);
% set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
% set(gca,'YTick',[1 35 69])
% set(gca,'XTick',[34 67 length(E)])
% set(gca,'YTickLabel',[num2str(floor(log10(s_vec(9)))),'|',num2str(ceil(log10(s_vec(35)))),'|',num2str(ceil(log10(s_vec(69))))]);
% set(gca,'XTickLabel',[num2str(floor(log10(E(34)^2))),'|',num2str(floor(log10(E(67)^2))),'|',num2str(floor(log10(E(end)^2)))]);



%% Plots of Partition Numbers
figure;
N=25;
iS=1;
PN = (1-exp(-2*beta*Delta))/(1-exp(-2*beta*Delta*N))*(1-exp(-beta*Delta*N))^2/(1-exp(-beta*Delta))^2;
axes('FontSize',14)
semilogx(E.^2,cPrtNum(iS,:),'b','LineWidth',3);
hold on
semilogx(E.^2,qPrtNum(iS,:),'r','LineWidth',3);
semilogx(E.^2,PN*ones(1,length(E)),'--g','LineWidth',3);
hold off;
xlabel('\epsilon^2','FontSize',24);
ylabel('PN','FontSize',24,'interpreter','latex');
legend(['Stochastic      ';...
        'Quantum         ';...
        'PN[e^{-\beta E}]'],'Location','NorthWest');
axis([0,E(end)^2,15,30]);

figure;
PN = (1-exp(-2*beta*Delta))/(1-exp(-2*beta*Delta*N))*(1-exp(-beta*Delta*N))^2/(1-exp(-beta*Delta))^2;
axes('FontSize',14)
semilogx(s_vec,cPrtNum(:,end),'b','LineWidth',3);
hold on
semilogx(s_vec,qPrtNum(:,end),'r','LineWidth',3);
semilogx(s_vec,PN*ones(1,length(E)),'--g','LineWidth',3);
hold off;
xlabel('s','FontSize',24);
ylabel('PN_{infty}','FontSize',24,'interpreter','tex');
legend(['Stochastic      ';...
        'Quantum         ';...
        'PN[e^{-\beta E}]'],'Location','NorthWest');
axis([0,s_vec(end),15,30]);


%% Images of Partition Numbers
figure;
axes('FontSize',14)
imagesc(qPrtNum);axis image; colorbar
xlabel('\epsilon^2','FontSize',24);
ylabel('s','FontSize',24);
% set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
set(gca,'YTick',[1 35 69])
set(gca,'XTick',[34 67 length(E)])
set(gca,'YTickLabel',[num2str(floor(log10(s_vec(9)))),'|',num2str(ceil(log10(s_vec(35)))),'|',num2str(ceil(log10(s_vec(69))))]);
set(gca,'XTickLabel',[num2str(floor(log10(E(34)^2))),'|',num2str(floor(log10(E(67)^2))),'|',num2str(floor(log10(E(end)^2)))]);

figure;
axes('FontSize',14)
imagesc(cPrtNum);axis image; colorbar
xlabel('\epsilon^2','FontSize',24);
ylabel('s','FontSize',24);
% set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
set(gca,'YTick',[1 35 69])
set(gca,'XTick',[34 67 length(E)])
set(gca,'YTickLabel',[num2str(floor(log10(s_vec(9)))),'|',num2str(ceil(log10(s_vec(35)))),'|',num2str(ceil(log10(s_vec(69))))]);
set(gca,'XTickLabel',[num2str(floor(log10(E(34)^2))),'|',num2str(floor(log10(E(67)^2))),'|',num2str(floor(log10(E(end)^2)))]);


figure;
axes('FontSize',14)
imagesc(prtNum_W);axis image; colorbar
xlabel('\epsilon^2','FontSize',24);
ylabel('s','FontSize',24);
set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);


figure;
axes('FontSize',14)
imagesc(cPrtNum-qPrtNum);axis image; colorbar
xlabel('\epsilon^2','FontSize',24);
ylabel('s','FontSize',24);
set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);

%% Plots of Temperatures
figure;
iS=100;
axes('FontSize',14)
loglog(E.^2,T_eff_q(100,:),'k','LineWidth',2);
hold on;
loglog(E.^2,T_eff_c(100,:),'k','LineWidth',2);

% loglog(E,T_eff_c(99,:),'--r','LineWidth',3);
% loglog(E,T_eff_c(100,:),'r','LineWidth',3);

loglog(E.^2,T_eff_q(99,:),'--r','LineWidth',3);
loglog(E.^2,T_eff_c(99,:),'--b','LineWidth',3);
% loglog(E,T_eff_c(99,:),'--b','LineWidth',3);
loglog(E.^2,T_eff_q(70,:),'.-r','MarkerSize',7);
loglog(E.^2,T_eff_c(70,:),'.-b','MarkerSize',7);
loglog(E.^2,1/beta*ones(100,1),'g','LineWidth',3)
% loglog(E,T_eff_ss(99,:),'b');
% loglog(E,T_eff_ss(100,:),'.g','LineWidth',1);
legend(['Quantum    (s=1)   ';...
        'Stochastic (s=1)   ';...
        'Quantum    (s=0.85)';...
        'Stochastic (s=0.85)';...
        'Quantum    (s=0.01)';...
        'Stochastic (s=0.01)';...
        'T_B                '],'Location','NorthWest');
xlabel('\epsilon^2','FontSize',24);
ylabel('T','FontSize',24,'interpreter','latex');
axis([0, E(end)^2, 1,10^6])

figure;
axes('FontSize',14)
loglog(s_vec,T_eff_q(:,end),'r','LineWidth',3);
hold on
loglog(s_vec,T_eff_c(:,end),'b','LineWidth',3);
loglog(s_vec,1/beta*ones(100,1),'g','LineWidth',3)
xlabel('S','FontSize',24);
ylabel('T_{infty}','FontSize',24,'interpreter','tex');
legend(['Quantum    ';...
        'Stochastic ';...
        'T_B        '],'Location','NorthWest')
axis([s_vec(1),s_vec(end),1,10^6]);
%% T_infty
figure;
axes('FontSize',14)
loglog((s_vec),(T_eff_q(:,1)),'r','LineWidth',3);

axis([(s_vec(1)), (s_vec(end)),1, 10^6]);

xlabel('s', 'FontSize',24);
ylabel('T','FontSize',24);
hold on;
loglog((s_vec),(T_eff_q(:,2)),'r','LineWidth',6);
loglog((s_vec),(1/beta*ones(100,1)),'g','LineWidth',6)
ax1 = get(gca);
ax2 = axes('Position',ax1.Position,...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k','FontSize',14);
set(gca,'XTickLabel',[]);
line(log10(s_vec),prtNum_W(:,1),'Color','m','LineWidth',6,'LineStyle','--','Parent',ax2);
line(log10(s_vec),cPrtNum(:,1),'Color','b','LineWidth',2,'LineStyle','--','Parent',ax2);
line(log10(s_vec),qPrtNum(:,1),'Color','r','LineStyle','--','LineWidth',2,'Parent',ax2);

% line(log10(s_vec),prtNum_W(:,2),'Color','k','LineWidth',1,'LineStyle','--','Parent',ax2);
line(log10(s_vec),cPrtNum(:,2),'Color','b','LineWidth',6,'LineStyle','--','Parent',ax2);
line(log10(s_vec),qPrtNum(:,2),'Color','r','LineStyle','--','LineWidth',6,'Parent',ax2);

axis([log10(s_vec(1)), log10(s_vec(end)),0, 30]);

ylabel('PN','FontSize',24);

%  xlabel('log(s)','FontSize',24);
%  ylabel('log(T_{infty})','FontSize',24);
%  line(log10(s_vec),log10(T_eff_q),'Color','r','LineWidth',3,'Parent',ax2);
% axis([log10(s_vec(1)), log10(s_vec(end)),1, 4]);


%% Images of Temperatures
figure;
axes('FontSize',14)
h=imagesc((T_eff_q(1:end,:)),[17 50]);
% h=imagesc((T_eff_q(1:end-1,:)),[0 1000]);
set(h,'XData',[log10(E(1)^2) log10(E(end)^2)])
set(h,'YData',[log10(s_vec(1)) log10(s_vec(end))]);
axis tight; 

colorbar;
xlabel('\epsilon^2','FontSize',24);
ylabel('s','FontSize',24)
% 
% figure;
% axes('FontSize',14)
% h=imagesc(abs(T_eff_r(1:end,:)),[0 1000]);
% set(h,'XData',[log10(E(1)^2) log10(E(end)^2)])
% set(h,'YData',[log10(s_vec(1)) log10(s_vec(end))]);
% axis tight; 
% 
% colorbar;
% xlabel('\epsilon^2','FontSize',24);
% ylabel('s','FontSize',24)
% set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'YTickLabel',[num2str(log10(s_vec(1))),'|',num2str(log10(s_vec(length(s_vec)/2))),'|',num2str(s_vec(end),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
% 
% set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'YTickLabel',[num2str((s_vec(1))),'|',num2str((s_vec(length(s_vec)/2))),'|',num2str((s_vec(end)),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
% set(gca,'YTick',[1 35 69])
% set(gca,'XTick',[34 67 length(E)])
% set(gca,'YTickLabel',[num2str(floor(log10(s_vec(9)))),'|',num2str(ceil(log10(s_vec(35)))),'|',num2str(ceil(log10(s_vec(69))))]);
% set(gca,'XTickLabel',[num2str(floor(log10(E(34)^2))),'|',num2str(floor(log10(E(67)^2))),'|',num2str(floor(log10(E(end)^2)))]);



figure;
axes('FontSize',14)
h = imagesc((T_eff_c(1:end,:)),[17 50]);
set(h,'XData',[log10(E(1)^2) log10(E(end)^2)])
set(h,'YData',[log10(s_vec(1)) log10(s_vec(end))]);
axis tight; 

colorbar;
xlabel('\epsilon^2','FontSize',24);
ylabel('s','FontSize',24);
% set(gca,'YTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'XTick',[1 length(s_vec)/2 length(s_vec)])
% set(gca,'YTickLabel',[num2str(s_vec(1),'%.1e'),'|',num2str(s_vec(length(s_vec)/2),'%.1e'),'|',num2str(s_vec(end),'%.1e')]);
% set(gca,'XTickLabel',[num2str(E(1)^2,'%.1e'),'|',num2str((E(length(s_vec)/2)^2),'%.1e'),'|',num2str(E(end)^2,'%.1e')]);
% set(gca,'YTick',[1 35 69])
% set(gca,'XTick',[34 67 length(E)])
% set(gca,'YTickLabel',[num2str(floor(log10(s_vec(9)))),'|',num2str(ceil(log10(s_vec(35)))),'|',num2str(ceil(log10(s_vec(69))))]);
% set(gca,'XTickLabel',[num2str(floor(log10(E(34)^2))),'|',num2str(floor(log10(E(67)^2))),'|',num2str(floor(log10(E(end)^2)))]);

%% Figs. x_n
figure;
iS=301;
axes('FontSize',14)
stem(1:N-1, (randX(iS,:)),'b','Marker','^','LineWidth',2)
xlabel('n','FontSize',24)
ylabel('x_n','FontSize',24);