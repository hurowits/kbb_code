
% runs vs eps  (eps/gamma)
% runs for various Delta (i.e. change gamma/Delta)
% look for circumstances in which the "off diagonals" are large.

% clear all
% load X
N =25; %odd
Delta_vec = 1;%[0.1,1,10,100];
gamma_vec = 1000;%[0.001,1000 ];%[1,100];%logspace(-2,5,100);
beta_vec=0.1;%[1e-4,1e-2,1e-1,1/3,1];%[0.01,1,5];
%  s_vec = 1e-4;[1e-1,1e-2,1e-5,1e-7];

 Estrength_vec = logspace(-2,4,100);
%   Estrength_vec = 500
%[0.1, 1, 20, 500, 1000]

% Mpoints=length(Estrength_vec);
Mpoints=max([length(Estrength_vec),length(Delta_vec),length(beta_vec),length(s_vec)]);
alpha = zeros(Mpoints,1);
SupMatSize = ceil(N/2)*(1+N)-N;
Wdot = zeros(Mpoints,1);
rho = zeros(SupMatSize,Mpoints);
complex_rho = zeros((SupMatSize-N)/2,Mpoints);
rho_mat = zeros(N,N,Mpoints);
W = zeros(SupMatSize);
qWdot = zeros(Mpoints,length(s_vec));
cWdot = zeros(Mpoints,length(s_vec));
prtNum = zeros(Mpoints,length(s_vec));
rho_cl = zeros(N,Mpoints);
flag=zeros(Mpoints,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
for iGamma = 1:length(gamma_vec)
    for iDelta = 1:length(Delta_vec)
        for iBeta = 1:length(beta_vec)
            for iS=1:length(s_vec)
                
                x = randX(iS,:);
%                  a = 1;
%                  b = s_vec(iS);%1e-5;
%                  x = a*mod(1:N-1,2)+b*mod(0:N-2,2);
                %         x=ones(1,24);
                beta=beta_vec(iBeta);
                for q=1:Mpoints
                    
                    %% init model parameters
                    %% INPUT(w_a,w_b,N,Estrength,gamma,Delta)
                    wb=waitbar(q/Mpoints);
                    
                    Delta = Delta_vec(iDelta);
                    gamma = gamma_vec(iGamma);
                    
                    gammaPerp = gamma;
                    gammaUp = gamma*exp(-beta*Delta/2);
                    gammaDown = gamma*exp(beta*Delta/2);
                    
                    E_strength = Estrength_vec(q)*sqrt(mean(x));
                    
                    w = E_strength*sqrt(x);%matrix element. transition ration = matrixelement^2
                   
                    % w =sqrt(E_strength*( 1*mod(1:N,2)+10*mod(0:N-1,2) ));
                    %                 alpha(q,iS)=mean(w)^2/(mean(w.^2)+gamma^2+Delta^2);
                    alpha(q,iS)=geomean(w)/sqrt((gamma+mean(w))^2+Delta^2);
                    %
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
                    
                    fst_l=N+1;
                    sl_l = N+1;
                    for k=N-2:-2:1
                        complex_rho(sl_l:sl_l+k-1,q)=rho(fst_l:fst_l+k-1,q) +1i*rho(fst_l+k:fst_l+2*k-1,q);
                        fst_l=fst_l+2*k;
                        sl_l = sl_l + k;
                        
                    end
                    P = X(1:N).';
                    P_n(:,q) = P; %P in energy basis
                    %Reshape rho as matrix
                    rho_mat(:,:,q) = rho_vec2mat (complex_rho,P,N);
                    
                    %P in diagonal basis of rho
                    [v(:,:,q),D(:,:,q)] = eig(rho_mat(:,:,q));
                    P_d(:,q) = diag(D(:,:,q)); 
                    E_d(:,q) = abs(v(:,:,q)').^2*(1:N)';
                    
                   %P in basis that diagonalizes W
                   W_driving = diag(w,1)+diag(w,-1);
                   [r,D_r] = eig(W_driving);
                   P_r(:,q) = diag(inv(r)*rho_mat(:,:,q)*r);
                   E_r(:,q) = abs(r').^2*(1:N)'; 
                   dE_r(:,q) = sqrt(abs(r').^2*((1:N)'-E_r(:,q)));
                    % P in position basis
                    z = exp(1i*2*pi*(0:N-1)/N);
                    for iF = 1:N
                        F(iF,:) = z.^(iF-1);
                    end
                    F = F/sqrt(N);
                    rho_xx = inv(F)*rho_mat(:,:,q)*F;
                    
                    
                    P_x(:,q) = diag(rho_xx);
                    
                    %Effective Temperature
                    T_n = 1./log(P(1:N-1)./P(2:N));
                    T_eff_q(q) = 1./mean(1./T_n);
                    dT_eff_q(q) = 1./std(1./T_n);
                    %% Calculate heating rate
%                     figure(1);plot(abs(v(:,25,100)));
                    qPrtNum(q,iS) = sum(P.^2)^-1;
%                     S_tr(q,iS) = sqrt(sum(abs(complex_rho(:)).^2));
                    S_tr(q,iS) = sqrt(sum(X(N+1:end).^2));
                    % Bath normalized heating rate (as in the classical graph):
                    
                    qWdot(q,iS) = -sum(diff(P).*w(1:N-1).^2)/(Delta*gamma);
                    if (sum(P<0))
                        flag(q)=1;
                    end
                    
                    %%stochastic FGR
                    Pcl = linsolve([W_cl;ones(1,N)],[zeros(N,1);1]).';
                    cPrtNum(q,iS) = sum(Pcl.^2)^-1;
                    rho_cl(:,q)=Pcl;
                    %         cWdot(q,qq) =-sum(diff(Pcl).*w(1:N-1).^2)/(Delta*gamma);
                    cWdot(q,iS) =-sum(diff(Pcl).*w(1:N-1).^2)/(Delta*gamma);
                     %Effective Temperature
                    T_n =1./ log(Pcl(1:N-1)./Pcl(2:N));
                    T_eff_c(q) = 1./mean(1./T_n);
                    %                 a=w(1)^2 ;
                    %                 AA = ((a+gamma*exp(-beta*Delta/2))/(a+gamma*exp(beta*Delta/2))).^(0:N-1);
                    %                 p0 = sum(AA)^(-1);
                    %                 P_cl_an(:,q) = AA*p0;
                    %                 caWdot(q) = -sum(diff(P_cl_an(:,q).')*a)/(Delta^2*gamma);
                    %
                    %
                end
                close(wb);
                % END OF MAIN LOOP
                % OUTPUT(qWdot,cWdot,rho,alpha)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%% FIGURES
                
                %Fig1 - $rho_cl$ and $rho_qm$ profiles for the last run
                
                %Fig2a - heating rate vs run
                %fig2b - error = -log(pN/p0)  vs run
                %fig2c - max(rho offdiagonal)/p0  vs run
                
                
                %%%%%%%%%%%%%
                %      figure;
                %      plot(rho(1:N,1),'-r')
                %      hold on
                % %     %
                %        plot(abs(complex_rho(N+1:end,end)),'-.b') %quantum numeric
                %     plot(rho(1:N,end),'r') %quantum numeric
                %     plot(P_cl_an(:,end),'b'); %classical analytic
                %     xlabel('n','FontSize',24);
                %     ylabel('p_n','FontSize',24);
                %     hold off
                % title('Quantum(red) vs. Stochastic (blue) Probabilities')
                
               
                %         figure(1);imagesc(log10(abs(rho_mat)+abs(rho_mat).'+diag(P)));axis image;colormap gray;colorbar;
                %     title('Density Matrix');
                %          print(gcf, '-depsc2', ['rho_mat_',num2str(s),'_',num2str(beta),'.eps']);
                
            
            
            figure;
            axes('FontSize',14)
            E=Estrength_vec*sqrt(mean(x));
            loglog(E,cWdot(:,iS),'b-','LineWidth',1);
            
            hold on
            %     loglog(E,caWdot,'b','LineWidth',1);
            loglog(E,qWdot(:,iS),'r-','LineWidth',1);
            
            
            loglog(E,E.^2*mean(x)*beta/gamma,'g--','LineWidth',1);    %LRT
            plot(E,ones(length(E),1)*beta,'g--','LineWidth',1); %bath limited saturation
            
            
            plot(ones(length(E),1)*sqrt((gamma)/mean(x)),linspace(1e-8,100,length(E)),'k--','LineWidth',1); %cross over to SLRT
            plot(ones(length(E),1)*sqrt((gamma)/min(x)),linspace(1e-8,100,length(E)),'k--','LineWidth',1); %crossover to saturation
            
%             plot(ones(length(E),1)*(mean(sqrt(x))-sqrt(mean(sqrt(x))^2*(1+4*gamma^2)-4*geomean(sqrt(x))^2*(gamma^2+Delta^2)))/2/(geomean(sqrt(x))^2-mean(sqrt(x))^2),linspace(1e-8,100,length(E)),'g--','LineWidth',1); %crossover to saturation
            
            
            text(sqrt((gamma)/mean(x))/10,1e-4,'LRT','FontSize',14,'interpreter','latex');
            text(sqrt((gamma)/min(x))/10,1e-4,'SLRT','FontSize',14,'interpreter','latex');
            text(sqrt((gamma)/min(x))*10^0.5,1e-4,'Saturation','FontSize',14,'interpreter','latex');
            %     loglog(E,cWdot(:,2:end),'LineWidth',3);
            %     loglog(E,qWdot(:,2:end),'--','LineWidth',3);
            
            %                  plot(E,prtNum,'LineWidth',1);
            %              loglog(E,S_tr(:,2:end),'LineWidth',1);
            %         loglog(E,alpha(:,1),'.','LineWidth',1);
%             loglog(E,cPrtNum(:,iS),'b','LineWidth',1);
%             loglog(E,qPrtNum(:,iS),'r','LineWidth',1);
%             loglog(E,S_tr(:,iS),'k','LineWidth',1);
%             loglog(E,alpha(:,1),'k','LineWidth',2);
            %                 loglog(E,alpha,'o','LineWidth',1);
            hold off;
            %         grid;
            xlabel('\epsilon','FontSize',24);
            ylabel('Rate of heating','FontSize',24,'interpreter','latex');
            legend(['Stochastic  ';...
                    'Quantum     '],'Location','SouthEast');
            axis([0,E(end),1e-8,1e2]);
            
            d_sorted = sortrows([E_d,P_d]);
            r_sorted = sortrows([E_r,P_r]);
            figure;plot(1:N,Pcl,...
                        1:N,P_n,...
                        d_sorted(:,1),d_sorted(:,2),...
                        1:N,P_x,...
                        r_sorted(:,1),r_sorted(:,2));
%                 figure;plot(Pcl,1:N,...
%                 P_n,1:N,...
%                 d_sorted(:,2),d_sorted(:,1),...
%                 P_x,1:N);hold on
%                 errorbar(r_sorted(:,2),r_sorted(:,1),dE_r/2,'m');
            
            

            legend(['Stochastic      ';...
                    'Momentum        ';...
                    '\rho Diagonal   ';...
                    'Position        ';...
                    'Driving Diagonal']);
                title(['\epsilon = ',num2str(Estrength_vec*sqrt(mean(x)))]);
%             print(gcf,
%             '-depsc2',['regimes','_',num2str(beta),'_',num2str(Delta),'_',num2str(s_vec(iS) ),'.eps']);
%             saveas(gcf,['regimes2','_',num2str(gamma),'_',num2str(s_vec(iS) ),'.fig']);
%                close all
            end
        end
    end
end