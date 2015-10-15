%% Initialization
clear all
close all
N =3;
beta_vec  = 1;
W = zeros(N^2);
W_cl = zeros(N);
beta = beta_vec;

w_beta =1000;
gamma = w_beta;
Energy = [0,0,0];
V_0 = [0.1,1,10];
Estrength_vec = linspace(0,30,1e3); %"Epsilon"

for iE = 1:length(Estrength_vec)
    
W = zeros(N^2);
W_cl = zeros(N);
W_Im = zeros(N);
W_Re = zeros(N);

    Estrength = Estrength_vec(iE);
    V = Estrength * V_0;
    %%  "W_cl"
    for ii = 1:N
        for jj = 1:N
            if (ii~=jj)
                if((ii==1 && jj==2)||(ii==2 && jj==1)) kk=3;
                elseif((ii==1 && jj==3)||(ii==3 && jj==1))kk=2;
                elseif((ii==2 && jj==3)||(ii==3 && jj==2))kk=1;
                end
                w(ii,jj) = V(kk)^2 + w_beta*exp(0.5*beta*(Energy(ii)-Energy(jj)))/cosh(0.5*beta*(Energy(ii)-Energy(jj)));
                w(jj,ii) = V(kk)^2 + w_beta*exp(-0.5*beta*(Energy(ii)-Energy(jj)))/cosh(0.5*beta*(Energy(ii)-Energy(jj)));
            end
            
        end
        
    end
    
    
    W_cl(1,1) = -(w(1,2)+w(1,3));
    W_cl(2,2) = -(w(2,1)+w(2,3));
    W_cl(3,3) = -(w(3,1)+w(3,2));
    
    W_cl(1,2) = w(2,1);
    W_cl(1,3) = w(3,1);
    W_cl(2,3) = w(3,2);
    
    W_cl(2,1) = w(1,2);
    W_cl(3,1) = w(1,3);
    W_cl(3,2) = w(2,3);

    
    
    %% W_R
    W_Re(1,1) = -(V(1)^2 + V(2)^2);
    W_Re(2,2) = -(V(2)^2 + V(3)^2);
    W_Re(3,3) = -(V(3)^2 + V(1)^2);
    
    W_Re(1,2) = V(1)*V(3);
    W_Re(1,3) = V(2)*V(3);
    W_Re(2,3) = V(1)*V(2);
    
    
    W_Re(2,1) = V(1)*V(3);
    W_Re(3,1) = V(2)*V(3);
    W_Re(3,2) = V(1)*V(2);
    
    W_Re = 0.5*W_Re;
    
    %% W_I
    
    W_Im(1,1) = -(V(1)^2 + V(2)^2+4*V(3)^2);
    W_Im(2,2) = -(4*V(1)^2 + V(2)^2 + V(3)^2);
    W_Im(3,3) = -(V(1)^2+4*V(2)^2 + V(3)^2);
    
    W_Im(1,2) = 3*V(1)*V(3);
    W_Im(1,3) = 3*V(2)*V(3);
    W_Im(2,3) = 3*V(1)*V(2);
    
    W_Im(2,1) = 3*V(1)*V(3);
    W_Im(3,1) = 3*V(2)*V(3);
    W_Im(3,2) = 3*V(1)*V(2);
    W_Im = 0.5*W_Im;
    
    %% W_Delta
    W_D = zeros(N);
    W_D(1,1) = Energy(1)-Energy(2);
    W_D(2,2) = Energy(2)-Energy(3);
    W_D(3,3) = Energy(3)-Energy(1);
    %% Lambda
    Lambda(1,:) = V(1)*V(2) * [-1,-1,2];
    Lambda(2,:) = V(2)*V(3) * [2,-1,-1];
    Lambda(3,:) = V(1)*V(3) * [-1,2,-1];
    Lambda = Lambda*0.5;
    %% WWW
%     W_cl=2*W_cl;Lambda=2*Lambda;W_Im=2*W_Im;W_Re=2*W_Re;
    W = [W_cl     ,    2*Lambda'         , zeros(3);...
        Lambda    , -gamma*eye(3) + W_Re , W_D;...
        zeros(3)  ,      -W_D            , -gamma*eye(3) + W_Im];
    
    %% Solve steady state equation
    NormCondition = zeros(1,N^2);
    NormCondition(1,1:N)=ones(1,N);
    B = zeros(N^2+1,1);
    B(end)=1;
    
    % X -> {p,re[rho]...,im[rho]...}
    %     [X] = linsolve([WW;NormCondition],B);
    [X,R] = linsolve([W;NormCondition],B);
    RHO(iE,:) = X; 
    
    %stochastic current formula with P from QME
%     sqI_12(iE) = X(1)*w(1,2) - X(2)*w(2,1);
%     sqI_23(iE) = X(2)*w(2,3) - X(3)*w(3,2);
%     sqI_31(iE) = X(3)*w(3,1) - X(1)*w(1,3);
    
    %quantum current + thermal current
    
   qI_12(iE) = -V(3)^2*(X(1)-X(2)) - V(1)*V(3)*X(6) + V(2)*V(3)*X(5);
   qI_23(iE) = -V(1)^2*(X(2)-X(3)) - V(1)*V(2)*X(4) + V(1)*V(3)*X(6);
   qI_31(iE) = -V(2)^2*(X(3)-X(1)) + V(1)*V(2)*X(4) - V(2)*V(3)*X(5);
   
    tI_12(iE) =  X(1)*(w(1,2)-V(3)^2)-X(2)*(w(2,1)-V(3)^2);
    tI_23(iE) =  X(2)*(w(2,3)-V(1)^2)-X(3)*(w(3,2)-V(1)^2);
    tI_31(iE) =  X(3)*(w(3,1)-V(2)^2)-X(1)*(w(1,3)-V(2)^2);
    
    %total current:
    qtI_12(iE) = qI_12(iE) - tI_12(iE);
    qtI_23(iE) = qI_23(iE) - tI_23(iE);
    qtI_31(iE) = qI_31(iE) - tI_31(iE);
    
    %stochastic current formula with P (stochastic)
    [X_cl,R] = linsolve([W_cl;[1 1 1]],[0;0;0;1]);
    P_cl(iE,:) = X_cl;
    cI_12(iE) = X_cl(1)*w(1,2) - X_cl(2)*w(2,1);
    cI_23(iE) = X_cl(2)*w(2,3) - X_cl(3)*w(3,2);
    cI_31(iE) = X_cl(3)*w(3,1) - X_cl(1)*w(1,3);
    
   %driving current
    cdI_12(iE) = (X_cl(1) - X_cl(2))*V(3)^2;
    cdI_23(iE) = (X_cl(2) - X_cl(3))*V(1)^2;
    cdI_31(iE) = (X_cl(3) - X_cl(1))*V(2)^2;
    
    %thermal current
    ctI_12(iE) =  X_cl(1)*(w(1,2)-V(3)^2)-X_cl(2)*(w(2,1)-V(3)^2);
    ctI_23(iE) =  X_cl(2)*(w(2,3)-V(1)^2)-X_cl(3)*(w(3,2)-V(1)^2);
    ctI_31(iE) =  X_cl(3)*(w(3,1)-V(2)^2)-X_cl(1)*(w(1,3)-V(2)^2);
    
    SMF = log(w(1,2)*w(2,3)*w(3,1))-log(w(2,1)*w(3,2)*w(1,3));
    
    %Theoretical solution
    I(iE) = (w(1,2)*w(2,3)*w(3,1) - w(2,1)*w(3,2)*w(1,3));
    I(iE) = I(iE) / (w(3,1)*w(2,1) + w(3,2)*w(2,1)+w(2,3)*w(3,1) + ...
        w(2,1)*w(1,3) + w(1,2)*w(2,3) + w(2,3)*w(1,3)+...
        w(3,1)*w(1,2)+w(1,3)*w(3,2)+w(3,2)*w(1,2));
end
%% Figures 
figure;
axes('FontSize',14)
plot(Estrength_vec.^2,RHO(:,1),Estrength_vec.^2,RHO(:,2),Estrength_vec.^2,RHO(:,3));
xlabel('\epsilon^2');
ylabel('p_n (quantum)');
legend('p_1','p_2','p_3');
% 
% 
% figure;
% axes('FontSize',14)
% plot(Estrength_vec,RHO(:,4),Estrength_vec,RHO(:,5),Estrength_vec,RHO(:,6));
% xlabel('\epsilon^2');
% ylabel('\rho^x');
% legend('\rho_{12}^x','\rho_{23}^x','\rho_{31}^x');
% 
% figure;
% axes('FontSize',14)
% plot(Estrength_vec,RHO(:,7),Estrength_vec,RHO(:,8),Estrength_vec,RHO(:,9));
% xlabel('\epsilon^2');
% ylabel('\rho^y ');
% legend('\rho_{12}^y','\rho_{23}^y','\rho_{31}^y');
% 
% % 
% 
figure;
axes('FontSize',14)
plot(Estrength_vec.^2,P_cl(:,1),Estrength_vec.^2,P_cl(:,2),Estrength_vec.^2,P_cl(:,3));
xlabel('\epsilon^2');
ylabel('p_n (stochastic)');
legend('p_1','p_2','p_3');

figure;
axes('FontSize',14)
plot(Estrength_vec.^2,qI_12,Estrength_vec.^2,qI_31,Estrength_vec.^2,qI_23);
xlabel('\epsilon^2');
ylabel('Quantum Current');
legend('I_{12}','I_{31}','I_{23}');
title('I_q(rho) - Driving current');


figure;
axes('FontSize',14)
plot(Estrength_vec.^2,tI_12,Estrength_vec.^2,tI_23,Estrength_vec.^2,tI_31);
xlabel('\epsilon^2');
ylabel('Thermal Current');
legend('I_{12}^{bath}','I_{31}^{bath}','I_{23}^{bath}');

figure;
axes('FontSize',14)
plot(Estrength_vec.^2,cdI_12,Estrength_vec.^2,cdI_23,Estrength_vec.^2,cdI_31);
xlabel('\epsilon^2');
ylabel('stochastic driving Current');
legend('I_{12}','I_{23}','I_{31}');



figure;
axes('FontSize',14)
plot(Estrength_vec.^2,ctI_12,Estrength_vec.^2,ctI_23,Estrength_vec.^2,ctI_31);
xlabel('\epsilon^2');
ylabel('stochastic Thermal Current');
legend('I_{12}^{bath}','I_{23}^{bath}','I_{31}^{bath}');


% % 
% figure;
% axes('FontSize',14)
% plot(Estrength_vec.^2,qtI_12,Estrength_vec.^2,qtI_23,Estrength_vec.^2,qtI_31);
% xlabel('\epsilon^2');
% ylabel('I');
% title('I_q + I_t');
% 
% 
% figure;
% axes('FontSize',14)
% plot(Estrength_vec.^2,-cI_12,Estrength_vec.^2,-I);
% xlabel('\epsilon^2');
% ylabel('I');
% legend(['stochastic ';...
%         'theoretical']);
    
    
    
figure;
axes('FontSize',14)
plot(Estrength_vec.^2,I,'og',Estrength_vec.^2,cI_12,'b',Estrength_vec.^2,-qtI_12,'r');
hold on
plot(Estrength_vec.^2,cdI_12,Estrength_vec.^2,cdI_23,Estrength_vec.^2,cdI_31);
plot(Estrength_vec.^2,ctI_12,Estrength_vec.^2,ctI_23,Estrength_vec.^2,ctI_31);
xlabel('\epsilon^2');
ylabel('I');
legend(['Theoretical';...
        'stochastic ';...
        'quantum    ']);
    
    
figure;
axes('FontSize',14)
plot(Estrength_vec.^2,cdI_12,'r',Estrength_vec.^2,cdI_23,'g',Estrength_vec.^2,cdI_31,'b');hold on
plot(Estrength_vec.^2,ctI_12,'.r',Estrength_vec.^2,ctI_23,'.g',Estrength_vec.^2,ctI_31,'.b');
plot(Estrength_vec.^2,cI_12,'-m');
xlabel('\epsilon^2');
ylabel('I');
legend(['D I_{12}';...
        'D I_{23}';... 
        'D I_{31}';...
        'T I_{12}';...
        'T I_{23}';...
        'T I_{31}';...
        'total   ']);