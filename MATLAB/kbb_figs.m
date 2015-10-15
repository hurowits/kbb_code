% load Temperature2
% load Wdot2

A = (T_eff_c);
B = (T_eff_q);
% C = (T_r);
% imagesc(A,[0,50])
% imagesc(B,[0,50])
E = Estrength_vec;
[X,Y] = meshgrid((E.^2),(s_vec));
%%
figure;
axes('FontSize',14)
hold on
contour(X,Y,A,[10.4:0.02:50]);
contour(X,Y,log10(A),[1.69:0.1:10],'k')
axis([-2 log10(E(100)^2) -6 0])
xlabel('log \epsilon^2','FontSize',24);
ylabel('log s','FontSize',24);
hold off
%%
figure;
axes('FontSize',14)
hold on
% contour(X,Y,B,[10.4:0.02:50]);
contour(X,Y,B,[10.04:0.02:50]);
contour(X,Y,log10(B),[1.69:0.1:10],'k')
% axis([log10(E(25)^2) log10(E(100)^2) (s_vec(1)) (s_vec(100))])
axis([-2 log10(E(100)^2) -6 0])
xlabel('log \epsilon^2 ','FontSize',24);
ylabel('log s ','FontSize',24);
hold off


% [X,Y] = meshgrid(log10(E.^2),log10(s_vec(1:end-1)));
%%
figure;
axes('FontSize',14)
hold on
% contour(X,Y,T_r,[10.4:0.02:50]);
contour(X,Y,T_mix,[104:0.2:500]);
contour(X,Y,log10(T_mix),[1.7:0.1:10],'k')
axis([log10(E(25)^2) log10(E(100)^2) (s_vec(1)) (s_vec(100))])

xlabel('log \epsilon^2 ','FontSize',24);
ylabel('log s ','FontSize',24);
hold off
%%
T_est = (N-1)./dispersion2(:,end) .*T_mix(:,end);

figure;
%  s_vec=exp(-s_vec.^2);
sig = sqrt(-log(s_vec));
sig=s_vec;
axes('FontSize',14)
loglog(sig,T_eff_q(:,end),'r','LineWidth',2);
hold on

plot(sig,T_mix(:,end),'b','LineWidth',2);
plot(sig, T_est,'g','LineWidth',2); 
plot(sig,1/beta*ones(length(s_vec),1),'k','LineWidth',2);

plot(sig,dispersion2(:,end),'c','LineWidth',2);
plot(sig,N*Delta*ones(length(s_vec),1),'--k','LineWidth',2);

xlabel('\sigma ','FontSize',24);
ylabel('log T ','FontSize',24);
  axis([1e-2 1e2 1e-1 1e6]);
legend(['T_{eff}         ';...
        'T_{mix}         ';...
        'expected T_{eff}';...
        'T_B             ';...
        'span[E_r]       ';...
        'span[E_n]       ']);
grid on;

hold off
% save Temperatures T_eff_q T_mix dispersion2 s_vec Estrength_vec beta N Delta
%%
figure;
iS=246;
axes('FontSize',14)
stem(1:N-1, Estrength_vec^2*(randX(iS,:)/mean(randX(iS,:))),'b','Marker','.','LineWidth',2);
hold on
plot(1:N-1,gamma*ones(N-1,1),'r','LineWidth',1);
% axis([1 N-1 0 1e4])
xlabel('n','FontSize',24)
ylabel('x_n','FontSize',24);
%%
iS=[73,146];[239,246,257]
figure;plot(E_r(:,iS),P_r(:,iS),'.');title(['sigma=',num2str(s_vec(iS))])
figure;plot(P_n_q(:,iS),'.');title(['sigma=',num2str(s_vec(iS))])
%%
%  s_vec=exp(-s_vec.^2);
T_est = (N)./dispersion2(:,end) .*T_mix;
figure;
axes('FontSize',14)

hold on

plot(s_vec,T_eff_q(:,end),'r','LineWidth',2);
plot(s_vec,T_mix(:,end),'b','LineWidth',2);
plot(s_vec, T_est,'g','LineWidth',2);
plot(s_vec,1/beta*ones(length(s_vec),1),'k','LineWidth',2);
plot(s_vec,dispersion2(:,end),'c','LineWidth',2);
plot(s_vec,N*Delta*ones(length(s_vec),1),'--k','LineWidth',2);

xlabel('\sigma ','FontSize',24);
ylabel('T ','FontSize',24);

set(gca,'xscale','log')
set(gca,'yscale','log')
axis([0.01 100 1 1000000])

legend(['T_{?inf}           ';...
        'T_{mix}            ';...
        'expected T_{inf}   ';...
        'T_B                ';...
        '\Delta[E_r]        ';...
        '\Delta[E_n]        ']);

grid on;
hold off