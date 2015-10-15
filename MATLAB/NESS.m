%Calculate and plot NESS heating rate and temperature 
close all
clear all
% Model parameters
gamma = 1; %bath coupling
eps = logspace(-2,2,100); %driving intensity
N = 1000; %number of levels

s = [0.001,0.1,0.3,1];%sparsity
textLegend = char(zeros(size(s,2),256));
textLegend(:,1:2)=repmat('s=',size(s,2),1);
for qq = 1:size(s,2)
   
    m = 1; %mean of log-normal distribution
    mu = log(sqrt(s(qq))*m);
    sigma = sqrt(-log(s(qq)));
    
    x = lognrnd(mu,sigma,1,N);
    w = eps.'*x;
    
    dotW(:,qq) = mean(w./(w+gamma),2);
    T(:,qq) = (mean(1./(w+gamma),2)).^-1;
    textLegend(qq,3:3+size(num2str(s(qq)),2)-1) = num2str(s(qq));
end
%reference curve (w not random)
% w = eps.';    
% dotW(:,qq+1) = w./(w+gamma);
% T(:,qq+1) = (1./(w+gamma)).^-1;
% 
%   

figure(1);
axes('FontSize',20);

plot(eps,dotW(:,1),'-r','LineWidth',2);
hold on
plot(eps,dotW(:,2),'-.g','LineWidth',2);
plot(eps,dotW(:,3),':b','LineWidth',2);
plot(eps,dotW(:,4),'k','MarkerSize',20)
hold off

grid
axis([0,eps(end),0,1.2])
xlabel('\epsilon','FontSize',27)
ylabel('dW/dt','FontSize',27)

legend(textLegend,'Location','SouthEast')

print(gcf, '-depsc2', '/Users/daniel/Desktop/Thesis/code/slrt/figures/dotWvsEps_logn.eps')

figure(2);
axes('FontSize',20);
loglog(eps,T(:,1),'-r','LineWidth',2);
hold on
loglog(eps,T(:,2),'-.g','LineWidth',2);
loglog(eps,T(:,3),':b','LineWidth',2);
loglog(eps,T(:,4),'k','MarkerSize',20)
hold off

grid
 axis([0,eps(end),0,100])
xlabel('\epsilon','FontSize',27)
ylabel('T/T_B','FontSize',27)
legend(textLegend,'Location','NorthWest')
 
print(gcf, '-depsc2', '/Users/daniel/Desktop/Thesis/code/slrt/figures/TvsEps_logn.eps')

