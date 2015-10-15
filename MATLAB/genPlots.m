pathOutput = '/Users/daniel/Desktop/University/Thesis/notes/dc'

a = 0.001;
b = 0.1;
N =25;
eps = logspace(-2,4,100);
gamma = 1;
Delta = 1;

[qW_dot,cW_dot,rho,alpha] = ness_bimodal(a,b,N,eps,gamma,Delta); %quantum, numerical solution

figure;
axes('FontSize',14)
loglog(eps,cW_dot,'b','LineWidth',1);
hold on
loglog(eps,qW_dot,'r','LineWidth',1);
   loglog(eps,abs(alpha),'g');
hold off;
grid;
xlabel('\epsilon / \gamma','FontSize',24);
ylabel('dW/dt','FontSize',24);

axis([0,eps(end),0,2])
% title(['N = ',num2str(N)]);
legend(['classical';'quantum  '],'Location','SouthEast');
%  print(gcf, '-depsc2', [pathOutput,'W_eps_cvq.eps']);
%   figure;loglog(eps,qW_dot./cW_dot);
%  xlabel('\epsilon / \gamma','FontSize',27);
%  ylabel('Cl / Qu','FontSize',27);
%  print(gcf, '-depsc2', 'W_eps_zero_lambda.eps')