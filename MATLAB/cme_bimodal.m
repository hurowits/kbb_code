function [W_dot,T] = cme_bimodal(a,b,N,eps,gamma)



x = (a*mod(1:N-1,2)+b*mod(0:N-2,2));
w = eps.'*x;

W_dot = mean(w./(w+gamma),2);
T = (mean(1./(w+gamma),2)).^-1;

