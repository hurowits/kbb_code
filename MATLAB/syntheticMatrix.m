% (2) 
% Generate in Matlab 20x20 matrices {E_n} and {Vnm} 
% where Vnm has finite bandwidth. 
% Its in-band elements are generated as  
% 
% element = c1*random1 + c2*random*random2  
% 
% c2 >> c1 are constants.
% random1, random2 = random numbers in the range [-1,1]
% random = 1 with probability $p$ otherwise 0. 
% 
% REMEMBER to have Vnm=Vmn and without loss of 
% generality set the diagonal elements zero.
% function [gamma,V] = syntheticMatrix(N,c1,c2,p,bandLimit)
N = 20; 
E = diag(1:N);
V = zeros(N);

bandLimit = 5;
mask = 1-(triu(ones(N),bandLimit)+ tril(ones(N),-bandLimit));

c1 = 0.01; 
c2 = 1;
p=0.5;
random = (rand(N)>p);
random1 = (rand(N)-rand(N))/2;
random2 = (rand(N)-rand(N))/2;

lowerTriangle = tril(c1*random1 + c2*random.*random2);
V = (lowerTriangle+lowerTriangle.')*eye(N);
V = V - diag(diag(V));
V = V.*mask;
figure(1);imagesc(V);axis image;title('V');
drivingPowerSpectrum = ones(N); %white noise
gamma = drivingPowerSpectrum.*abs(V).^2;

%dpm/dt = gamma*(pm-pn) ==> dp/dt = g*p
g = gamma - gamma*(ones(N)-eye(N));

%Numerical solution of differential equation
% %Define a matrix PHI(t) such that p(t) = PHI(t)*p(0)
% 
%Initial conditions are random, since we are  looking for steady state
 P0 = rand(N,1);
 P0 = P0/sum(P0);
 P0 = zeros(N,1); P0(1)=1;
 t0 = 0;
 numIter = 100;
 P = zeros(N,numIter);
dt = 0.1;
T = 100;
numIter = round(T/dt);
for q= 0:numIter-1
    t = t0+q*dt;
    PHI = expm(g*t);
    P(:,q+1) = PHI*P0;
 
end
figure(2);plot(linspace(0,T,numIter),P);
% end