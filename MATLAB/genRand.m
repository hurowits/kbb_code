%Generate lognormal numbers
function randX = genRand(s_vec,randInd,N,strMethod)
% function randX = genRand(sigma_vec,randInd,N)
% close all
N=N-1;
DELTA = 1/N;

% m=1;
mu = 0;


%  s_vec=fliplr([1e-5,0.01,0.1,0.3,0.999999]);
% s_vec = exp(-15);logspace(-2,-eps,100);%[1-eps,0.5,0.1,0.01,0.001,0.0001,1e-7,1e-9]
% s_vec = exp(-0.01);
% exp(-15);[1-eps, exp(-2), exp(-15)];
% sigma_vec = [eps, 4];
% s_vec=1-eps;
% randInd = randperm(N);

for iS=1:length(s_vec)
    if(strcmp(strMethod,'sigma'))
        sigma = s_vec(iS);
    else
        s=s_vec(iS);
        
        sigma = sqrt(-log(s));
    end
    
    %     sigma = 300;%13.10.2010
    %     sigma = sigma_vec(iS);
    Zl(1) = mu;
    Zr(1) = mu;
    for q = 1:N/2
        Zr(q+1) = Zr(q)+1/(N*normpdf(Zr(q),mu,sigma));
    end
    for q = 1:N/2-1
        Zl(q+1) = Zl(q)-1/(N*normpdf(Zl(q),mu,sigma));
    end
    Z = unique([Zl Zr]);
    X = exp(Z);
    % figure;hist(Z);
    % figure;hist(X);title(num2str(s));
    % sparsity(iS) = mean(X)^2/mean(X.^2);
    randX(iS,:) = X(randInd);
    
end
%    save X.mat randX randInd
% end