function [x] = inexact_alm_cnm_1D(y, g,  n, lambda, tol, maxIter)
%%% This matlab code implements the (inexact) augmented Lagrange multiplier method for
%%% min_{x}  |A(x)|_* + 0.5*n*lambda |P_Omega(x - y)|_2^2
%%% A(x) -- linear operator that returns the convolutional matrix of x
%%% | |_* -- the nuclear norm of a matrix
%% parameters:
% y -- m x 1 vector of observations (required input),with the missing entries being filled with zero
% g -- m x 1 vector with 1 and 0 denoting the observed and missing entries respecitvely (required input) 
% n -- kernel size (required input) 
% lambda -- weight parameter
%      - DEFAULT 1000
% tol -- tolerance for stopping criterion.
%     -- DEFAULT 1e-7 
% maxIter -- maximum number of iterations
%         -- DEFAULT 10000
%
m = length(y);

if nargin < 6 || isempty(maxIter)
    maxIter = 10000;
end

if nargin < 5 || isempty(tol)
    tol = 1e-7;
end

if nargin < 4 || isempty(lambda)
    lambda = 1000;
end

fnorm = sqrt(n)*norm(y,2);

% initialize
W = zeros(m,n); %Lagrange multiplier 
x = y;

%parameters
obsrate = sum(g>0.5)/n;
mu = obsrate/fnorm; % this one can be tuned
rho = 1.05; % this one should be chosen carefully

iter = 0;
%% loop
Ax = cconv1mtx(x,n);
while iter < maxIter       
   iter = iter + 1;   
    %update L
    temp = Ax + W./mu;
    [U,S,V] = svd(temp, 'econ');
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    diagS = max(0,diagS - 1/mu);
    if svp < 0.5 %svp = 0
        svp = 1;
    end
    diagS = diagS(1:svp);
    U = U(:,1:svp);
    V = V(:,1:svp);
    L =U*bsxfun(@times,V',diagS);  
        
    %update x
    temp = adj_1D(L - W./mu)./n;
    temp = lambda*y+mu*temp;
    diagA = lambda*g + mu;
    x = temp./diagA;
    Ax = cconv1mtx(x,n);
     
    H = Ax - L;
    %% stop Criterion    
    stopCriterion = norm(H, 'fro') / fnorm;
    if iter == 1 || mod(iter,100)==0 || stopCriterion < tol
        disp(['Iteration' num2str(iter) ', mu ' num2str(mu) ', stopCriterion ' num2str(stopCriterion)]);
    end
    
    if stopCriterion < tol
        break;
    else
        W = W + mu*H;
        mu = min(mu*rho,10^10);
    end    
end
end

