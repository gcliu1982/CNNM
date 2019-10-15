function [X] = inexact_alm_cnm_3D(Y, G,  ksize, lambda, tol, maxIter)
%%% This matlab code implements the (inexact) augmented Lagrange multiplier method for
%%% min_{x}  |A(X)|_* + 0.5*ksize(1)*ksize(2)*ksize(3)*lambda |P_Omega(x - y)|_2^2
%%% A(X) -- linear operator that returns the convolutional matrix of X
%%% | |_* -- the nuclear norm of a matrix
%% parameters:
% Y -- m x n x l vector of observations (required input),with the missing entries being filled with zero
% G -- m x n x l vector with 1 and 0 denoting the observed and missing entries respecitvely (required input) 
% ksize -- kernel size (required input) 
% lambda -- weight parameter
%      - DEFAULT 1000
% tol -- tolerance for stopping criterion.
%     -- DEFAULT 1e-7 
% maxIter -- maximum number of iterations
%         -- DEFAULT 1000
%
[m, n, l] = size(Y);
k1 = ksize(1);
k2 = ksize(2);
k3 = ksize(3);

if nargin < 6 || isempty(maxIter)
    maxIter = 10000;
end

if nargin < 5 || isempty(tol)
    tol = 1e-7;
end

if nargin < 4 || isempty(lambda)
    lambda = 1000;
end

fnorm = sqrt(k1*k2*k3)*norm(Y(:));

% initialize
W = zeros(m*n*l,k1*k2*k3); %Lagrange multiplier 
X = Y;

%parameters
obsrate = sum(G(:)>0.5)/(m*n*l);
mu = obsrate/fnorm; % this one can be tuned
rho = 1.05; % this one should be chosen carefully

iter = 0;
%% loop
Ax = cconv3mtx(X,ksize);
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
    temp = adj_3D((L - W./mu)./(k1*k2*k3),[m n l],ksize);
    temp = lambda*Y+mu*temp;
    diagA = lambda*G + mu;
    X = temp./diagA;
    Ax = cconv3mtx(X,ksize);
    
    H = Ax - L;
    %% stop Criterion    
    stopCriterion = norm(H(:)) / fnorm;
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

