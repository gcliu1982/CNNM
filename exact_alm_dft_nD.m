function [X] = exact_alm_dft_nD(Y, G, lambda, tol, maxIter)
%%% This matlab code implements the (exact) augmented Lagrange multiplier method for
%%% min_{X}  |A(X)|_1 + 0.5*lambda |P_Omega(X - Y)|_F^2
%%% A(X) -- linear operator that returns the DFT of a signal
%% parameters:
% Y -- m x n x l tensor of observations (required input),with the missing entries being filled with zero
% G -- m x n x l tensor with 1 and 0 denoting the observed and missing entries respecitvely (required input) 
% lambda -- weight parameter
%      - DEFAULT 1000
% tol -- tolerance for stopping criterion.
%     -- DEFAULT 1e-7 
% maxIter -- maximum number of iterations
%         -- DEFAULT 1000
%
sy = size(Y);
m = prod(sy);

if nargin < 5 || isempty(maxIter)
    maxIter = 1000;
end

if nargin < 4 || isempty(tol)
    tol = 1e-7;
end

if nargin < 3 || isempty(lambda)
    lambda = 1000;
end

maxIter_primal = 10000;
fnorm = norm(Y(:));

% initialize
W = zeros(sy); %Lagrange multiplier 
X = Y;
L = zeros(sy);
Fx = fftn(X);

%parameters
tolProj = 1e-6 * fnorm;
obsrate = sum(G(:)>0.5)/m;
mu = obsrate/fnorm; % this one can be tuned
rho = 5;          % this one can be tuned

iter = 0;
%% loop
while iter < maxIter       
    iter = iter + 1;
    
    %% solve the primal problem by alternating minimization
    primal_iter = 0;
    while primal_iter < maxIter_primal
        primal_iter = primal_iter + 1;
        temp_L = L;
        
        %update L
        temp = Fx + W./mu;
        L = shrink_l2(temp, 1/mu);
        
        %update X
        temp = ifftn(mu.*L - W);
        temp = lambda.*Y + temp;
        diagA = lambda*G + mu;
        X = temp./diagA;
        Fx = fftn(X);
        if norm(L(:) - temp_L(:)) < tolProj
            break;
        end
    end
    H = Fx - L;
    %% stop Criterion    
    stopCriterion = norm(H(:)) / fnorm;
    disp(['Iteration' num2str(iter) '(' num2str(primal_iter) '), mu ' num2str(mu) ', stopCriterion ' num2str(stopCriterion)]);
    
    
    if stopCriterion < tol
        break;
    else
        W = W + mu*H;
        mu = mu*rho;
    end    
end
X = real(X);
end

