function [a] = shrink_l2(b,lambda)
%SHRINK_L2 此处显示有关此函数的摘要
%   此处显示详细说明
bnorm = abs(b);
a = ((bnorm - lambda)./bnorm).*b;
inds = bnorm <= lambda;
a(inds) = 0;

end

