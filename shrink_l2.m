function [a] = shrink_l2(b,lambda)
%SHRINK_L2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
bnorm = abs(b);
a = ((bnorm - lambda)./bnorm).*b;
inds = bnorm <= lambda;
a(inds) = 0;

end

