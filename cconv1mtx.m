function [A] = cconv1mtx(x,n)
%CCONV1MTX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
m = length(x);
A = zeros(m,n);
xi = x;
for i = 1:n
    A(:,i) = xi;
    if i < n
        xi = circshift(xi,1);
    end
end

end

