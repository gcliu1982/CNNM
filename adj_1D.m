function [x] = adj_1D(M)
%CONJA_1D �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[m,n] = size(M);
x = zeros(m,1);
for i = 1 : n
    x = x + circshift(M(:,i),1-i);
end

end

