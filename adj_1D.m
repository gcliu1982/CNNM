function [x] = adj_1D(M)
%CONJA_1D 此处显示有关此函数的摘要
%   此处显示详细说明
[m,n] = size(M);
x = zeros(m,1);
for i = 1 : n
    x = x + circshift(M(:,i),1-i);
end

end

