function [A] = cconv1mtx(x,n)
%CCONV1MTX 此处显示有关此函数的摘要
%   此处显示详细说明
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

