function [A] = cconv2mtx(X,ksize)
%CCONV1MTX 此处显示有关此函数的摘要
%   此处显示详细说明
[m, n] = size(X);
A = zeros(m*n,ksize(1)*ksize(2));
Xj = X;
for j = 1:ksize(2)
    Xi = Xj;
    for i = 1:ksize(1)
        A(:,(j-1)*ksize(1)+i) = Xi(:);
        if i < ksize(1)
            Xi = circshift(Xi,1,1);
        end
    end
    if j < ksize(2)
        Xj = circshift(Xj,1,2);
    end
end

end

