function [A] = cconv3mtx(X,ksize)
%CCONV1MTX 此处显示有关此函数的摘要
%   此处显示详细说明
[m, n, l] = size(X);
A = zeros(m*n*l,ksize(1)*ksize(2)*ksize(3));
Xk = X;
for k = 1:ksize(3)
    Xj = Xk;
    for j = 1:ksize(2)
        Xi = Xj;
        for i = 1:ksize(1)
            A(:,((k-1)*ksize(2)+(j-1))*ksize(1)+i) = Xi(:);
            if i < ksize(1)
                Xi = circshift(Xi,1,1);
            end
        end
        if j < ksize(2)
            Xj = circshift(Xj,1,2);
        end
    end
    if k < ksize(3)
        Xk = circshift(Xk,1,3);
    end
end

end

