function [X] = adj_3D(M, isize, ksize)
%CONJA_1D 此处显示有关此函数的摘要
%   此处显示详细说明
X = zeros(isize(1),isize(2),isize(3));
for k = 1: ksize(3)
    for j = 1: ksize(2)
        for i = 1 : ksize(1)
            Xi = M(:,((k-1)*ksize(2)+(j-1))*ksize(1)+i);
            Xi = reshape(Xi,isize(1),isize(2),isize(3));
            Xi = circshift(Xi, 1 - i, 1);
            Xi = circshift(Xi, 1 - j, 2);
            Xi = circshift(Xi, 1 - k, 3);
            X = X + Xi;
        end
    end
end
end

