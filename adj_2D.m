function [X] = adj_2D(M, isize, ksize)
%CONJA_1D 此处显示有关此函数的摘要
%   此处显示详细说明
X = zeros(isize(1),isize(2));
for j = 1: ksize(2)
    for i = 1 : ksize(1)
        Xi = M(:,(j-1)*ksize(1)+i);
        Xi = reshape(Xi,isize(1),isize(2));
        Xi = circshift(Xi, 1 - i, 1);
        Xi = circshift(Xi, 1 - j, 2);
        X = X + Xi;
    end
end

end

