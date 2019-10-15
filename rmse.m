function [e] = rmse(x,ref)
%RMSE 此处显示有关此函数的摘要
%   此处显示详细说明
x = x(:);
ref = ref(:);

e = sqrt(sum((x - ref).^2)/length(x));

end

