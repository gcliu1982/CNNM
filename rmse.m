function [e] = rmse(x,ref)
%RMSE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
x = x(:);
ref = ref(:);

e = sqrt(sum((x - ref).^2)/length(x));

end

