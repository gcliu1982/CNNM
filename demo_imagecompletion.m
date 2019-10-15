function [] = demo_imagecompletion()
%EXP3 此处显示有关此函数的摘要
%   此处显示详细说明
close all;
X0 = imread('lena.png');
m = 200;
X0 = imresize(X0,[m m]);
X0 = im2double(X0);
pval = max(X0(:));

ms = 5:5:195;

Omega = ones(m,m);
Omega(:,ms) = 0;
Omega(ms,:) = 0;

inds = Omega>0.1;
disp(['missing rate = ' num2str(sum(~inds(:))/(m*m))]);
Y = X0.*Omega;
subplot(1,3,1);
imshow(Y);
title('input');

%run DFT_L1
tic;
X_hat = inexact_alm_dft_nD(Y,Omega);
toc;
subplot(1,3,2);
imshow(X_hat);
title('DFT_L1');
acc = psnr(X_hat(~inds),X0(~inds),pval);
disp(['PSNR (DFT_L1) = ' num2str(acc)]);


% run_CNNM
ksize = [13 13];
tic;
X_hat = inexact_alm_cnm_2D(Y, Omega,ksize);
toc;
subplot(1,3,3);
imshow(X_hat);
title('CNNM');
acc = psnr(X_hat(~inds),X0(~inds),pval);
disp(['PSNR (CNNM) = ' num2str(acc)]);

end