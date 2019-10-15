function [] = demo_videoprediction()
%EXP3 此处显示有关此函数的摘要
%   此处显示详细说明
close all;
ids = 105:166;

disp('reading images ...');
dpath = 'highway';
m = 50;
n = 50;
l = length(ids);
D = zeros(m, n, l);
for i = 1:l
    id = ids(i);
    ipath = [dpath '/' num2str(id) '.png'];
    I = imread(ipath);
    I = im2double(I);
    D(:,:,i) = I;
end
mids = 57:62; %unknown future data
nbmiss = length(mids);
Omega = ones(m,n,l);
Omega(:,:,mids) = 0;
gt = D(:,:,mids); %ground truth
pval = 1;
Y = D.*Omega;

figure;
%show ground truth
for i = 1:nbmiss
    subplot(3,nbmiss,i);
    imshow(gt(:,:,i));
end

%run DFT_L1
tic;
X_hat = inexact_alm_dft_nD(Y, Omega);
X_hat = X_hat(:,:,mids);
toc;
for i = 1:nbmiss
    acc = psnr(X_hat(:,:,i),gt(:,:,i),pval);
    disp([num2str(i) 'st frame, PSNR = ' num2str(acc)]);
    subplot(3,nbmiss,nbmiss+i);
    imshow(X_hat(:,:,i));
end

%run CNNM
ksize = [13 13 13];
tic;
X_hat = inexact_alm_cnm_3D(Y, Omega,ksize);
X_hat = X_hat(:,:,mids);
toc;
for i = 1:nbmiss
    acc = psnr(X_hat(:,:,i),gt(:,:,i),pval);
    disp([num2str(i) 'st frame, PSNR = ' num2str(acc)]);
    subplot(3,nbmiss,2*nbmiss+i);
    imshow(X_hat(:,:,i));
end

end

