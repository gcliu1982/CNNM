function [] = demo_1Dtimeseries()
%DEMO1D 此处显示有关此函数的摘要
%   此处显示详细说明
close all;
%% data
x1 = [83.3 83.5 83.2 82.6 82.2 82.1 81.7 82.2 81.6 82.1 82.7 82.8 81.5 82.2 82.3 82.1 83.6 82.7 82.5 81.5 82.1 82.2 82.6 83.3 83.1 83.3 83.7 82.9 82.3 81.8 81.6 80.9 81 81.3 81.4 80.2 80 80.85 80.83 81.1 80.7 81.1 80.83 80.82 81.5 81.6 81.5 81.6 81.8 81.1 80.5 80 80.7 81.3 80.7 80 81.1 81.87 81.91 81.3 81 80.5 80.6 79.8 79.6 78.49 78.49 79.6 80.6 82.3 81.2 79.1 78.6 78.7 78 78.6 78.7 78.6 79.7 80 79.3 79 80.2 81.5 80.8 81 80.96 81.1 80.8 79.7 80 81.6 82.7 82.1 81.7 81.5];
x2 = [101 82 67 35 31 7 20 93 154 126 85 68 39 23 10 24 83 132 131 118 90 67 60 47 41 21 16 6 4 7 15 34 45 43 48 42 28 10 8 3 0 1 5 12 14 35 46 41 30 24 16 7 4 2 9 17 36 50 64 67 71 48 28 9 13 57 122 138 103 86 65 37 24 11 15 40 62 99 125 96 67 65 54 39 21 7 4 23 55 94 96 77 59 44 47 31 16 7 38 74];

%%
figure;
for i = 1:2
    if i == 1
        x0 = x1';
    elseif i == 2
        x0 = x2';
    end
    
    pval = max(abs(x0));
    m = length(x0);
    t = (1:m)';
    disp(['sequence length = ' num2str(m)]);
    
    mb = round(m*0.1); %select 10% as the unknown future data
    miss = m-mb+1:m;
    omega = ones(m,1);
    omega(miss) = 0;
    inds = omega>0.1;
    
    t1 = t(inds);
    y1 = x0(inds);

    %show data
    subplot(2,3,(i-1)*3+1);
    plot(t1,y1,'-ko','MarkerSize',3);
    hold on;
    t1 = t(~inds);
    y1 = x0(~inds);
    plot(t1,y1,'-g*','MarkerSize',3);
    hold off;
    title('future entries 10%');
    if i == 2
        xlabel('data');
    end
    y = x0.*omega;
    
    %run DFT_L1
    x = inexact_alm_dft_nD(y,omega);
    subplot(2,3,(i-1)*3+2);
    t1 = t(inds);
    y1 = x(inds);
    plot(t1,y1,'-ko','MarkerSize',3);
    hold on;
    t1 = t(~inds);
    y1 = x(~inds);
    plot(t1,y1,'-r*','MarkerSize',3);
    hold off;
    xlabel('DFT_{L1}');
    acc = psnr(x(~inds),x0(~inds),pval);
    err = rmse(x(~inds),x0(~inds));
    disp(['DFT_L1(PSNR) = ' num2str(acc)]);
    disp(['DFT_L1(RMSE) = ' num2str(err)]);

    %run CNNM
    if i==1
        ksize = 74;
    elseif i==2
        ksize = 56;
    end
    x = inexact_alm_cnm_1D(y,omega,ksize);
    subplot(2,3,(i-1)*3+3);
    t1 = t(inds);
    y1 = x(inds);
    plot(t1,y1,'-ko','MarkerSize',3);
    hold on;
    t1 = t(~inds);
    y1 = x(~inds);
    plot(t1,y1,'-r*','MarkerSize',3);
    hold off;
    xlabel('CNNM');
    acc = psnr(x(~inds),x0(~inds),pval);
    err = rmse(x(~inds),x0(~inds));
    disp(['CNNM(PSNR) = ' num2str(acc)]);
    disp(['CNNM(RMSE) = ' num2str(err)]);
end


end

