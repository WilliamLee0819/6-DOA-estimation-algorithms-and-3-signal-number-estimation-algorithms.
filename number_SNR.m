clear; clc; close all;

%% 参数设置
fc = 1e9;               % 载波频率 (Hz)
c = 3e8;                % 光速 (m/s)
lambda = c/fc;          % 波长 (m)
d = lambda/2;           % 阵元间距 (半波长)
M = 8;                  % 阵元数量
N = 64;                 % 快拍数
theta_true = [-40, 0, 20];  % 真实DOA角度 (度)
K_true = length(theta_true); % 真实信号源数量
SNR_range = -20:0.5:20;   % 信噪比范围 (dB)
num_trials = 1000;       % 每个SNR的蒙特卡洛试验次数

%% 初始化存储变量
K_est_MDL = zeros(num_trials, length(SNR_range));
K_est_AIC = zeros(num_trials, length(SNR_range));
K_est_EDC = zeros(num_trials, length(SNR_range));
%% 主循环
for snr_idx = 1:length(SNR_range)
    SNR = SNR_range(snr_idx);
    
    for trial = 1:num_trials
        %% 生成阵列接收信号
        A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_true) / lambda); % 导向矩阵
        S = (randn(K_true, N) + 1j * randn(K_true, N)) / sqrt(2);       % 信号源 (复高斯)
        noise = (randn(M, N) + 1j * randn(M, N)) * sqrt(10^(-SNR/10)/2); % 噪声
        X = A * S + noise;                                               % 接收信号
        
        %% MDL方法估计信号源数量
        R_hat = (X * X') / N;                % 样本协方差矩阵
        [~, D] = eig(R_hat);                 % 特征分解
        eig_values = sort(diag(D), 'descend'); % 降序排列特征值
        
        % 调用MDL函数
        K_est_MDL(trial, snr_idx) = MDL(N, M, eig_values);
        
        % 调用AIC函数
        K_est_AIC(trial, snr_idx) = AIC(N, M, eig_values);

        % 调用EDC函数
        K_est_EDC(trial, snr_idx) = EDC(N, M, eig_values);
      
    end
end

%% 计算RMSE
rmse_mdl = sqrt(mean((K_est_MDL - K_true).^2));
rmse_aic = sqrt(mean((K_est_AIC - K_true).^2));
rmse_edc = sqrt(mean((K_est_EDC - K_true).^2));

%% 绘制成功率曲线
success_rate_mdl = mean(K_est_MDL == K_true, 1);
success_rate_aic = mean(K_est_AIC == K_true, 1);
success_rate_edc = mean(K_est_EDC == K_true, 1);
figure;

plot(SNR_range, success_rate_aic*100, 's-', 'LineWidth', 1.5, 'DisplayName', 'AIC');
hold on;
plot(SNR_range, success_rate_mdl*100, 'o-', 'LineWidth', 1.5, 'DisplayName', 'MDL');

plot(SNR_range, success_rate_edc*100, '*-', 'LineWidth', 1.5, 'DisplayName', 'EDC');
xlabel('SNR (dB)');
ylabel('正确估计概率 (%)');
title(['不同SNR下信号源数量估计方法成功率']);
grid on;
legend;
set(gca, 'FontSize', 12);
ylim([0 100]);

%% AIC算法实现
function [n] = AIC(N,M,eigenvalue)
    aicvalue = zeros(1,M);
    for i = 0 : M-1
        sumeigvalue = 0;
        multipeigvalue = 1;
        for j = i+1 : M
            sumeigvalue = sumeigvalue + eigenvalue(j);
            multipeigvalue = multipeigvalue * eigenvalue(j);
        end
        gama = (sumeigvalue * (1/(M-i))) / (multipeigvalue^(1/(M-i)));
        
        aicvalue(i+1) = 2*N*(M-i)*log(gama) + 2*i*(2*M-i);
    end
    
    [~,index] = min(aicvalue);
    n = index -1;
end
 
%% MDL算法实现
function [n] = MDL(N,M,eigenvalue)
    aicvalue = zeros(1,M);
    for i = 0 : M-1
        sumeigvalue = 0;
        multipeigvalue = 1;
        for j = i+1 : M
            sumeigvalue = sumeigvalue + eigenvalue(j);
            multipeigvalue = multipeigvalue * eigenvalue(j);
        end
        gama = (sumeigvalue * (1/(M-i))) / (multipeigvalue^(1/(M-i)));
        
        aicvalue(i+1) = N*(M-i)*log(gama) + 0.5*i*(2*M-i)*log(N);
    end
    
    [~,index] = min(aicvalue);
    n = index -1;
end


%% EDC算法实现
function [n] = EDC(N,M,eigenvalue)
    aicvalue = zeros(1,M);
    for i = 0 : M-1
        sumeigvalue = 0;
        multipeigvalue = 1;
        for j = i+1 : M
            sumeigvalue = sumeigvalue + eigenvalue(j);
            multipeigvalue = multipeigvalue * eigenvalue(j);
        end
        gama = (sumeigvalue * (1/(M-i))) / (multipeigvalue^(1/(M-i)));
        
        aicvalue(i+1) = N*(M-i)*log(gama) + 0.5*i*(2*M-i)*log(log(N));
    end
    
    [~,index] = min(aicvalue);
    n = index -1;
end
