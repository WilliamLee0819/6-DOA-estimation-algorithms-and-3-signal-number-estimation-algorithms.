clear; clc; close all;

%% 参数设置
fc = 1e9;               % 载波频率 (Hz)
c = 3e8;                % 光速 (m/s)
lambda = c/fc;          % 波长 (m)
d = lambda/2;           % 阵元间距 (半波长)
M = 8;                  % 阵元数量
N = 64;                 % 快拍数
theta_true = [-40,0, 20];  % 真实DOA角度 (度)
K = length(theta_true); % 真实信号源数量
SNR = 10;               % 信噪比 (dB)

%% 生成阵列接收信号
A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_true) / lambda); % 导向矩阵
S = (randn(K, N) + 1j * randn(K, N)) / sqrt(2);                  % 信号源 (复高斯)
noise = (randn(M, N) + 1j * randn(M, N)) * sqrt(10^(-SNR/10)/2); % 噪声
X = A * S + noise;                                               % 接收信号

%% MDL方法估计信号源数量
R_hat = (X * X') / N;                % 样本协方差矩阵
[U, D] = eig(R_hat);                 % 特征分解

eig_values = sort(diag(D), 'descend'); % 降序排列特征值
eigenvalue=eig_values;
%% AIC算法实现
AIC = zeros(1,M);
    for i = 0 : M-1
        sumeigvalue = 0;
        multipeigvalue = 1;
        for j = i+1 : M
            sumeigvalue = sumeigvalue + eigenvalue(j);
            multipeigvalue = multipeigvalue * eigenvalue(j);
        end
        gama = (sumeigvalue * (1/(M-i))) / (multipeigvalue^(1/(M-i)));
        
        AIC(i+1) = 2*N*(M-i)*log(gama) + 2*i*(2*M-i);
    end
    

%% MDL算法实现
MDL = zeros(1,M);
    for i = 0 : M-1
        sumeigvalue = 0;
        multipeigvalue = 1;
        for j = i+1 : M
            sumeigvalue = sumeigvalue + eigenvalue(j);
            multipeigvalue = multipeigvalue * eigenvalue(j);
        end
        gama = (sumeigvalue * (1/(M-i))) / (multipeigvalue^(1/(M-i)));
        
        MDL(i+1) = N*(M-i)*log(gama) + 0.5*i*(2*M-i)*log(N);
    end
    

%% EDC算法实现
EDC = zeros(1,M);
    for i = 0 : M-1
        sumeigvalue = 0;
        multipeigvalue = 1;
        for j = i+1 : M
            sumeigvalue = sumeigvalue + eigenvalue(j);
            multipeigvalue = multipeigvalue * eigenvalue(j);
        end
        gama = (sumeigvalue * (1/(M-i))) / (multipeigvalue^(1/(M-i)));
        
        EDC(i+1) = N*(M-i)*log(gama) + 0.5*i*(2*M-i)*log(log(N));
    end
    


% 找到最小值对应的k值
[~, K_est_AIC] = min(AIC);
K_est_AIC = K_est_AIC - 1; % 转换为0-based索引

[~, K_est_MDL] = min(MDL);
K_est_MDL = K_est_MDL - 1;

[~, K_est_EDC] = min(EDC);
K_est_EDC = K_est_EDC - 1;

%% 可视化
figure('Position', [100, 100, 800, 600]);

% AIC准则
subplot(3,1,1);
plot(0:M-1, AIC, '-o', 'LineWidth', 1.5);
hold on;
plot(K_est_AIC, AIC(K_est_AIC+1), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('候选信号源数量 k');
ylabel('AIC值');
title(['AIC准则']);
grid on;
legend('AIC值', '最小值', 'Location', 'best');

% MDL准则
subplot(3,1,2);
plot(0:M-1, MDL, '-s', 'LineWidth', 1.5);
hold on;
plot(K_est_MDL, MDL(K_est_MDL+1), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('候选信号源数量 k');
ylabel('MDL值');
title(['MDL准则']);
grid on;
legend('MDL值', '最小值', 'Location', 'best');

% EDC准则
subplot(3,1,3);
plot(0:M-1, EDC, '-*', 'LineWidth', 1.5);
hold on;
plot(K_est_EDC, EDC(K_est_EDC+1), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('候选信号源数量 k');
ylabel('EDC值');
title(['EDC准则']);
grid on;
legend('EDC值', '最小值', 'Location', 'best');

%% 显示估计结果
fprintf('真实信号数: %d\n', K);
fprintf('AIC估计信号数: %d\n', K_est_AIC);
fprintf('MDL估计信号数: %d\n', K_est_MDL);
fprintf('EDC估计信号数: %d\n', K_est_EDC);