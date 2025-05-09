clear all; 
clc;

%% 参数设置
p = 2;                          % 信号源数量，表示有两个信号源
M = 5;                          % 阵元数，表示天线阵列的元素数目
fc = 1e9;                        % 载波频率，1 GHz（10^9 Hz）
DOA_true = sort([-20,40]);       % 真实DOA角度（度），这里是-20°和40°，并进行排序
fs = 4*fc;                       % 采样率，设定为载波频率的4倍
N = 64;                          % 快拍数，表示每次实验中采集的样本数

snr_values = -10:5:20;           % 信噪比的范围，-10到20dB，步长为5
c = 3e8;                         % 光速，3×10^8 m/s
d = 0.15;                        % 阵元间距，假设为15cm
lambda = c/fc;                   % 波长，光速除以载波频率

k = 2*pi/lambda;                 % 波数

theta_scan = -90:0.1:90;         % 角度扫描范围，从-90°到90°，步长为0.1°
num_trials = 100;                % 蒙特卡洛实验次数，进行100次实验来平均化结果
RMSE = zeros(length(snr_values), 1); % 初始化存储每个信噪比值下的平均RMSE（均方根误差）

%% 迭代不同信噪比
for idx = 1:length(snr_values)   % 循环遍历每一个信噪比值
    snr = snr_values(idx);       % 当前的信噪比
    temp_rmse = zeros(num_trials, 1); % 临时存储每次实验的RMSE，初始化为零
    
    for trial = 1:num_trials     % 进行num_trials次蒙特卡洛实验
        %% 生成接收信号
        A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(DOA_true) / lambda);  % 计算阵列响应矩阵，基于真实的DOA角度
        S = sqrt(2)*(randn(p, N) + 1j*randn(p, N)); % 生成复高斯白噪声信号源，p为信号源数目，N为快拍数
        X = awgn(A * S, snr, 'measured');  % 将信号添加白噪声（使用指定的信噪比）
        
        %% 计算样本协方差矩阵
        R = X * X' / N;  % 计算协方差矩阵，X是接收到的信号，N为快拍数
        [U, D] = eig(R);
        [~, order] = sort(diag(D), 'descend');
        U = U(:, order);
        Un = U(:, p+1:end);
        %% DML算法
        theta_scan = -90:0.5:90; % 减小扫描步长提高精度
        f_dml = zeros(size(theta_scan));
        for i = 1:length(theta_scan)
            a = exp(-1j*k*d*(0:M-1)'*sind(theta_scan(i)));
            P_A = a*pinv(a);
            f_dml(i) = real(trace(P_A*R));
        end
        [~, locs] = findpeaks(f_dml, 'SortStr','descend','NPeaks',p);
        est_DOA_dml = sort(theta_scan(locs));
        if length(est_DOA_dml) < p
              random_values = -90 + (90 + 90) * rand(1, p - length(est_DOA_dml));
              est_DOA_dml = [est_DOA_dml, random_values];
        end
        est_DOA_dml=sort(est_DOA_dml);
                     
        %% 计算RMSE（所有角度已排序，直接匹配）
        rmse_dml = sqrt(mean((DOA_true - est_DOA_dml).^2));  % 计算估计DOA与真实DOA之间的均方根误差
  
        temp_rmse(trial, :) = [rmse_dml];  % 存储每次实验的RMSE值
    end
    
    %% 计算平均RMSE
    RMSE(idx, :) = mean(temp_rmse, 1);  % 对所有实验的RMSE值进行平均，得到每个信噪比下的平均RMSE
end

%% 绘图
figure;                   % 创建新图形窗口
hold on;                  % 保持当前图形，不覆盖
plot(snr_values, RMSE(:, 1), '-o', 'DisplayName', 'DML', 'LineWidth', 1.5);  % 绘制信噪比与RMSE的关系曲线，使用圆点标记
xlabel('SNR (dB)');       % x轴标签：信噪比（dB）
ylabel('RMSE');           % y轴标签：均方根误差（RMSE）
grid on;                  % 打开网格
legend;                   % 显示图例
title('不同信噪比下DML方法RMSE评估');  % 图表标题：DOA估计的RMSE
xlim([-10, 20]);          % 设置x轴的范围为-10到20dB
box on;                   % 开启坐标框
