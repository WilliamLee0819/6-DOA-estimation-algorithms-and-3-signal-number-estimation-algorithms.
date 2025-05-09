clear all; 
clc;

%% 参数设置
p = 3;                          % 信号源数量
M = 8;                         % 阵元数
fc = 1e9;                       % 载波频率
DOA_true = sort([-40,0, 20]);      % 真实DOA角度(度)并预排序
fs = 4*fc;                      % 采样率
N = 64;                         % 快拍数
snr_values = -20:2:30;          % 信噪比[dB]
c = 3e8;                        % 光速
d = 0.15;                       % 阵元间距
lambda = c/fc;                  % 波长
k = 2*pi/lambda;                % 波数
theta_scan = -90:0.1:90;        % 角度扫描范围
num_trials = 10000;               % 蒙特卡洛实验次数
RMSE = zeros(length(snr_values), 6); % 存储平均RMSE

%% 迭代不同信噪比
for idx = 1:length(snr_values)
    snr = snr_values(idx);
    temp_rmse = zeros(num_trials, 6); % 临时存储每次实验的RMSE
    
    for trial = 1:num_trials
        %% 生成接收信号
        A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(DOA_true) / lambda);
        S = sqrt(2)*(randn(p, N) + 1j*randn(p, N));
        X = awgn(A * S, snr, 'measured');
        
        %% 计算样本协方差矩阵
        R = X * X' / N;
        
        %% MUSIC算法
        [U, D] = eig(R);
        [~, order] = sort(diag(D), 'descend');
        U = U(:, order);
        Un = U(:, p+1:end);
        
        P_music = zeros(size(theta_scan));
        for i = 1:length(theta_scan)
            a = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_scan(i)) / lambda);
            P_music(i) = 1 / abs(a' * (Un * Un') * a);
        end
        [~, peaks_idx] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', p);
        est_DOA_music = sort(theta_scan(peaks_idx)); % 确保估计角度排序
        if length(est_DOA_music) < p
              random_values = -90 + (90 + 90) * rand(1, p - length(est_DOA_music));
              est_DOA_music = [est_DOA_music, random_values];
        end
        est_DOA_music=sort(est_DOA_music);
        %% Root-MUSIC算法
        Gn = Un * Un';
        coe = zeros(1, 2*M-1);

        % 正确构造多项式系数
        for i = -(M-1):(M-1)
         coe(i + M) = sum(diag(Gn, i));  % 注意索引调整为i+M
        end

        % 求根并筛选
        r = roots(coe);
        r = r(abs(r) < 1);  % 放宽筛选条件提高稳定性

        % 选择最接近单位圆的p个根
        [~, I] = sort(abs(abs(r) - 1));

        Theta = r(I(1:p));
        % 转换为角度
        theta_rootmusic = asin(angle(Theta)/pi)/(pi/180);
        theta_rootmusic = sort(theta_rootmusic).';
        if length(theta_rootmusic) < p
              random_values = -90 + (90 + 90) * rand(1, p - length(theta_rootmusic));
              theta_rootmusic = [theta_rootmusic, random_values];
        end
        theta_rootmusic=sort(theta_rootmusic);
        %% LS-ESPRIT算法
        U_s = U(:, 1:p);
        U1 = U_s(1:end-1, :);
        U2 = U_s(2:end, :);
        Phi_ls = (U1' * U1) \ (U1' * U2);
        theta_ls = sort ((asind(-angle(eig(Phi_ls)) * lambda / (2*pi*d))))'; % 排序输出
        
        %% TLS-ESPRIT算法
        C = [U1, U2];
        [~, ~, V] = svd(C);
        V12 = V(1:p, p+1:2*p);
        V22 = V(p+1:2*p, p+1:2*p);
        Phi_tls = -V12 / V22;
        theta_tls = sort ((asind(-angle(eig(Phi_tls)) * lambda / (2*pi*d))))';
        
        %% Capon算法
        R_inv = inv(R);
        P_capon = zeros(size(theta_scan));
        for i = 1:length(theta_scan)
            a = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_scan(i)) / lambda);
            P_capon(i) = 1 / real(a' * R_inv * a);
        end
        [~, peaks_idx] = findpeaks(P_capon, 'SortStr', 'descend', 'NPeaks', p);
        est_DOA_capon = sort(theta_scan(peaks_idx)); % 确保估计角度排序
        if length(est_DOA_capon) < p
              random_values = -90 + (90 + 90) * rand(1, p - length(est_DOA_capon));
              est_DOA_capon = [est_DOA_capon, random_values];
        end
        est_DOA_capon=sort(est_DOA_capon);
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
        rmse_music = sqrt(mean((DOA_true - est_DOA_music).^2));
        rmse_rootmusic = sqrt(mean((DOA_true - theta_rootmusic).^2));
        rmse_ls = sqrt(mean((DOA_true - theta_ls).^2));
        rmse_tls = sqrt(mean((DOA_true - theta_tls).^2));
        rmse_capon = sqrt(mean((DOA_true - est_DOA_capon).^2));
        rmse_dml = sqrt(mean((DOA_true - est_DOA_dml).^2));
        temp_rmse(trial, :) = [rmse_music, rmse_rootmusic, rmse_ls, rmse_tls, rmse_capon,rmse_dml];
    end
    
    %% 计算平均RMSE
    RMSE(idx, :) = mean(temp_rmse, 1);
end

%% 绘图
figure;
hold on;
plot(snr_values, RMSE(:, 1), '-o', 'DisplayName', 'MUSIC', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 2), '-*', 'DisplayName', 'ROOT-MUSIC', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 3), '-s', 'DisplayName', 'LS-ESPRIT', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 4), '-^', 'DisplayName', 'TLS-ESPRIT', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 5), '-d', 'DisplayName', 'Capon', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 6), '->', 'DisplayName', 'DML', 'LineWidth', 1.5);
xlabel('SNR (dB)'); ylabel('RMSE'); grid on;
legend; title(' ');
xlim([-20, 20]);
box on;

axes('Position',[0.55 0.3 0.25 0.25]); % 定义子图位置和大小
box on;  %显示子图的外框线
hold on;
plot(snr_values, RMSE(:, 1), '-o', 'DisplayName', 'MUSIC', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 2), '-*', 'DisplayName', 'ROOT-MUSIC', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 3), '-s', 'DisplayName', 'LS-ESPRIT', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 4), '-^', 'DisplayName', 'TLS-ESPRIT', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 5), '-d', 'DisplayName', 'Capon', 'LineWidth', 1.5);
plot(snr_values, RMSE(:, 6), '->', 'DisplayName', 'DML', 'LineWidth', 1.5);
grid on;

