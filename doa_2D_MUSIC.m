clear all;
clc;

%% 参数设置
p = 2;                          % 信源数
Mx = 8; My = 8;                 % x和y方向阵元数
fc = 1e9;                       % 载频
DOA_az = [30, -20];             % 方位角（度）
DOA_el = [45, 30];              % 俯仰角（度）
fs = 4*fc;                      % 采样率
N = 1024;                       % 快拍数
snr = [15, 10];                 % 信噪比(dB)
c = 3e8;                        % 光速
d = 0.15;                       % 阵元间距
lambda = c/fc;                  % 波长

%% 信号模型
t = (0:N-1)/fs;
s1 = sqrt(2)*cos(2*pi*fc*t);    % 信号1
s2 = sqrt(2)*cos(2*pi*(fc+5e6)*t); % 信号2
S = [s1; s2];                   % 信号矩阵

%% 二维阵列流形矩阵构建
A = zeros(Mx*My, p);
for k = 1:p
    % x方向相位差
    phase_x = 2*pi*d/lambda * sind(DOA_el(k)) * cosd(DOA_az(k)) * (0:Mx-1);
    a_x = exp(-1j*phase_x(:));
    
    % y方向相位差
    phase_y = 2*pi*d/lambda * sind(DOA_el(k)) * sind(DOA_az(k)) * (0:My-1);
    a_y = exp(-1j*phase_y(:));
    
    % Kronecker积得到二维导向矢量
    A(:,k) = kron(a_x, a_y);
end

%% 接收信号生成
X = awgn(A*S, snr(1));          % 添加高斯白噪声

%% 协方差矩阵计算
R = X*X'/N;                     % 样本协方差矩阵
[U, D] = eig(R);                % 特征分解
eig_values = diag(D);
[~, idx] = sort(eig_values, 'descend');
Un = U(:, idx(p+1:end));        % 噪声子空间

%% 二维空间谱扫描
az_scan = -90:1:90;             % 方位角扫描范围
el_scan = 0:1:90;               % 俯仰角扫描范围
P_music = zeros(length(el_scan), length(az_scan));

for i = 1:length(el_scan)
    for j = 1:length(az_scan)
        % 当前扫描角度对应的导向矢量
        phase_x = 2*pi*d/lambda * sind(el_scan(i)) * cosd(az_scan(j)) * (0:Mx-1);
        a_x = exp(-1j*phase_x(:));
        
        phase_y = 2*pi*d/lambda * sind(el_scan(i)) * sind(az_scan(j)) * (0:My-1);
        a_y = exp(-1j*phase_y(:));
        
        a = kron(a_x, a_y);
        
        % MUSIC谱计算
        P_music(i,j) = 1/abs(a'*(Un*Un')*a);
    end
end

%% 结果可视化
P_music = 10*log10(P_music/max(P_music(:))); % 归一化为dB

% 三维谱峰图
figure;
mesh(az_scan, el_scan, P_music);
xlabel('方位角(°)'); ylabel('俯仰角(°)'); zlabel('空间谱(dB)');
title('二维MUSIC空间谱');
hold on;

% 标记真实DOA
for k = 1:p
    plot3(DOA_az(k), DOA_el(k), 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
legend('空间谱', '真实DOA');
colorbar;

%% 峰值搜索（自动提取前p个峰值）
[el_idx, az_idx] = findpeaks2D(P_music, p);
estimated_az = az_scan(az_idx);
estimated_el = el_scan(el_idx);
disp('估计结果:');
disp(['方位角: ', num2str(estimated_az)]);
disp(['俯仰角: ', num2str(estimated_el)]);

%% 自定义二维峰值查找函数
function [row_idx, col_idx] = findpeaks2D(matrix, num_peaks)
    [rows, cols] = size(matrix);
    peak_values = zeros(1, num_peaks);
    row_idx = zeros(1, num_peaks);
    col_idx = zeros(1, num_peaks);
    
    for k = 1:num_peaks
        [max_val, max_pos] = max(matrix(:));
        [r, c] = ind2sub(size(matrix), max_pos);
        
        peak_values(k) = max_val;
        row_idx(k) = r;
        col_idx(k) = c;
        
        % 将当前峰值周围归零避免重复检测
        r_start = max(1, r-3); r_end = min(rows, r+3);
        c_start = max(1, c-3); c_end = min(cols, c+3);
        matrix(r_start:r_end, c_start:c_end) = -inf;
    end
end