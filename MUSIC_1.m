clear all; 
clc;

% 参数设置
p = 2;                  % 信源数目
M = 5;                  % 阵元数目
fc = 1e9;               % 载波频率 (Hz)
DOA_true = [-20, 40];   % 真实DOA角度 (度)
fs = 4 * fc;            % 采样频率 (Hz)
N = 64;                 % 快拍数
snr = 10;               % 信噪比 (dB)
c = 3e8;                % 光速 (m/s)
d = 0.15;               % 阵元间距 (m)
lambda = c / fc;        % 波长 (m)

% 生成接收信号
A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(DOA_true) / lambda);
S = sqrt(2)*(randn(p, N) + 1j*randn(p, N));
X = awgn(A * S, snr, 'measured');

% MUSIC算法实现
R = X*X'/N;
[U, D] = eig(R);  
[D_sort, idx] = sort(diag(D), 'descend');  % 特征值降序排序
U = U(:, idx);                            % 对应特征向量排序

Un = U(:, p+1:end);  % 噪声子空间

theta_scan = -90:0.1:90;
P_music = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_scan(i)) / lambda);
    P_music(i) = 1 / (a' * (Un * Un') * a);
end

% 峰值检测
P_music = abs(P_music);
[peaks, locs] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', p);
estimated_DOA = sort(theta_scan(locs));  % 估计的DOA角度

% 结果可视化
P_music = 10 * log10(P_music / max(P_music)); 
plot(theta_scan, P_music, 'LineWidth', 1.5, 'DisplayName', 'MUSIC');
xlabel('角度 (度)');
ylabel('空间谱 (dB)');
title('MUSIC空间谱');
xline(DOA_true(1), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(DOA_true(2), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
grid on;
legend;

% 显示估计结果
disp('真实DOA:');
disp(DOA_true);
disp('估计DOA:');
disp(estimated_DOA);