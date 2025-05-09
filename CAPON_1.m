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

R=X*X'/N;
[U, D] = eig(R);  
eig_values = sort(diag(D), 'descend');
Un = U(:, 1:M-p); 

theta_scan = -90:0.1:90;
%% Capon算法
R_inv = inv(R);
P_capon = zeros(size(theta_scan));
for i = 1:length(theta_scan)
  a = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_scan(i)) / lambda);
  P_capon(i) = 1 / real(a' * R_inv * a);
end

P_capon = abs(P_capon);
P_capon = 10 * log10(P_capon / max(P_capon)); 

plot(theta_scan, P_capon, 'LineWidth', 1.5, 'DisplayName', 'CAPON');
xlabel('角度');
ylabel('空间谱 (dB)');
title('CAPON空间谱');
xline(-20, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % 在 x=0 处画黑色虚线
xline(40, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % 在 x=0 处画黑色虚线
grid on;
hold on;
legend;