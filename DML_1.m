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
k = 2*pi/lambda;        % 波数
% 生成接收信号
A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(DOA_true) / lambda);
S = sqrt(2)*(randn(p, N) + 1j*randn(p, N));
X = awgn(A * S, snr, 'measured');

R=X*X'/N;
[U, D] = eig(R);  
eig_values = sort(diag(D), 'descend');
Un = U(:, 1:M-p); 

theta_scan = -90:0.1:90;
%% DML算法

f_dml = zeros(size(theta_scan));
for i = 1:length(theta_scan)
   a = exp(-1j*k*d*(0:M-1)'*sind(theta_scan(i)));
   P_A = a*pinv(a);
   f_dml(i) = real(trace(P_A*R));
end

f_dml = abs(f_dml);
f_dml = 10 * log10(f_dml / max(f_dml)); 

plot(theta_scan, f_dml, 'LineWidth', 1.5, 'DisplayName', 'DML');
xlabel('角度');
ylabel('空间谱 (dB)');
title('DML空间谱');
xline(-20, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % 在 x=0 处画黑色虚线
xline(40, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % 在 x=0 处画黑色虚线
grid on;
hold on;
legend;