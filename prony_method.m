clc; clear; 

%% Signal Parameters
Fs = 64;       % Sampling frequency (Hz)
T = 1 / Fs;      % Sampling period (seconds)
L = 64;        % Number of samples
% t = (0:L-1) * T * 2 * pi - pi;  % Time vector
t = (0:L-1) * T ;  % Time vector
T

k1 = 5.3;
k2 = 12.4;
% k1 = 1;
% k2 = 2;
% Signal: sum of two non-harmonic sinusoids
x = 0.7 * sin(k1 * t) + sin(k2 * t);
dx_true = 0.7 * k1 * cos(k1 * t) + k2 * cos(k2 * t); % True derivative

% x = 0.7 * exp(1i*k1 * t) + exp(1i*k2 * t);
% dx_true = 0.7 * 1i*k1 * exp(1i*k1 * t) + 1i*k2 * exp(1i*k2 * t); % True derivative

%% === Prony Method (Based on Plonka's Book) ===
M = 32; % Number of exponentials assumed known
N = 2*M;

% Step 1: Construct Hankel matrix and RHS
H = zeros(M, M);
rhs = zeros(M, 1);
for i = 1:M
    H(i, :) = x(i:i+M-1);
    rhs(i) = -x(M+i);
end

% Step 2: Solve for Prony polynomial coefficients
p = H \ rhs;
poly_coeffs = [1; flip(p)];
% the polynomial is C(1)*X^N + ... + C(N)*X + C(N+1).
z = roots(poly_coeffs);     % Roots (z_j = e^{phi_j})
z

% Step 3: Compute phi_j
r = z./abs(z);
phi = (log(abs(z)) + 1i * imag(log(r)));
phi = phi / T;
% phi = log(z) / T;           % phi_j = ln(z_j) / T
phi

% Step 4: Solve Vandermonde system to find coefficients c_j
V = zeros(M, M);
for i = 1:M
    V(i, :) = z.'.^(i-1);
end
c = V \ x(1:M).';            % Coefficients c_j
c

kc =  real(phi)
ks =  imag(phi)
uc =  2*real(c)
us = -2*imag(c)

% Step 5: Reconstruct derivative using h'(t) = sum c_j * phi_j * e^{phi_j t}
 x_prony = zeros(size(t));
dx_prony = zeros(size(t));
for j = 1:M
     x_prony =  x_prony + c(j) * exp(phi(j) * t);
    dx_prony = dx_prony + c(j) * phi(j) * exp(phi(j) * t);
end

%% === GFS-style Cosine-Sine Fit (Non-harmonic) ===
freqs = [k1, k2];
Phi = [cos(freqs(1)*t') sin(freqs(1)*t') cos(freqs(2)*t') sin(freqs(2)*t')];
A_gfs = Phi \ x';

dx_gfs = -freqs(1) * A_gfs(1) * sin(freqs(1) * t) + ...
          freqs(1) * A_gfs(2) * cos(freqs(1) * t) - ...
          freqs(2) * A_gfs(3) * sin(freqs(2) * t) + ...
          freqs(2) * A_gfs(4) * cos(freqs(2) * t);

%% Compute Fourier Transform using fft()
Y = fft(x);

% Compute two-sided spectrum and convert to single-sided spectrum
P2 = abs(Y / L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1); % Compensate for power on both sides

% Frequency axis for single-sided spectrum
f = Fs * (0:(L/2)) / L;
% f = linspace(0,Fs/2,Fs/2+1);

%% === Plot Comparison of Derivatives ===
figure(2); clf;
plot(t, dx_true, 'k', 'LineWidth', 2); hold on;
plot(t, real(dx_prony), 'r--', 'LineWidth', 2);
plot(t, dx_gfs, 'o', 'LineWidth', 2);
legend('True Derivative', 'Prony (Plonka)', 'GFS', 'Location', 'best');
xlabel('Time t'); ylabel('f''(t)');
title('Comparison of Derivative Methods');
grid on; set(gca, 'FontSize', 14);
ylim([-20,20]);

%% === Error Analysis ===
err_prony = norm(dx_true - real(dx_prony)) / norm(dx_true);
err_gfs = norm(dx_true - dx_gfs) / norm(dx_true);
fprintf('Relative Error (Prony): %.2e\n', err_prony);
fprintf('Relative Error (GFS):   %.2e\n', err_gfs);

fs = 15;
lw = 2;
% Plot the original signal
figure(4); clf;
subplot(2, 1, 1);
h = plot(t, x,'-','LineWidth',lw); hold on;
h = plot(t, x_prony,'o','LineWidth',lw);
% title('Time-domain signal','Interpreter','latex');
% title('$f(t)=0.7\sin(5t)+1.0\sin(12t)$','Interpreter','latex');
title('$f(t)=0.7\sin(5.3t)+1.0\sin(12.4t)$','Interpreter','latex');
xlabel('time, t (s)','Interpreter','latex');
ylabel('$f(t)$','Interpreter','latex');
set(gca,'FontSize',fs);
legend('original','Prony')
axis tight;

grid on;

% Plot the Fourier spectrum
subplot(2, 1, 2);
h = plot(f, P1,'-o','LineWidth',lw);
title('Single-sided amplitude spectrum of $f(t)$','Interpreter','latex');
xlabel('frequency, $\omega (Hz)$','Interpreter','latex');
ylabel('$|\hat{f}(\omega)|$','Interpreter','latex');
set(gca,'FontSize',fs);
axis tight;
grid on;
% saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/f_periodic.jpg');
% saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/f_aperiodic.jpg');

