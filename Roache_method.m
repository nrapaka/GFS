% Roache's Method for First Derivative via FFT with Reduction to Periodicity

% Parameters
M = 128;             % Number of points
x = linspace(0,1,M); % Domain
L = 1;               % Domain length
dx = x(2) - x(1);
k = [0:M/2 -M/2+1:-1] * (2*pi/L); % FFT wavenumbers

% Test function: f(x) = exp(-2x) * sin(2*pi*x)
f = exp(-2*x) .* sin(2*pi*x);

% Exact derivative for validation
f_exact = -2*exp(-2*x).*sin(2*pi*x) + 2*pi*exp(-2*x).*cos(2*pi*x);

% Roache's Method for First Derivative via FFT with Reduction to Periodicity

% Parameters
M = 128;             % Number of points
x = linspace(0,1,M); % Domain
L = 1;               % Domain length
dx = x(2) - x(1);
k = [0:M/2 -M/2+1:-1] * (2*pi/L); % FFT wavenumbers

% Degree of reducing polynomial (user-defined)
N = 0; % Change this as needed

% Test function: f(x) = exp(-2x) * sin(2*pi*x)
f = exp(-2*x) .* sin(2*pi*x);

% Exact derivative for validation
f_exact = -2*exp(-2*x).*sin(2*pi*x) + 2*pi*exp(-2*x).*cos(2*pi*x);

% Compute boundary derivatives analytically (up to N)
uj = zeros(1,N+1); % uj0 to ujN
syms xs
f_sym = exp(-2*xs) * sin(2*pi*xs);
for n = 0:N
    uj(n+1) = double(subs(diff(f_sym, xs, n), 1)) - double(subs(diff(f_sym, xs, n), 0));
end

% Compute polynomial coefficients for g(x)
a = zeros(1,N+1);
a(N+1) = uj(N+1) / factorial(N); % Highest order coefficient
for n = N-1:-1:1
    coeff_sum = 0;
    for j = n+1:N
        coeff_sum = coeff_sum + a(j+1) * nchoosek(j, n);
    end
    a(n+1) = (uj(n+1) - coeff_sum) / n;
end
a(1) = double(subs(f_sym, 0)); % Constant term

% Evaluate g(x) and its derivative analytically
g = zeros(size(x));
gp = zeros(size(x));
for i = 1:N+1
    g = g + a(i) * x.^(i-1);
    if i > 1
        gp = gp + a(i)*(i-1)*x.^(i-2);
    end
end

% Residual function
f_res = f - g;

% FFT-based derivative of the residual
fhat = fft(f_res);
df_res = real(ifft(1i*k .* fhat));

% Total derivative
df_total = df_res + gp;

% Plot
figure;
plot(x, f_exact, 'k', 'LineWidth', 1.2); hold on;
plot(x, df_total, 'r--');
legend('Exact df/dx', 'Roache''s Method');
title(['Roache''s Method for First Derivative (N = ', num2str(N), ')']);
xlabel('x'); ylabel('df/dx');
grid on;

% Plot
figure;
plot(x, f_exact, 'k', 'LineWidth', 1.2); hold on;
plot(x, df_total, 'r--');
legend('Exact df/dx', 'Roache''s Method');
title('Roache''s Method for First Derivative');
xlabel('x'); ylabel('df/dx');
grid on;
