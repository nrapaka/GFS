function [df, df_smooth, df_singular, f_smooth, f_singular] = ddx_Eckhoff(x, f, fj, Q, M, df_ex)
% ECKHOFF_DERIVATIVE Compute first derivative using Eckhoff's method
%   [df, gamma, A] = eckhoff_derivative(f, N, Q, M, gamma_guess)
%
% Inputs:
%   f - function handle for the function to differentiate
%   N - number of Fourier modes (even)
%   Q - order of reconstruction (typically 1-3)
%   M - number of discontinuities (default 2 for endpoints)
%   gamma_guess - initial guess for discontinuity locations (default [-pi, pi])
%
% Outputs:
%   df - derivative approximation at grid points
%   gamma - detected discontinuity locations
%   A - jump amplitudes at discontinuities (A = f(x+) - f(x-); @boundaries: A = f(1) - f(end))

N = size(x);
dx = x(2) - x(1);   % mesh spacing
L = N*dx; % domain length

% Set defaults
if isempty(M)
    M = 1; % Default to 1 discontinuity at endpoints
end

gamma = x(1); %linspace(-pi, pi, M+1);
% gamma_guess = x(end)+dx; %linspace(-pi, pi, M+1);

A = reshape(-fj,M,Q+1);

% % Create k vector
% k = ifftshift(-N/2:N/2-1)';
% 
% % Compute Fourier coefficients
% f_hat = fft(ifftshift(f));

% % Step 1: Find discontinuity locations (gamma)
% gamma = find_discontinuities(f_hat, k, M, gamma);

% % Step 2: Compute jump amplitudes (A)
% A = compute_jump_amplitudes(f_hat, k, gamma, Q, N, f, x);

% Step 3: Reconstruct the smooth part
[f_smooth,f_singular] = reconstruct_smooth_part(f, x, gamma, A, Q);
% [f_smooth,~] = reconstruct_smooth_part(f, x, gamma, A, Q);
% [f_smooth_hat,f_discontinuous_hat] = reconstruct_smooth_part_spectrum(f_hat, k, gamma, A, Q);
% f_smooth = fftshift(ifft(f_smooth_hat));
% f_singular = fftshift(ifft(f_discontinuous_hat));


% Step 4: Compute derivative
[df, df_smooth, df_singular] = compute_derivative(f_smooth, gamma, A, Q, x);

% add the last element for periodic part for plotting compatibility
f_smooth = add_last_element(f_smooth);
df_smooth = add_last_element(df_smooth);


% M = floor(N/2)*2;
% kn = 2*pi/L * ifftshift(-M/2:M/2-1).'; % modes
% % kn = 2*pi/L * ifftshift(-M/2:M/2).'; % modes
% % kn(1) = 0;
% % use this for even derivatives
% % kn = 2*pi/L * [0:M/2-1,-M/2,-M/2+1:-1].'; 
% % use this for odd derivatives: nullify odd-ball wave number
% % kn = 2*pi/L * [0:M/2-1, 0  ,-M/2+1:-1].'; 
% 
% uk = fft(ifftshift(u(1:M)));   % FFT: ifftshift rearranges x=-L/2:L/2-dx to x=0:L-dx
% % uk(M/2+1)=0;    % nullify odd-ball wave number
% ux = fftshift(ifft(1i*kn.*uk));  % iFFT
% ux = [ux;ux(1)];
% ux = real(ux);
% uF = fftshift(ifft(uk));  % iFFT
% uF = [uF;uF(1)];



% Plot results
figure(225); clf;
subplot(2,1,1);
plot(x, f, 'k-', 'LineWidth', 2); hold on;
plot(x, f_smooth, 'b-', 'LineWidth', 2);
plot(x, f_singular, 'r--', 'LineWidth', 2);
% plot(x, fs, 'bo', 'LineWidth', 2);
legend('original', 'smooth', 'discontinuous')
axis tight;
xlabel('x');
ylabel('u(x)');
title('Decomposition with Eckhoff method');
grid on;
% return;
subplot(2,1,2);  hold on;
plot(x, df_ex, 'k-', 'LineWidth', 2);
plot(x, df, 'm--', 'LineWidth', 2); 
plot(x, df_smooth, 'b-', 'LineWidth', 2); 
plot(x, df_singular, 'r--', 'LineWidth', 2); 
legend('analytical', 'Echoff', 'du-smooth','du-singular');
axis tight;
xlabel('x');
ylabel("u'(x)");
title('Derivative computed with Eckhoff''s method');
grid on;


end

function [f_smooth,f_singular] = reconstruct_smooth_part(f, x, gamma, A, Q)
% RECONSTRUCT_SMOOTH_PART Compute Fourier coefficients of smooth part
    
    M = size(A, 1);
    xx = x;
    xx(end) = xx(end) - 1e-14;
    % Compute derivative of singular part
    f_singular = zeros(size(x));
    for j = 1:M
        for n = 0:Q
            % The derivative of V_n(x; gamma_j) is V_{n-1}(x; gamma_j)
            f_singular = f_singular + A(j,n+1)*V_function(n, xx - gamma(j));
        end
    end

    % add a constant to match the first node of singular part with the
    % original part (does not affect derivatives)
    f_singular = f_singular + f(1)-f_singular(1);
    
    f_smooth = f - f_singular;
    f_smooth = remove_last_element(f_smooth);
end

function f_periodic = add_last_element(f)
    f_periodic = [f;f(1)];
end

function f_periodic = remove_last_element(f)
    f_periodic = f(1:end-1);
end

function [df,df_smooth,df_singular] = compute_derivative(f_smooth, gamma, A, Q, x)
% COMPUTE_DERIVATIVE Compute the derivative using the reconstruction
    
    N = length(f_smooth);
    N = floor(N/2)*2;
    M = size(A, 1);
    
    % Create k vector
    k = ifftshift(-N/2:N/2-1)';

    f_smooth_hat = fft(ifftshift(f_smooth));
    
    % Compute derivative of smooth part via Fourier differentiation
    df_smooth = fftshift(ifft(1i*k.*(f_smooth_hat)));
    
    xx = x;
    xx(end) = xx(end) - 1e-14;
    
    % Compute derivative of singular part
    df_singular = zeros(size(x));
    for j = 1:M
        for n = 0:Q
            % The derivative of V_n(x; gamma_j) is V_{n-1}(x; gamma_j)
%             if n > 0
                df_singular = df_singular + A(j,n+1)*V_function(n-1, xx - gamma(j));
%             end
        end
    end
    
    % Total derivative
    df = real(add_last_element(df_smooth) + df_singular);
end

function v = V_function(n, xi)
% V_FUNCTION The special function V_n(xi) from Eckhoff's paper
% This is the periodic extension of the Bernoulli polynomial
    
    % Bring xi from [-π,π) into [0, 2π)
    xi = mod(xi + 2*pi, 2*pi);  % [-π,π) -> [π,2π) + [0,π)
    % xi = mod(xi + pi, 2*pi) - pi;
%     xi/(2*pi)
    
    v = - (2*pi)^n/factorial(n+1) * bernoulli(n+1,xi/(2*pi));
    
%     % Bernoulli polynomial implementation (simplified)
%     % For a complete implementation, we'd use the exact Bernoulli polynomials
%     switch n
%         case 0
%             % V_0(x) is a sawtooth function
%             v = -xi/(2*pi) - 0.5;
%             v(xi > 0) = v(xi > 0) + 1;
%         case 1
%             % V_1(x) is piecewise quadratic
%             v = -xi.^2/(4*pi) + xi/2 - pi/6;
%             v(xi > 0) = v(xi > 0) + xi(xi > 0) - pi;
%         otherwise
%             % For higher orders, we'd implement the full Bernoulli polynomials
%             % This is a placeholder - would need proper implementation
%             v = zeros(size(xi));
%             warning('Higher order V_n functions not fully implemented');
%     end
% 
end

function gamma = find_discontinuities(f_hat, k, M, gamma_guess)
% FIND_DISCONTINUITIES Locate discontinuities using Eckhoff's method
    
    % Use Prony's method to find roots of characteristic polynomial
    % This is a simplified version - see Eckhoff's paper for full details
    
    % For M=2 (endpoints), we know the locations
    if M == 1
        gamma = [-pi];
        % gamma = [-1];
        return;
    end

    if M == 2
        gamma = [-pi; pi];
        % gamma = [-1; 1];
        return;
    end
    
    % For general case, we'd implement the full algorithm from the paper
    % Here we just return the initial guess as a placeholder
    gamma = gamma_guess(:);
    
    % Note: The full implementation would solve the nonlinear system
    % described in Eckhoff's paper to find the singularity locations
end

function A = compute_jump_amplitudes(f_hat, k, gamma, Q, N, f, x)
% COMPUTE_JUMP_AMPLITUDES Calculate jump amplitudes at discontinuities
    
    M = length(gamma);
    off = 8;    % additional modes for least squares

    % Construct the system matrix
    L = floor(N/2) - M*(Q+1) - off;
    U = floor(N/2) - 1;
    n_eq = U - L + 1;
    
    % Preallocate
    lhs = zeros(n_eq, M*(Q+1));
    Ctilde = zeros(n_eq, 1);
    
    % Create the system of equations
    idx = 1;
    for m = L:U
        for j = 1:M
            for n = 0:Q
                lhs(idx, (j-1)*(Q+1)+n+1) = exp(-1i*k(m)*gamma(j))/(1i*k(m))^(n);
            end
        end
        % f_hat from FFT needs to be divided by N to get the actual Fourier coefficients
        Ctilde(idx) = (2*pi*1i*k(m)) * (f_hat(m)/N);  
        idx = idx + 1;
    end
    
    % Solve the system (using least squares if overdetermined)
    A_vec = lhs \ Ctilde;
%     lhs
%     Ctilde
%     k
%     f_hat
    
    % Reshape into M x (Q+1) matrix
    A = reshape(A_vec, Q+1, M)';
    A
    eps = 1e-8;
    Jvec = (f(x(1))-f(x(end))); %f(gamma+eps)-f(gamma-eps);
%     Jvec = 0.5*(f(-1)-f(1)); %f(gamma+eps)-f(gamma-eps);
    A = repmat(Jvec,M,1)
end

function [f_smooth_hat,f_discontinuous_hat] = reconstruct_smooth_part_spectrum(f_hat, k, gamma, A, Q)
% RECONSTRUCT_SMOOTH_PART Compute Fourier coefficients of smooth part
    
    M = size(A, 1);
    N = length(f_hat);
    f_smooth_hat = f_hat;
    f_aperiodic = zeros(size(f_hat));
    f_discontinuous_hat = zeros(size(f_hat));
    
    % Subtract the singular part from the Fourier coefficients
    for m = 1:length(k)
        if k(m) ~= 0
            sum_term = 0;
            for j = 1:M
                for n = 0:Q
                    sum_term = sum_term + A(j,n+1)*exp(-1i*k(m)*gamma(j))/(1i*k(m))^(n+1);
                end
            end
            f_discontinuous_hat(m) = sum_term/(2*pi)*N;
            f_smooth_hat(m) = f_hat(m) - sum_term/(2*pi)*N;
        end
    end
end

