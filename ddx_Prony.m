function [ux_prony, u_prony] = ddx_Prony(x,u,M,ux_ex)
    %% === Prony Method (Based on Plonka's Book) ===
    % M = 2; % Number of exponentials assumed known
    N = 2*M;
    dx = x(2)-x(1);
    a = x(1);
%     x = x-a;
    
    % Step 1: Construct Hankel matrix and RHS
    H = zeros(M, M);
    rhs = zeros(M, 1);
    for i = 1:M
        H(i, :) = u(i:i+M-1);
        rhs(i) = -u(M+i);
    end
    
    % Step 2: Solve for Prony polynomial coefficients
    p = H \ rhs;
    poly_coeffs = [1; flip(p)];
    % the polynomial is C(1)*X^N + ... + C(N)*X + C(N+1).
    z = roots(poly_coeffs);     % Roots (z_j = e^{phi_j})
    
    % Step 3: Compute phi_j
    r = z./abs(z);
    phi = (log(abs(z)) + 1i * imag(log(r)));
    phi = phi / dx;
%     phi = log(z) / dx;           % phi_j = ln(z_j) / dx
    
    % Step 4: Solve Vandermonde system to find coefficients c_j
    V = zeros(M, M);
    for i = 1:M
        V(i, :) = z.'.^(i-1);
    end
    c = V \ u(1:M);            % Coefficients c_j
%     c = c .* exp(-phi*x(1));
    
%     kc =  real(phi)
%     ks =  imag(phi)
%     uc =  2*real(c)
%     us = -2*imag(c)
    
    % Step 5: Reconstruct derivative using h'(x) = sum c_j * phi_j * e^{phi_j x}
     u_prony = zeros(size(x));
    ux_prony = zeros(size(x));
    for j = 1:M
         u_prony =  u_prony + c(j) * exp(phi(j) * x);
        ux_prony = ux_prony + c(j) * phi(j) * exp(phi(j) * x);
    end

%     x = x+a;
%     u_prony
    % Plot results
    figure(226); clf;
    subplot(2,1,1);
    plot(x(1:end-1), u(1:end-1), 'k-', 'LineWidth', 2); hold on;
    plot(x(1:end-1), u_prony(1:end-1), 'b-o', 'LineWidth', 2);
    legend('original', 'Prony');
    axis tight;
    xlabel('x');
    ylabel('u(x)');
    title('Approximation with Prony method');
    grid on;
    % return;
    subplot(2,1,2);  hold on;
    plot(x(1:end-1), ux_ex(1:end-1), 'k-', 'LineWidth', 2);
    plot(x(1:end-1), ux_prony(1:end-1), 'm--', 'LineWidth', 2); 
    legend('analytical', 'Prony');
    axis tight;
    xlabel('x');
    ylabel("u'(x)");
    title('Derivative computed with Prony''s method');
    grid on;
    
end

