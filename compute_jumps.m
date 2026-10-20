function jumps = compute_jumps(x_vals, u_vals, max_order, stw)
    % MATLAB Function to Compute Jumps of Derivatives of u(x) at Endpoints
    % Input:
    %   x_vals - Array of x values (including endpoints)
    %   u_vals - Array of corresponding u(x) values
    %   max_order - Maximum order of derivatives to compute
    %   stw - Number of points to use for polynomial fitting at each endpoint
    % Output:
    %   jumps - Array of jumps at the endpoints for each derivative up to max_order

    % Define endpoints for the interval -pi <= x <= pi
    x_start = -pi;
    x_end = pi;

    % Fit polynomials using the first and last 'stw' points
    poly_coeffs_start = polyfit(x_vals(1:stw), u_vals(1:stw), stw - 1);
    poly_coeffs_end = polyfit(x_vals(end-stw+1:end), u_vals(end-stw+1:end), stw - 1);

    % Initialize the array to store jumps
    jumps = zeros(1, max_order);

    % Compute derivatives and evaluate jumps at endpoints
    for n = 1:max_order
        % Compute the nth derivative coefficients for both ends
        derivative_coeffs_start = polyder(poly_coeffs_start, n);
        derivative_coeffs_end = polyder(poly_coeffs_end, n);

        % Evaluate the nth derivative at endpoints
        u_start = polyval(derivative_coeffs_start, x_start);
        u_end = polyval(derivative_coeffs_end, x_end);

        % Compute the jump
        jumps(n) = u_end - u_start;
    end
end

function derivative_coeffs = polyder(coeffs, n)
    % Compute the nth derivative of a polynomial given its coefficients
    derivative_coeffs = coeffs;
    for k = 1:n
        degree = length(derivative_coeffs) - 1;
        derivative_coeffs = derivative_coeffs(1:end-1) .* (degree:-1:1);
    end
end
