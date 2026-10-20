function compute_derivative_jump_with_plots()
    %{
    Compute the derivative jump at the domain endpoints using polynomial
    approximation, compare with analytical values, and plot the fitted
    polynomials alongside the analytical function.
    %}
    clc; clear all;
    % Define the function and its derivative
    f = @(x) sin(10.3 * x);             % Highly oscillatory function
    df = @(x) 10.3 * cos(10.3 * x);       % Analytical derivative of f(x)
    
    % Polynomial degree for approximation
    degree = 3;

    N = 10000;
    Nb = degree + 2;
    % Generate data
    x = linspace(0, 2*pi, N);       % Domain
    y = f(x);                         % Function values
        
    % Compute the derivative jump
    jump = compute_derivative_jump(x, y, degree);
    
    % Analytical derivative values at the boundaries
    analytical_left = df(x(1));       % Derivative at the left boundary
    analytical_right = df(x(end));   % Derivative at the right boundary
    analytical_jump = analytical_right - analytical_left;
    
    % Display results
    disp('Computed Results:');
    fprintf('Left Derivative: %.6f\n', jump.left_derivative);
    fprintf('Right Derivative: %.6f\n', jump.right_derivative);
    fprintf('Jump: %.6f\n', jump.value);
    
    disp('Analytical Results:');
    fprintf('Left Derivative: %.6f\n', analytical_left);
    fprintf('Right Derivative: %.6f\n', analytical_right);
    fprintf('Jump: %.6f\n', analytical_jump);

    fprintf('Error in jump: %.6f\n',jump.value-analytical_jump);
    
    % Plot the function and local polynomial fits
    figure(1); clf;
    
    % Plot analytical function
    subplot(2, 1, 1);
    plot(x, y, 'b-', 'DisplayName', 'Analytical Function');
    hold on;
    
    % Fit and plot local polynomial at left endpoint
    num_points = Nb; %degree + 2;
    x_left = x(1:num_points);
    y_left = y(1:num_points);
    p_left = polyfit(x_left, y_left, degree);
    x_fit_left = linspace(x(1), x(num_points), Nb);
    y_fit_left = polyval(p_left, x_fit_left);
    plot(x_fit_left, y_fit_left, 'r--o', 'DisplayName', 'Left Fit');
    
    % Fit and plot local polynomial at right endpoint
    x_right = x(end-num_points+1:end);
    y_right = y(end-num_points+1:end);
    p_right = polyfit(x_right, y_right, degree);
    x_fit_right = linspace(x(end-num_points+1), x(end), Nb);
    y_fit_right = polyval(p_right, x_fit_right);
    plot(x_fit_right, y_fit_right, 'g--o', 'DisplayName', 'Right Fit');
    
    % Formatting
    legend('Location', 'Best');
    title('Function and Local Polynomial Fits');
    xlabel('x');
    ylabel('f(x)');
    grid on;
    
    % Plot derivatives
    subplot(2, 1, 2);
    bar(1, jump.left_derivative, 0.4, 'FaceColor', 'r', 'DisplayName', 'Left Derivative');
    hold on;
    bar(2, jump.right_derivative, 0.4, 'FaceColor', 'g', 'DisplayName', 'Right Derivative');
    bar(3, jump.value, 0.4, 'FaceColor', 'b', 'DisplayName', 'Jump');
    set(gca, 'XTickLabel', {'Left', 'Right', 'Jump'});
    
    % Formatting
    title('Derivative Comparisons');
    ylabel('Derivative Value');
    grid on;
    legend('Location', 'Best');
end

function jump = compute_derivative_jump(x, y, degree)
    %{
    Compute derivative jump at the domain endpoints using polynomial approximation.
    The jump is defined as (derivative at the right endpoint - derivative at the left endpoint).
    
    Arguments:
        x: Vector of x-values (independent variable).
        y: Vector of y-values (dependent variable).
        degree: Degree of the polynomial used for approximation.
    
    Returns:
        jump: A struct containing:
              - left_derivative: The derivative at the left endpoint.
              - right_derivative: The derivative at the right endpoint.
              - value: The jump (right_derivative - left_derivative).
    %}
    
    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Number of points for local fitting
    num_points = degree + 2; % Use degree + 2 points for better accuracy
    
    % Validate input
    if length(x) < num_points
        error('Not enough data points for the specified polynomial degree.');
    end
    
    % Left endpoint
    x_left = x(1:num_points);
    y_left = y(1:num_points);
    p_left = polyfit(x_left, y_left, degree); % Fit polynomial at left
    left_derivative = polyval(polyder(p_left), x(1)); % Compute derivative at x(1)
    
    % Right endpoint
    x_right = x(end-num_points+1:end);
    y_right = y(end-num_points+1:end);
    p_right = polyfit(x_right, y_right, degree); % Fit polynomial at right
    right_derivative = polyval(polyder(p_right), x(end)); % Compute derivative at x(end)
    
    % Compute jump
    jump.left_derivative = left_derivative;
    jump.right_derivative = right_derivative;
    jump.value = right_derivative - left_derivative;
end
