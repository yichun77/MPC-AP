function y = glucose_model(params, data)
    % Extract parameters
    p = struct('t_max_IA', params(1), ...
               'S_I',     params(2), ...
               'X_b',     params(3), ...
               't_max_G', params(4), ...
               'A_G',     params(5), ...
               'K',       params(6), ...
               'G_b',     params(7), ...
               'W',       70);  % Assume 70kg for example

    time = data.time;
    insulin = data.insulin;
    meal_times = data.meal_times;
    meal_amounts = data.meal_amounts;

    % Define input functions
    u_I = @(t) interp1(time, insulin, t, 'linear', 'extrap');
    u_G = @(t) sum(meal_amounts .* (t == meal_times));

    % Initial condition
    y0 = [0; 0; 0; 0; p.G_b];
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

    % Simulate the model
    [~, y_out] = ode45(@(t, y) glucose_insulin_model(t, y, p, u_I(t), u_G(t), t), ...
                       time, y0, options);

    y = y_out(:, 5); % Return simulated glucose
end
