clear; clc; close all;

%% Glucose-Insulin Dynamics Model
function dydt = glucose_insulin_model(t, y, params, u_I, u_G_tj, tj)
    % Parameters
    t_max_IA = params.t_max_IA;
    MCR_I = 0.017;
    W = params.W;
    t_max_G = params.t_max_G;
    A_G = params.A_G;
    V_G = 0.16;
    S_I = params.S_I;
    X_b = params.X_b;
    K = params.K;
    G_b = params.G_b;

    % State variables
    x1 = y(1); % Insulin absorption compartment 1
    x2 = y(2); % Insulin absorption compartment 2
    a1 = y(3); % Meal absorption compartment 1
    a2 = y(4); % Meal absorption compartment 2
    G  = y(5); % Glucose concentration

    % Insulin absorption sub-model
    dx1dt = -x1 / t_max_IA + u_I / 60;
    dx2dt = (x1 - x2) / t_max_IA;
    X = 1000 * x2 / (t_max_IA * MCR_I * W);

    % Meal absorption sub-model
    delta = (t == tj); % Approximation of Dirac delta function
    da1dt = -a1 / t_max_G + delta * u_G_tj;
    da2dt = (a1 - a2) / t_max_G;
    U_M = 5.556 * A_G * a2 / (t_max_G * V_G * W);

    % Glucose dynamics
    dGdt = -S_I * (X - X_b) + U_M - K * (G - G_b);

    dydt = [dx1dt; dx2dt; da1dt; da2dt; dGdt];
end

%% Utility function
function result = f3(G, time_scale)
    C3 = 1000;
    Vg = 10;
    result = G / (C3 * Vg);
end

%% Simplified MCMC Parameter Estimation
% Assumes experimental data: time_vec, CGM_data, insulin_data, meal_times, meal_amounts

% Initial parameter settings
params_init = struct('t_max_IA', 75, 'S_I', 0.005, 'X_b', 12, ...
                     't_max_G', 45, 'A_G', 0.8, 'K', 0.004, 'G_b', 6.5);
param_names = fieldnames(params_init);
n_iter = 1000;

% Initialize MCMC chains
chain = struct();
for p = 1:length(param_names)
    chain.(param_names{p}) = zeros(n_iter, 1);
end

% Evaluate initial posterior
current_params = params_init;
current_logpost = calculate_log_posterior(current_params, @glucose_insulin_model, time_vec, CGM_data);

% MCMC iterations
for i = 1:n_iter
    proposed_params = current_params;
    for p = 1:length(param_names)
        proposed_params.(param_names{p}) = current_params.(param_names{p}) + randn * 0.1;
    end
    
    proposed_logpost = calculate_log_posterior(proposed_params, @glucose_insulin_model, time_vec, CGM_data);

    % Metropolis-Hastings acceptance rule
    if exp(proposed_logpost - current_logpost) > rand
        current_params = proposed_params;
        current_logpost = proposed_logpost;
    end
    
    % Save to chain
    for p = 1:length(param_names)
        chain.(param_names{p})(i) = current_params.(param_names{p});
    end
end

%% Model Simulation
function G_sim = simulate_model(params, time_vec, insulin_profile, meal_times, meal_amounts)
    y0 = [0; 0; 0; 0; params.G_b]; % Initial conditions
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

    % Define input functions
    u_I = @(t) interp1(time_vec, insulin_profile, t, 'linear', 'extrap');
    u_G = @(t) sum(meal_amounts .* (t == meal_times)); % Simplified meal input

    % Solve ODE
    [~, y_out] = ode45(@(t, y) glucose_insulin_model(t, y, params, u_I(t), u_G(t), t), ...
                       time_vec, y0, options);

    G_sim = y_out(:, 5); % Extract glucose trajectory
end

%% Posterior Probability Calculation
function logpost = calculate_log_posterior(params, model_func, time_vec, CGM_data)
    % Prior distributions (assumed normal)
    prior_means = struct('t_max_IA', 75, 'S_I', 0.005, 'X_b', 12, ...
                         't_max_G', 45, 'A_G', 0.8, 'K', 0.004, 'G_b', 6.5);
    prior_sds = struct('t_max_IA', 10, 'S_I', 0.001, 'X_b', 2, ...
                       't_max_G', 5, 'A_G', 0.1, 'K', 0.0005, 'G_b', 0.5);

    % Calculate log-prior
    log_prior = 0;
    param_names = fieldnames(params);
    for p = 1:length(param_names)
        log_prior = log_prior + log(normpdf(params.(param_names{p}), ...
                                            prior_means.(param_names{p}), ...
                                            prior_sds.(param_names{p})));
    end

    % Simulate glucose response and compute log-likelihood
    G_sim = simulate_model(params, time_vec, insulin_data, meal_times, meal_amounts);
    residuals = CGM_data - G_sim;
    log_likelihood = -0.5 * sum((residuals / 0.1).^2); % Assuming CV = 10%

    % Return log-posterior
    logpost = log_prior + log_likelihood;
end
