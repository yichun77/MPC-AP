% Add mcmcstat to path
addpath(genpath('mcmcstat'))

% Load your experimental data
data.time = time_vec;
data.CGM = CGM_data;
data.insulin = insulin_data;
data.meal_times = meal_times;
data.meal_amounts = meal_amounts;

% Parameter structure
params = {
    {'t_max_IA', 75, 10, 30, 120}
    {'S_I', 0.005, 0.001, 0.001, 0.02}
    {'X_b', 12, 2, 6, 20}
    {'t_max_G', 45, 5, 20, 70}
    {'A_G', 0.8, 0.1, 0.4, 1.2}
    {'K', 0.004, 0.0005, 0.001, 0.01}
    {'G_b', 6.5, 0.5, 4, 10}
};

% Setup options
options.nsimu = 5000;
model.ssfun = @ssfun;
model.sigma2 = 0.1^2;  % 10% CV
model.N = length(CGM_data);

% Run MCMC
results = mcmcrun(model, data, params, options);

% Plot results
mcmcplot(results);

best_params = median(results.chain);
G_sim = glucose_model(best_params, data);

figure;
plot(data.time, data.CGM, 'ko', 'DisplayName', 'Observed CGM');
hold on;
plot(data.time, G_sim, 'r-', 'LineWidth', 2, 'DisplayName', 'Simulated Glucose');
legend;
xlabel('Time (min)');
ylabel('Glucose (mmol/L)');
title('MCMC Fitting Result');
