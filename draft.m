clear; clc; close all;

%% 1. 定义葡萄糖-胰岛素调节模型
function dydt = glucose_insulin_model(t, y, params, u_I, u_G_tj, tj)
    % 参数解包
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
    
    % 状态变量
    x1 = y(1);  % 胰岛素吸收第一室
    x2 = y(2);  % 胰岛素吸收第二室
    a1 = y(3);  % 膳食吸收第一室
    a2 = y(4);  % 膳食吸收第二室
    G = y(5);   % 葡萄糖浓度
    
    % 胰岛素吸收子模型
    dx1dt = -x1/t_max_IA + u_I/60;
    dx2dt = (x1 - x2)/t_max_IA;
    X = 1000 * x2 / (t_max_IA * MCR_I * W);
    
    % 膳食吸收子模型
    delta = (t == tj);  % Dirac delta函数近似
    da1dt = -a1/t_max_G + delta * u_G_tj;
    da2dt = (a1 - a2)/t_max_G;
    U_M = 5.556 * A_G * a2 / (t_max_G * V_G * W);
    
    % 葡萄糖动力学
    dGdt = -S_I*(X - X_b) + U_M - K*(G - G_b);
    
    dydt = [dx1dt; dx2dt; da1dt; da2dt; dGdt];
end

function result = f3(G, time_scale)
    C3 = 1000;
    Vg = 10;
    result = G / (C3 * Vg);
end


%% 2. 简化的MCMC参数估计示例
% 假设已有实验数据：time_vec, CGM_data, insulin_data, meal_times, meal_amounts
% 定义先验分布和参数范围
params_init = struct('t_max_IA', 75, 'S_I', 0.005, 'X_b', 12, ...
                    't_max_G', 45, 'A_G', 0.8, 'K', 0.004, 'G_b', 6.5);
param_names = fieldnames(params_init);
n_iter = 1000;  % 简化迭代次数

% MCMC链初始化
chain = struct();
for p = 1:length(param_names)
    chain.(param_names{p}) = zeros(n_iter,1);
end
current_params = params_init;
current_logpost = calculate_log_posterior(current_params, @glucose_insulin_model, time_vec, CGM_data);

for i = 1:n_iter
    % 参数提议
    proposed_params = current_params;
    for p = 1:length(param_names)
        proposed_params.(param_names{p}) = current_params.(param_names{p}) + randn*0.1;
    end
    
    % 计算后验概率
    proposed_logpost = calculate_log_posterior(proposed_params, @glucose_insulin_model, time_vec, CGM_data);
    
    % Metropolis-Hastings接受准则
    if exp(proposed_logpost - current_logpost) > rand
        current_params = proposed_params;
        current_logpost = proposed_logpost;
    end
    
    % 保存链
    for p = 1:length(param_names)
        chain.(param_names{p})(i) = current_params.(param_names{p});
    end
end

%% 3. 模型仿真
function G_sim = simulate_model(params, time_vec, insulin_profile, meal_times, meal_amounts)
    y0 = [0; 0; 0; 0; params.G_b];  % 初始条件
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    
    % 创建胰岛素输入函数
    u_I = @(t) interp1(time_vec, insulin_profile, t);
    
    % 创建膳食输入函数
    u_G = @(t) sum(meal_amounts.*(t == meal_times));  % 简化处理
    
    % 解微分方程
    [t_out, y_out] = ode45(@(t,y) glucose_insulin_model(t, y, params, u_I(t), u_G(t), t),...
                          time_vec, y0, options);
    
    G_sim = y_out(:,5);
end

%% 辅助函数：计算对数后验概率
function logpost = calculate_log_posterior(params, model_func, time_vec, CGM_data)
    % 先验分布（示例：正态分布）
    prior_means = struct('t_max_IA', 75, 'S_I', 0.005, 'X_b', 12, ...
                       't_max_G', 45, 'A_G', 0.8, 'K', 0.004, 'G_b', 6.5);
    prior_sds = struct('t_max_IA', 10, 'S_I', 0.001, 'X_b', 2, ...
                      't_max_G', 5, 'A_G', 0.1, 'K', 0.0005, 'G_b', 0.5);
    
    % 计算先验概率
    log_prior = 0;
    param_names = fieldnames(params);
    for p = 1:length(param_names)
        log_prior = log_prior + log(normpdf(params.(param_names{p}),...
                                          prior_means.(param_names{p}),...
                                          prior_sds.(param_names{p})));
    end
    
    % 计算似然（假设测量误差为高斯分布）
    G_sim = simulate_model(params, time_vec, insulin_data, meal_times, meal_amounts);
    residuals = CGM_data - G_sim;
    log_likelihood = -0.5 * sum((residuals/0.1).^2);  # 假设CV=10%
    
    logpost = log_prior + log_likelihood;
end