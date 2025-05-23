clear; clc; close all;

%% Parameters and Constantes
data = load("cam01.mat");  % 01-08

W = data.data.Weight_KG;  % Body weight (KG)
MCR_I = 0.017;  % Metabolic clearance rate of effective insulin (l/kg/min)
V_G = 0.16;  % plasma glucose pool size (l/kg)

%% Submodels (3)
% Insulin absorption and action
syms x1(t) x2(t) X(t) t_maxIA uI(t)
iaa1 = diff(x1, t) == -(1/t_maxIA)*x1+(uI/60);
x1Sol(t) = dsolve(iaa1);
iaa2 = diff(x2, t) == (1/t_maxIA)*(x1-x2);
x2Sol(t) = dsolve(iaa2);
%X = @(t) (1000*x2Sol) / (t_maxIA*MCR_I*W);
X(t) = (1000*x2Sol) / (t_maxIA*MCR_I*W);
X = simplify(X(t));

% Meal absorption dynamics
syms a1(t) a2(t) t_maxG deltatj(t) u_G(t_j) U_M(t) A_G
mad1 = diff(a1, t) == -(1/t_maxG)*a1+deltatj*(u_G(t_j));
a1Sol(t) = dsolve(mad1);
mad2 = diff(a2, t) == (1/t_maxG)*(a1-a2);
a2Sol(t) = dsolve(mad2);
%U_M = @(t) (5.556*A_G*a2Sol)/(t_maxG*V_G*W);
U_M(t) = (5.556*A_G*a2Sol)/(t_maxG*V_G*W);
U_M = simplify(U_M(t));

% Glucose dynamics
syms G(t) S_I X_b K G_b
gd = diff(G, t) == -S_I*(X-X_b)+U_M-K*(G-G_b);
gdSol(t) = dsolve(gd);

%% Parameter Estimation
% t_maxIA, t_maxG, A_G, S_I, X_b, K, G_b


