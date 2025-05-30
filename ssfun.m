function ss = ssfun(params, data)
    G_sim = glucose_model(params, data);
    G_obs = data.CGM;
    ss = sum((G_obs - G_sim).^2); % Sum of squares
end
