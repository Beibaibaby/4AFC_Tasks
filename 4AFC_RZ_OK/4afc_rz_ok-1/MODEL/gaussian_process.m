function [predicted_y, variance_predict] = gaussian_process( X, y, X_star,...
    sigma_prior, sigma_noise, scale_factor, vari_too )

K = squared_exp_K( X, X, scale_factor, sigma_prior );
L = chol( K + sigma_noise^2*eye(size(X,1)),'lower' );
alpha = L'\(L\y);

K_star =  squared_exp_K( X, X_star, scale_factor, sigma_prior );

predicted_y = K_star'*alpha;


if vari_too == 1
    
    v = L\K_star;
    K_xx = squared_exp_K( X_star, X_star, scale_factor, sigma_prior );
    variance_predict = K_xx - v'*v;
    
end
