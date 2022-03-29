function [mu, sigma] = cal_response_dist( mean_surr, s2_surr, mean_center, s2_center, k, radius, x,...
    lambda_center, lambda_surr, duration, stimulus )

stimulus = stimulus/sum(stimulus);
mu_center = mean( mean_center( abs(x) < radius ) );
mu_surr = mean( mean_surr( abs(x) < radius ) );
% % 
% mu_center = mean_center'*stimulus;
% mu_surr = mean_surr'*stimulus;


mu_center = mu_center * ( 1 - exp( -duration/lambda_center ) );
mu_surr = mu_surr * ( 1 - exp( -duration/lambda_surr ) );

mu = mu_center - mu_surr;

s2_center = diag( s2_center );
s2_surr = diag( s2_surr );

s2_max = var( mean_center( abs(x) < radius ) ) + var( mean_surr( abs(x) < radius ) ) + ...
    mean( s2_center( abs(x) < radius ) ) + mean( s2_surr( abs(x) < radius ) );

sigma = k * sqrt( s2_max );
sigma = sigma*(1./sqrt(duration));