function K = squared_exp_K( X1, X2, scale_factor, sigma_prior )
% covariance matrix = exp( -0.5*( X1 - X2 )^2);
% 
% default: scale_factor = [1,1,1...]; sigma_prior = 1; sigma_noise = 0;

for i = 1:size( X1, 2 )
    
    temp_x1 = X1(:,i);
    temp_x2 = X2(:,i);
    
    x1 = repmat( temp_x1, 1, length(temp_x2) );
    x2 = repmat( temp_x2', length(temp_x1), 1 );
    
    if i == 1
        square_diff = ((x1-x2)/scale_factor(i)).^2;
    else
        square_diff = square_diff + ((x1-x2)/scale_factor(i)).^2;
    end
    
end

K = sigma_prior^2*exp( - 0.5 * square_diff );