function [mse, predicts, mu] = cal_mse_supp_curve( params, data, fixed_params )


% data.size, data.low, data.high
% fixed_params.signal_max

contrasts = data.contrast;
sizes = data.sizes;
low_th = data.low;
high_th = data.high;

%signal_max = fixed_params.signal_max;
space_range = fixed_params.space_range;
space_steps = fixed_params.space_steps;
sensory_sigma = fixed_params.sensory_sigma;

space_l_center = params(1);
space_l_surr = params(2);

prior_sigma_low = params(3);
prior_sigma_high = params(4);
signal_max = params(7);

lambda_center = params(5);
lambda_surr = params(6);
%signal_max = 5;
%k = params(8);
k = 0.5;%0.5
duration = 1:300;

%%%

x = linspace( -0.5*space_range, 0.5*space_range, space_steps );
x = x';
predicts = zeros( length(contrasts), length(sizes) );

for radi_index = 1:length(sizes)
    
    radius = sizes(radi_index);
    for cont_index = 1:length(contrasts)
        
        contrast = contrasts(cont_index);
        dist_from_zero = abs(x);
        
        stimulus = 0.5*(cos( pi*dist_from_zero/radius) + 1);
        stimulus( dist_from_zero > radius ) = 0;
        stimulus = signal_max*stimulus*contrast/100;
        measurements = stimulus(:);
        
        if cont_index == 1
            prior_sigma = prior_sigma_low;
        else
            prior_sigma = prior_sigma_high;
        end
        
        [predicted_center, s2_center] = gaussian_process( x, measurements, x,...
            prior_sigma, sensory_sigma, [space_l_center], 1 );
        
        [predicted_surr, s2_surr] = gaussian_process( x, measurements, x,...
            prior_sigma, sensory_sigma, [space_l_surr], 1 );
        
        [mu, sigma] = cal_response_dist( predicted_surr, s2_surr, predicted_center, s2_center, k, radius, x,...
            lambda_center, lambda_surr, duration, stimulus );
        
        temp_thres = find( mu - sigma > 0 );
        if isempty(temp_thres)
            predicts(cont_index, radi_index) = max(duration);
        else
            predicts(cont_index, radi_index) = temp_thres(1);
        end
        
    end
end

data_all = [ low_th; high_th ];
diff_all = data_all(:) - predicts(:);

mse = diff_all'*diff_all;

if prior_sigma_low < 0 || prior_sigma_high <0 || params(6) < 20
 mse = mse + 100000;
end

mse
