clear
clc
% 
% data.contrast = [2, 50];
% data.sizes = [1, 1.5, 3, 6, 11];
% data.low = [70, 50, 30, 15, 10];
% data.high = [ 10, 15, 30, 50, 70]; 
% 
data.contrast = [3, 30];
data.sizes = [0.7, 1.3, 2.7, 4, 5];
data.low = [92, 50, 42, 39, 40];
data.high = [ 32, 37, 80, 80, 130]; 

fixed_params.space_range = 100;
fixed_params.space_steps = 200;
fixed_params.sensory_sigma = 1;

% params0(1) = 2;% space_l_center
% params0(2) = 30;% space_l_surr
% params0(3) = 1;% prior_sigma_low
% params0(4) = 20;% prior_sigma_high
% params0(5) = 5;% signal_max
% params0(6) = 70;%lambda_center
% params0(7) = 30;%lambda_surr

params0(1) = 2;% space_l_center
params0(2) = 50;% space_l_surr
params0(3) = 1;% prior_sigma_low
params0(4) = 20;% prior_sigma_high
params0(7) = 5;% signal_max
params0(5) = 70;%lambda_center
params0(6) = 30;%lambda_surr
% params0(8) = 0.5;

options = optimset( 'MaxIter', 100000 );
options = optimset( 'MaxFunEvals', 100000 );

best_mse = 10000000;

for i = 1:8
%while 1    
    
    params_best = fminsearch( @(params) cal_mse_supp_curve( params, data, fixed_params ), params0, options );
    [mse, predicts, mus] = cal_mse_supp_curve( params_best, data, fixed_params );
    %params_best
    if mse < best_mse
        best_mse = mse;
        params_best_sofar = params_best;
        best_mus = mus;
        
        params0 = params_best + randn( size(params0))*diag(0.1*params_best_sofar);
        while params0(6) < 20
            params0 = params_best + randn( size(params0))*diag(0.1*params_best_sofar);
        end
        
    else
        
        
    end
end

[mse, predicts, model_all] = cal_mse_supp_curve_lookat( params_best_sofar, data, fixed_params );
hold off
plot( predicts' );
hold on
plot( data.low );
plot( data.high );






