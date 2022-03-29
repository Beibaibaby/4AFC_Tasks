function [mu_2AFC, mu_4AFC] = fit_and_plot_two_plus_four


load('LS_curves100.mat')
%data = data_all;
% load('OK_mask2_100_4.mat')
% data_all = [data; data_all];

data_all = data_all(101:end,:);

size( data_all )
options = optimset('MaxIter', 1e+5,'MaxFunEvals', 1e+5);

%LB = eps * ones(size(params0));%[ 0.001, 0.0001, -1, 0.00001];
%UB = 100 * ones(size(params0));%[ 200, inf, 100, 500];
%%%    curr_mu = params(1) .* size.^ (-2) + params(2) .* xs .^ params(3);

params0 = [log(0.015), 2 ];
%params0 = [0.5, -2, 0.1, 2, 2];


LB = [ 0.000001, 0.0001];
UB = [ 0.1, 0.1];

correct_2AFC = zeros( length(data_all), 1 );
correct_4AFC = zeros( length(data_all), 1 );

for i = 1:length( data_all )
    
    if data_all(i,4) == 1
        
        if data_all(i,6) == 79
            correct_4AFC(i) = 1;
            correct_2AFC(i) = 1;
        elseif data_all(i,6) == 80
            correct_2AFC(i) = 1;
        end
        
    elseif data_all(i,4) == 2
        
        if data_all(i,6) == 80
            correct_4AFC(i) = 1;
            correct_2AFC(i) = 1;
        elseif data_all(i,6) == 79
            correct_2AFC(i) = 1;
        end
        
    elseif data_all(i,4) == 3
        
        if data_all(i,6) == 81
            correct_4AFC(i) = 1;
            correct_2AFC(i) = 1;
        elseif data_all(i,6) == 82
            correct_2AFC(i) = 1;
        end
        
    elseif data_all(i,4) == 4
        
        if data_all(i,6) == 82
            correct_4AFC(i) = 1;
            correct_2AFC(i) = 1;
        elseif data_all(i,6) == 81
            correct_2AFC(i) = 1;
        end
    
    end
    
    
end

data_all = [data_all, correct_2AFC, correct_4AFC];

correct_trials  = data_all( correct_4AFC == 1, : );
incorrect_trials  = data_all( correct_4AFC == 0, : );

plot( correct_trials(:,2), correct_trials(:,3), 'b.' );
hold on
plot( incorrect_trials(:,2), incorrect_trials(:,3), 'r.' );



%for mask = 1:2
    size_all = unique( data_all(:,2) );
    
    for size_index = 1:length( size_all )
       
        % doi = data_all( data_all(:,7) == mask, : );
        doi = data_all( data_all(:,2) == size_all(size_index), : );
        
%         best_params_2AFC = fmincon(@(params) cal_likeli_MS( params, doi, 2 ), params0,[],[],[],[],LB,UB,[], options);
%         best_params_4AFC = fmincon(@(params) cal_likeli_MS( params, doi, 4 ), params0,[],[],[],[],LB,UB,[], options);
%         
        best_params_2AFC = fminsearch(@(params) cal_likeli_MS( params, doi, 2 ), params0, options);
        best_params_4AFC = fminsearch(@(params) cal_likeli_MS( params, doi, 4 ), params0, options);
        
        
        mu_2AFC(size_index) = best_params_2AFC(1);
        mu_4AFC(size_index) = best_params_4AFC(1);
        
        sigma_2AFC(size_index) = best_params_2AFC(2);
        sigma_4AFC(size_index) = best_params_4AFC(2);
        
        
    end
    
%end
mu_2AFC
mu_4AFC
sigma_2AFC
sigma_4AFC

% [trial, size_index, size, duration, stimulus_direction,
% correct(incomplete), response, mask_direction, 2AFC, 4AFC];

function log_likeli = cal_likeli_MS( params, data, n_alternatives )

durations = data(:,3);
durations = log(durations);
mu = (params(1));
sigma = params(2);

if n_alternatives == 2
    correct_or_not = data(:,8);
elseif n_alternatives == 4
    correct_or_not = data(:,9);
end
    

P_cdf = 0.5*erfc( -(durations-mu)./(sqrt(2)*sigma));% normcdf( durations, mu, sigma );
if n_alternatives == 2
    P_correct = 0.5 + 0.5*P_cdf;
elseif n_alternatives == 4
    P_correct = 0.25 + 0.75*P_cdf;
end
    
P_correct( P_correct == 1 ) = 1-eps;

P_correct( correct_or_not == 0 ) = 1 - P_correct( correct_or_not == 0 );

log_likeli = -sum( log( P_correct ));

