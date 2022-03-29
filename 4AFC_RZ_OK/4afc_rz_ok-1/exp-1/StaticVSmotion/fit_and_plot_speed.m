function thresholds = fit_and_plot_speed


options = optimset('MaxIter', 1e+5,'MaxFunEvals', 1e+5);
params0 = [log(0.015), 1 ];

subj = 'LS_B';

index = 1;
for block = 1:3
    filename = strcat(subj,'_block',num2str(block),'.mat');
    load(filename);
    
    if index == 1
        data = data_all;
        index = index + 1;
    else
        data = [ data; data_all ];
    end
end

condition_list = unique( data(:, 6 ) );

for condi_index = 1:length( condition_list )
    
    condi = condition_list(condi_index);
    data_condi = data( data(:,6) == condi, :);
    mean( data_condi ) 
    duration_list = unique( data_condi(:,2) );
    hit_rate = zeros( length( duration_list ), 1 );
    for dura_index = 1:length( duration_list )
        
        dura = duration_list( dura_index );
        data_dura = data_condi( data_condi(:,2) == dura, : );
        hit_rate(dura_index) = mean( data_dura(:,5) );
    end
    switch condi_index
        case 1
    plot( (duration_list), hit_rate, 'r' )
    hold on
        case 2
            plot( (duration_list), hit_rate, 'g' )
        case 3
            plot( (duration_list), hit_rate, 'b' )
    end
    
    
    best_params = fminsearch(@(params) cal_likeli_MS( params, data_condi, 2 ), params0, options);
    
    thresholds(condi_index) = best_params(1);
    
    
    
end



function log_likeli = cal_likeli_MS( params, data, n_alternatives )

durations = data(:,2);
durations = log(durations);
mu = (params(1));
sigma = params(2);


correct_or_not = data(:,5);
    

P_cdf = 0.5*erfc( -(durations-mu)./(sqrt(2)*sigma));% normcdf( durations, mu, sigma );
if n_alternatives == 2
    P_correct = 0.5 + 0.5*P_cdf;
elseif n_alternatives == 4
    P_correct = 0.25 + 0.75*P_cdf;
end
    
P_correct( P_correct == 1 ) = 1-eps;

P_correct( correct_or_not == 0 ) = 1 - P_correct( correct_or_not == 0 );

log_likeli = -sum( log( P_correct ));



