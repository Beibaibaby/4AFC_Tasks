function minus_log_likeli = cal_likeli_SS( params, data, n_alternatives )

% params 1-3 are for SS function, params 4 is sigma. 
% data : stimulus type, size, duration, response, correct.
% size, duration, stimulus, correct_not, response, 2 or 4    

%types = unique(data(:,1));
%types      =[1 2 3 4 5 6];
%sizes     = exp(linspace( log(0.7), log(8), 6 ));
sizes = unique( data(:,1) );
for i = 1:length(sizes)
    
    doi = data( data(:,1) == sizes(i), :);
    size = sizes(i);
    if length(params) == 4
        curr_mu = params(1) .* size.^ (-2) + params(2) .* size .^ params(3);
        log_likeli_duration(i) = cal_likeli_duration( doi, curr_mu, params(4), n_alternatives);
    elseif length(params) == 5
        curr_mu = params(1) .* size.^ params(2) + params(3) .* size .^ params(4);
        log_likeli_duration(i) = cal_likeli_duration( doi, curr_mu, params(5), n_alterantives);
    end
    
    
end

minus_log_likeli = -sum(log_likeli_duration);

function log_likeli = cal_likeli_duration( data, mu, sigma, n_alterantives )

% size, duration, stimulus, correct_not, response, 2 or 4    

durations = data(:,2);
durations = log(durations);
mu = log(mu);
correct_or_not = data(:,4);

P_cdf = 0.5*erfc( -(durations-mu)./(sqrt(2)*sigma));% normcdf( durations, mu, sigma );
if n_alterantives == 2
    P_correct = 0.5 + 0.5*P_cdf;
elseif n_alterantives == 4
    P_correct = 0.25 + 0.75*P_cdf;
end
    
P_correct( P_correct == 1 ) = 1-eps;

P_correct( correct_or_not == 0 ) = 1 - P_correct( correct_or_not == 0 );

log_likeli = sum( log( P_correct ));