%clear all;close all;clc


load('FP_curves100.mat');
n_alternatives = 2;

color='kr';
%color='bg';



%load('OK_20_300.mat');
hold on

initial_trial = 1;
trialFit =200;

data_all = data_all(initial_trial:trialFit, :);
data_all(:,3) = data_all(:,3);
% 
% 
% if n_alternatives == 2
% data_all = data_all( data_all(:,7) == 2, : );
% elseif n_alternatives == 4
% data_all = data_all( data_all(:,7) == 1, : );
% end




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


% for i=1:trialFit
% [time_gauss,data_all(i,3)] = envelope( 0.5*data_all(i,3), frame_rate, which_envelope, amplitude );
% data_all(i,3)=data_all(i,3)/360;
% end

correct_data = data_all( data_all(:,5) == 1, :);
incorrect_data = data_all( data_all(:,5) == 0, :);


% responses = zeros( length(data_all), 1);
% responses(data_all(:,7) > 80) = 270;
% 
% correct_data = data_all( data_all(:,6) == responses, : );
% incorrect_data = data_all( data_all(:,6) ~= responses, : );




plotTrials =1;
% if plotTrials==1 
% loglog( correct_data(initial_trial:end,2), correct_data(initial_trial:end,3),[color(1) '.'],'MarkerSize',10);
% hold on
% loglog( incorrect_data(initial_trial:end,2), incorrect_data(initial_trial:end,3),[color(2) '.'],'MarkerSize',10)
% hold on
% end

options = optimset('MaxIter', 1e+5,'MaxFunEvals', 1e+5);

%LB = eps * ones(size(params0));%[ 0.001, 0.0001, -1, 0.00001];
%UB = 100 * ones(size(params0));%[ 200, inf, 100, 500];
%%%    curr_mu = params(1) .* size.^ (-2) + params(2) .* xs .^ params(3);

params0 = [0.2, 0.01, 2, 1];
%params0 = [0.01, -2, 0.5, 0.5, 2];

if length(params0) == 4
    LB = [ eps, eps, eps, eps];
    UB = [ 15, 15, 15, 12];
elseif length(params0) == 5
    LB = [ eps, -5, eps, eps, eps];
    UB = [ 5, 0, 5, 5, 2];
end


best_params = fmincon(@(params) cal_likeli_SS( params, data_all, n_alternatives ), params0,[],[],[],[],LB,UB,[], options)
%best_params = fminunc(@(params) cal_likeli_SS( params, data ), params0, options);

sizes_point = unique( data_all(:,2));
sizes = linspace( min(sizes), max(sizes), 100 );
if length(params0) == 4
    predict_y = best_params(1) .* sizes.^ (-2) + best_params(2) .* sizes .^ best_params(3);
elseif length(params0) == 5
    predict_y = best_params(1) .* sizes.^ best_params(2) + best_params(3) .* sizes .^ best_params(4);
end

if length(params0) == 4
    predict_y_point = best_params(1) .* sizes_point.^ (-2) + best_params(2) .* sizes_point .^ best_params(3);
elseif length(params0) == 5
    predict_y_point = best_params(1) .* sizes_point.^ best_params(2) + best_params(3) .* sizes_point .^ best_params(4);
end
    


hold on
subplot(1,2,1),plot( sizes, predict_y,'r' )
hold on
subplot(1,2,1),plot( sizes_point, predict_y_point,'ro' )
%xlim([0.2, 10.3]);
%ylim([0, 0.120])
predict_y
hold off