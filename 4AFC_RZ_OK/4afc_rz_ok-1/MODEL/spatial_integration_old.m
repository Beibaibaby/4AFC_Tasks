clear all
clc

% system parameters
space_l_center = 3; % covariance parameter space.
space_l_surr = 18; % covariance parameter space.
time_l_center = 40; % covariance parameter time.
time_l_surr = 60; % covariance parameter time.

prior_sigma = 0.2; %


% signal parameters
sensory_sigma = 1;
signal_max = 4;
stimulus_duration = 50;
radius = 2;
contrast = 0.1;

% house keeping parameters

space_range = 20;
time_range_max = 100; % ms.
time_range_min = -50; % ms.
space_steps = 30;
time_steps = 15;

dimention = 3;


%%%%%

if dimention == 2
    
    x = linspace( -0.5*space_range, 0.5*space_range, space_steps );
    y = linspace( -0.5*space_range, 0.5*space_range, space_steps );
    t = linspace( time_range_min, time_range_max, time_steps );
    
    [X, Y, T] = meshgrid( x, y, t );
    dist_from_zero = sqrt(X.^2 + Y.^2);
    
    stimulus = 0.5*(cos( pi*dist_from_zero/radius) + 1);
    stimulus( dist_from_zero > radius ) = 0;
    stimulus = signal_max*stimulus*contrast;
    stimulus( T < 0 | T > stimulus_duration ) = 0;
    
    XX = [X(:), Y(:), T(:)];
    measurements = stimulus(:);% + sensory_sigma * randn(length(stimulus(:)),1);
    
    [predicted_center] = gaussian_process( XX, measurements, XX,...
        prior_sigma, sensory_sigma, [space_l_center, space_l_center, time_l], 0 );
    
    predict_matrix_center = reshape( predicted_center, [space_steps, space_steps, time_steps] );
    
    %mesh( mean( predict_matrix( :,:, t >= 0 & t <= stimulus_duration ), 3 ));
    
    [predicted_surr] = gaussian_process( XX, measurements, XX,...
        prior_sigma, sensory_sigma, [space_l_surr, space_l_surr, time_l], 0 );
    
    predict_matrix_surr = reshape( predicted_surr, [space_steps, space_steps, time_steps] );
    
    mesh( mean( predict_matrix_center( :,:, t >= 0 & t <= stimulus_duration )...
        - predict_matrix_surr( :,:, t >= 0 & t <= stimulus_duration ), 3 ));
    
    
else
    
    x = linspace( -0.5*space_range, 0.5*space_range, space_steps );
    t = linspace( time_range_min, time_range_max, time_steps );
    
    [X, T] = meshgrid( x, t );
    dist_from_zero = abs(X);
    
    stimulus = 0.5*(cos( pi*dist_from_zero/radius) + 1);
    stimulus( dist_from_zero > radius ) = 0;
    stimulus = signal_max*stimulus*contrast;
    stimulus( T < 0 | T > stimulus_duration ) = 0;
    
    XX = [X(:), T(:)];
    measurements = stimulus(:);% + sensory_sigma * randn(length(stimulus(:)),1);
    
    [predicted_center] = gaussian_process( XX, measurements, XX,...
        prior_sigma, sensory_sigma, [space_l_center, time_l_center], 0 );
    
    predict_matrix_center = reshape( predicted_center, [time_steps, space_steps] );
    
    %mesh( mean( predict_matrix( :,:, t >= 0 & t <= stimulus_duration ), 3 ));
    
    [predicted_surr] = gaussian_process( XX, measurements, XX,...
        prior_sigma, sensory_sigma, [space_l_surr, time_l_surr], 0 );
    
    predict_matrix_surr = reshape( predicted_surr, [time_steps, space_steps] );
    
%     plot( x, mean( predict_matrix_center( t >= 0 & t <= stimulus_duration, : ))...
%         -mean( predict_matrix_surr( t >= 0 & t <= stimulus_duration, : ),1 ));
%     
    
    %figure
    
%     plot( t, predict_matrix_center( :, round( 0.5*space_steps ) )...
%         - predict_matrix_surr( :, round( 0.5*space_steps ) ) );
    
    plot( t, predict_matrix_center( :, round( 0.5*space_steps ) ));
    hold on
    plot( t, -predict_matrix_surr( :, round( 0.5*space_steps ) ));
    plot( t, predict_matrix_center( :, round( 0.5*space_steps ) )...
        - predict_matrix_surr( :, round( 0.5*space_steps ) ),'r' );
        


end








