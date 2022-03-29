clear all
clc

% system parameters
space_l_center = 3; % covariance parameter space.
space_l_surr = 15; % covariance parameter space.
%time_l_center = 40; % covariance parameter time.
%time_l_surr = 60; % covariance parameter time.

prior_sigma = 0.5; %

% signal parameters
sensory_sigma = 1;
signal_max = 5;
%stimulus_duration = 50;
radius = 10;
contrast = 1;

% house keeping parameters

space_range = 100;
%time_range_max = 100; % ms
%time_range_min = -50; % ms.
space_steps = 200;
%time_steps = 15;

%dimention = 3;


%%%%%
x = linspace( -0.5*space_range, 0.5*space_range, space_steps );
x = x';
%t = linspace( time_range_min, time_range_max, time_steps );

%[X, T] = meshgrid( x, t );
dist_from_zero = abs(x);

stimulus = 0.5*(cos( pi*dist_from_zero/radius) + 1);
stimulus( dist_from_zero > radius ) = 0;
stimulus = signal_max*stimulus*contrast;
%stimulus( T < 0 | T > stimulus_duration ) = 0;

%XX = [X(:), T(:)];


for i = 1:10


measurements = stimulus(:);% + sensory_sigma * randn(length(stimulus(:)),1);

[predicted_center, s2_center] = gaussian_process( x, measurements, x,...
    prior_sigma, sensory_sigma, [space_l_center], 1 );

% predict_matrix_center = reshape( predicted_center, [time_steps, space_steps] );

%mesh( mean( predict_matrix( :,:, t >= 0 & t <= stimulus_duration ), 3 ));

[predicted_surr, s2_surr] = gaussian_process( x, measurements, x,...
    prior_sigma, sensory_sigma, [space_l_surr], 1 );

% predict_matrix_surr = reshape( predicted_surr, [time_steps, space_steps] );

%     plot( x, mean( predict_matrix_center( t >= 0 & t <= stimulus_duration, : ))...
%         -mean( predict_matrix_surr( t >= 0 & t <= stimulus_duration, : ),1 ));
%

%figure

%     plot( t, predict_matrix_center( :, round( 0.5*space_steps ) )...
%         - predict_matrix_surr( :, round( 0.5*space_steps ) ) );

%plot( x, predicted_center);
hold on
%plot( x, -predicted_surr);
plot( x, predicted_center - predicted_surr,'r' );
%plot( x, predicted_center - predicted_surr + sqrt(diag(s2_center) + diag(s2_surr)),'r' );
%plot( x, predicted_center - predicted_surr - sqrt(diag(s2_center) + diag(s2_surr)),'r' );

mean_center(i) = mean( predicted_center( x < radius ) );

% plot( x, stimulus, 'k');
end


k = 0.5;
lambda_center = 60;
lambda_surr = 40;
duration = 10;

[mu, sigma] = cal_response_dist( predicted_surr, s2_surr, predicted_center, s2_center, k, radius, x,...
    lambda_center, lambda_surr, duration );




