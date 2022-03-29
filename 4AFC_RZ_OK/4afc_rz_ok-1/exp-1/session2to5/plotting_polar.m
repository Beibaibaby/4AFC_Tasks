clear
clc

for block_number = 1:10
    
    data_name = strcat( 'FP_high_block', num2str(block_number),'_5.mat');
    load( data_name );
    
    if block_number == 1
        data = data_all;
    else
        data = [data; data_all];
    end
    
end


duration_list = unique( data(:,2));
doi = data( data(:,2) == duration_list(6), : );
%doi = doi(1:450, :);

diffs = doi(:,3) - doi(:,4);
diffs( diffs < -180 ) = 360 + diffs( diffs < -180 );
diffs( diffs >= 180 ) = diffs( diffs >= 180 ) - 360;


diffs = diffs*pi/180;

[n, X] = hist( diffs, 30 );

%[n, X] = ksdensity( diffs, -0.5:0.1:(2*pi + 0.5), 'width', 0.15 );

theta = [X, X(1)];
rho = [n, n(1)];

r_max  = max(n);
h_fake = polar(theta, r_max*ones(size(theta)));
hold on;

set(h_fake, 'Visible', 'Off');

polar( theta, rho );
%compass( theta, rho );

%rose( theta, rho );

