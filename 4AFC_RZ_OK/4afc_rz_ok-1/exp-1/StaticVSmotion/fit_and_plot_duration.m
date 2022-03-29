%clear all;close all;
clc;
clear

% load('OK_motionRect_15_block1_5.mat')
% %load('OK_static_10_average_duration1_100.mat')
% data = data_all;
% load('OK_motionRect_15_block2_5.mat')
% 
% data = [ data; data_all ];


load('OK_constant_block1_3.mat')
data = data_all;
% load('SA_constant_block2_100.mat')
% data = [ data; data_all ];
%load('OK_static_10_average_duration1_100.mat')
%data = data_all;
% load('KS_station3Dot_15_block2_5.mat')
% % %
% data = [ data; data_all ];




data_all = data;
clear size

size( data_all )
%data_all = data_all( data_all(:,6) < 81, :);

correct_or_not = zeros(size(data_all,1),1);
for i = 1:size(data_all, 1)
    angle = data_all(i,3);
    response = data_all(i,4);
    
    if angle == 1
        if response == 79
            correct_or_not(i) = 1;
        elseif response == 80
            correct_or_not(i) = -1;
        else
            correct_or_not(i) = 0;
        end
        
    elseif angle == 2
        if response == 80
            correct_or_not(i) = 1;
        elseif response == 79
            correct_or_not(i) = -1;
        else
            correct_or_not(i) = 0;
        end
        
    elseif angle == 3
        if response == 81
            correct_or_not(i) = -1;
        elseif response == 82
            correct_or_not(i) = 1;
        else
            correct_or_not(i) = 0;
        end
    elseif angle == 4
        if response == 82
            correct_or_not(i) = -1;
        elseif response == 81
            correct_or_not(i) = 1;
        else
            correct_or_not(i) = 0;
        end
    end
       
end

data_all = [data_all, correct_or_not];
durations = unique( data_all(:,2) );
correct_ratio = zeros( 1, length(durations));
opposite_ratio = zeros( 1, length(durations));
ortho_ratio = zeros( 1, length(durations));


%first, parallel masking
for i = 1:length(durations)
    doi = data_all( data_all(:,2) == durations(i), : );
    %doi=doi(doi(:,8)==1,:);
    correct_ratio(i) = sum( doi(:,end) == 1 )/size(doi, 1);
    opposite_ratio(i) = sum( doi(:,end) == -1 )/size(doi, 1);
    ortho_ratio(i) = sum( doi(:,end) == 0 )/size(doi, 1);
    size(doi,1);
    
end
%plot parallel mask
%figure
hold on
plot( durations, correct_ratio, 'bo-' );
    hold on;
plot( durations, opposite_ratio, 'ro-' );
plot( durations, ortho_ratio, 'ko-' );


