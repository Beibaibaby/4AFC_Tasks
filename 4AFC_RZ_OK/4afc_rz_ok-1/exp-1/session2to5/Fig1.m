clear all
clc

%subject_initials = 'MM_low';
subject_initials = 'FP_high';

speed = 5;
for block = 1:10
    
    file_name = strcat(subject_initials,'_block',num2str(block),'_',num2str(speed),'.mat');
    load( file_name, 'data_all' );
    
    if block == 1
        data = data_all;
    else
        data = [data; data_all];
    end
    
end
data_high = data;
clear data

subject_initials = 'FP_low';

speed = 5;
for block = 1:10
    
    file_name = strcat(subject_initials,'_block',num2str(block),'_',num2str(speed),'.mat');
    load( file_name, 'data_all' );
    
    if block == 1
        data = data_all;
    else
        data = [data; data_all];
    end
    
end

data_low = data;


%data = data( data(:,3) > 45 & data(:,3) < 135, : );
%data = data_all;
durations_all = unique( double(data(:,2)) );
durations_all = durations_all( find(durations_all ~= 0 ));


%Xs = -170:20:170;
Xs = -175:10:175;


figure

for i = 1:length(durations_all)
    
    for cont = 1:2
        
        duration = durations_all(i);
        
        if cont == 1
            data = data_low;
        else
            data = data_high;
        end
        
        doi = data( data(:,2) == duration, : );
        diff_all = doi(:,3) - doi(:,4);
        diff_all = diff_all - 90;
        diff_all( diff_all < - 180 ) = 360 + diff_all( diff_all < - 180 );
        diff_all( diff_all >= 180 ) = diff_all( diff_all >= 180 ) - 360;
        
        %subplot(length( durations_all ), 1, i),histnorm( diff_all, Xs );
        subplot(length(durations_all), 2, i*2-2 + cont),hist( diff_all, Xs );
        xlim([-180, 180])
        %ylim([0, 0.03])
        
    end
    
end


