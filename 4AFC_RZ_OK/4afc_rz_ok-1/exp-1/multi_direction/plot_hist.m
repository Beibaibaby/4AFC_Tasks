subject_initials = 'OK';
speed = 5;
for block = 1:10
    
    file_name = strcat(subject_initials,'_block',num2str(block),'_',num2str(speed),'.mat');
    load( file_name );
    
    if block == 1
        data = data_all;
    else
        data = [data; data_all];
    end
    
end



%data = data( data(:,3) > 45 & data(:,3) < 135, : );
%data = data_all;
durations_all = unique( data(:,2) );

%Xs = -170:20:170;
Xs = -175:10:175;
for i = 1:length(durations_all)
    
    duration = durations_all(i);
    
    doi = data( data(:,2) == duration, : );
    diff_all = doi(:,3) - doi(:,4);
    diff_all = diff_all - 90;
    diff_all( diff_all < - 180 ) = 360 + diff_all( diff_all < - 180 );
    diff_all( diff_all >= 180 ) = diff_all( diff_all >= 180 ) - 360;
    
    %subplot(length( durations_all ), 1, i),histnorm( diff_all, Xs );
    subplot(length( durations_all ), 1, i),hist( diff_all, Xs );
    xlim([-180, 180])
    %ylim([0, 0.03])
    
end
