
%%%%%%%for 4-AFC motion program, using dots to get 4 AFC, stimulus can be motion or stationary %%%%%%%%%%%%

clc;clear all;close all;

warning('off','MATLAB:dispatcher:InexactMatch');
ListenChar(2);
% KEY Monitor Parameters
scale_factor   = 2.07;     % most important parameter - how many acrmin is one screen pixel? for SONY monitor at [1024 640] resolution, and 30.4in vieving distance, use scale_factor=2
frame_rate     = 360;        %
resolution     = [1280 720]; %
linearize 	   = 1;         % whether monitor is linearized. if=1, program will look for "MyGammaTable"
break_if_mismatch = 1;       % if 1, this will quit the program if monitor resolution does not match above settings (use 1 for the ACTUAL EXP, 0 for testing on a different monitor, eg laptop)


%% you need to change parameters here every time you run the code

subject_initials    ='LS_B';
block               = 3;  % 0 means test 
n_repeat = 90;
stair_condition = [1,2];
up_down_design = [1, 2; 1, 3];
stimulus_type = [3]; % 1 motion, 2: average, 3: enlongated
use_old_stair = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env_radi            = 11;    %deg, radius
contrast            = 100;      %
mask                = 0;  % 0 no mask, 1 dots, 2 texture
mask_duration       = 0.3;
speed               = 5;  %deg/sec


%%

% duration_list       = [logspace(log10(dura_range(1)), log10(dura_range(2)), 12)];
savefile            =1;
tme                 =clock;
filename = strcat(subject_initials,'_block',num2str(block),'.mat');
stair_name_save = strcat(subject_initials,'_curr_stair','.mat');
stair_name_load = strcat(subject_initials,'_curr_stair','.mat');

IsExist = exist(filename,'file');
if IsExist && block% ~= 0
    ListenChar(1);
    error('data file name exists')
end
if block==0
    savefile=0;
end

%%
%--------STIMULI parameters-------------------------------

SF                     = 1.2;          %cycles/deg

orientation            = 0;   %deg
spatial_envelope       = 2;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
which_envelope 	       = 2;
background             = 126;   % background intensity, in gray scale untis
white                  = background*2;
% n_total_dots           = 400;
% dot_size               = 60; % max 63.37
n_total_dots           = 10000;
dot_size               = 15; % max 63.37

mask_dot_size = 15;
mask_n_total_dots = 10000;

n_alternatives         = 4; % 2 or 4.
pixStep                = speed*60 / (frame_rate * scale_factor);        %moving speed, pixel/frame


%%
%-------other parameters;
l = 7;
ww = 4;
%---------------------------------------
%%


%% upper_limit = round(upper_limit/(1000/frame_rate));
H_ecc_fix       = 0;      H_ecc_fix = H_ecc_fix*60/scale_factor;
V_ecc_fix       = 0;      V_ecc_fix = V_ecc_fix*60/scale_factor;
H_ecc_stim      = 0;      H_ecc_stim = round(H_ecc_stim*60/scale_factor);
V_ecc_stim      = 0;      V_ecc_stim = round(V_ecc_stim*60/scale_factor);
amplitude = background*contrast/100;
%% data structure


%%%%%%%%%%%%%%%%%%%%%%%%%%

step_size = 0.2;% in log space
if block <= 1 || use_old_stair == 0
    for i = 1:length( stimulus_type )
        for j = 1:length( stair_condition )
            curr_state(i,j).duration = 0.05;
            curr_state(i,j).up_down = [0, 0];
        end
    end
else
    load( stair_name_load );
end

[stim, stair] = meshgrid( 1:length(stimulus_type), stair_condition );

condi_all = [ stim(:), stair(:) ];
condi_all = repmat( condi_all, n_repeat, 1 );
total_trials  = size( condi_all, 1 );
condi_all = condi_all( randperm( total_trials ), : );

data_all      = zeros( total_trials, 7 );% trial number, duration, direction, choice, correct, stimulus_type, stair_case
data_all(:,1) =(1:total_trials)';

direction_list = rem(randperm(total_trials),4)+1;
% direction_list = 4*ones( total_trials, 1 );


%%

try
    %open Screen windows,
    Screen('Preference', 'SkipSyncTests', 1);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens);
    %w=Screen('OpenWindow',screenNumber,0,[],8,2);
    w=Screen('OpenWindow',screenNumber,0,[],[],2);
    
    screen_rect = Screen('Rect',w);
    if linearize
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, background);Screen('Flip', w);Screen('FillRect',w, background);
    Screen('TextSize',w,30);Screen('TextFont',w,'Charcoal');
    sr_hor = round(screen_rect(3)/2); sr_ver = round(screen_rect(4)/2);
    
    
    % MAIN LOOP
    %HideCursor;
    tic; skip_first = 1;trial = 1;
    Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
    Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
    %     Screen('DrawText',w,['KEEP YOUR EYES ON THE FIXATION POINT above!'],sr_hor-300,sr_ver+50,0);
    %     Screen('DrawText',w,['press LEFT arrow for LEFTWARD motion'],sr_hor-300,sr_ver+300,0);
    %     Screen('DrawText',w,['press RIGHT arrow for RIGHTWARD motion'],sr_hor-300,sr_ver+350,0);
    %     Screen('DrawText',w,['press SPACE BAR to inititate each trial'],sr_hor-300,sr_ver+250,0);
    
    Screen('DrawText',w,[int2str(total_trials),'  trials. Press 0 key to start the experiment'],sr_hor-300,sr_ver-80,0);Screen('Flip',w);
    KbWait;
    
    GetChar;
    
    Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
    Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
    Screen('Flip', w);
    FlushEvents('keyDown');
    
    
    %housekeeping stuff
    stimulus_radius  = round(60*env_radi/scale_factor);
    
    f=(SF*scale_factor/60)*2*pi;
    orientation=orientation*pi/180;
    a=cos(orientation)*f; b=sin(orientation)*f;
    q1 = 1; q2 = 1;
    amplitude = background*contrast/100;
    
    
    % make the spatial envelope
    [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle=((stimulus_radius)^2-(x.^2+y.^2));
    for i=1:bps;
        for j =1:bps;
            if circle(i,j) < 0; circle(i,j) = 0;
            else
                circle(i,j) = 1;
            end;
        end;
    end;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/2)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;circle = (cos2D.*circle);
    end
    circleCOS = zeros(bps, bps, 2);
    circleCOS(:,:,2)=(1-circle)*white;
    circleCOS(:,:,1)=ones(bps,bps)*background;
    
    %create cosine mask to blend image below
    maskCOS=Screen('MakeTexture',w, circleCOS);
    
    %  create stimulus rectangles
    movie_rect = [0,0,bps,bps];
    scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
    scr_top = fix(screen_rect(4)/2)-round(bps/2);
    screen_rect_middle = movie_rect + [scr_left_middle, scr_top+V_ecc_fix, scr_left_middle, scr_top+V_ecc_fix];
    screen_patch = screen_rect_middle+[H_ecc_stim,V_ecc_stim,H_ecc_stim,V_ecc_stim];
    exit_index = 0;
    
    while trial <= total_trials
        
        if exit_index == 1;
            break;
        end
        
        t1 = GetSecs;
        curr_stim_index = condi_all(trial, 1);
        curr_stair = condi_all(trial, 2);
        %current_t = duration_list( t_conditions(trial) );%curr_stair_t( curr_stair_index );
        curr_stim = stimulus_type(curr_stim_index);
        current_t = curr_state( curr_stim_index, curr_stair ).duration;%curr_stair_t( curr_stair_index );
        
        %        [time_gauss,mv_length] = envelope(((time_sigmaP)/frame_rate),frame_rate,which_envelope,amplitude);
        % time_sigma = 0.5*current_t;
        [time_gauss,mv_length] = envelope_rect( current_t, frame_rate, amplitude );
        cmodulator = time_gauss/amplitude;
        %%%---------------------make mask----------------
        
        %bigPatch_X = ceil( bps + pixStep * mv_length);
        
        if mask == 1
            
            
            
            xy = diag([bps, bps]) * rand( 2, mask_n_total_dots );
            
            luminance = 2*mod( randperm( size(xy,2) ), 2 ) - 1 + background;
            color_arg = [luminance; luminance; luminance; 255*ones( 1, size(xy,2) )];
            
            Screen('DrawDots', w, xy , mask_dot_size, color_arg, [], 1);
            mask_im = Screen('GetImage', w, [0, 0, bps, bps], 'backBuffer', [], 1);
            mask_im = double(mask_im);
            %
        elseif mask == 2
            
            mask_thickness = 10;
            % rand_xy = 2*(randi( 2, ceil( bps/mask_thickness ) , 2 )-1) - 1;
            
            lum_order = ones( ceil( bps/mask_thickness ), 2 );
            lum_order( 1:2:end, : ) = 0;
            rand_xy = 2*lum_order - 1;
            
            mask_v = ones( bps, bps );
            mask_h = ones( bps, bps );
            
            for i = 1:bps
                mask_v(:,i) = rand_xy(ceil(i/mask_thickness), 1);
                mask_h(i,:) = rand_xy(ceil(i/mask_thickness), 2);
            end
            
            mask_im = 0.5*( mask_v + mask_h ) + background;
            
            
        end
        
        if mask ~= 0
            
            maskPatch = Screen('MakeTexture', w, (mask_im - background) * amplitude + background); %make texture of big matrix
            
        end
        
        
        %------------------------make the movie----------------------------
        % make a wide image patch with dots.
        frame=zeros(bps,bps,3);
        p = 1;
        if curr_stim == 3
            
            srcrect=[0 0 bps bps];
            smallPatch_X = ceil( dot_size + pixStep * mv_length );
            streak_x = 0.5*dot_size + (0:mv_length - 1)*pixStep;
            streak_y = 0.5*dot_size*ones(1,mv_length);
            
            Screen('FillRect',w, 0);
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawDots', w, [streak_x; streak_y] , dot_size, white, [], 2);
            
            white_streak = Screen('GetImage', w, [0, 0, smallPatch_X, dot_size], 'backBuffer', [], 1);
            image_streak = 255*ones( size(white_streak,1), size(white_streak,2), 2);
            image_streak(:,:,2) = white_streak;
            
            Screen('FillRect',w, background);
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            texture_streak = Screen('MakeTexture', w, image_streak );
            % try one by one first...
            
            xy = diag([bps, bps]) * rand( 2, round( n_total_dots ) );
            xy_matrix = [ xy; xy + repmat( [smallPatch_X; dot_size], 1, n_total_dots ) ];
            color_index = white*mod( randperm( n_total_dots ), 2 );
            alpha_index = 255*ones( 1, n_total_dots );
            
            color_matrix = [ color_index; color_index ; color_index; alpha_index ];
            
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawTextures', w, texture_streak, [0, 0, smallPatch_X, dot_size], xy_matrix, [], [], [], color_matrix );
            
            station_im = Screen('GetImage', w, [0.5*dot_size, 0.5*dot_size, 0.5*dot_size + bps, 0.5*dot_size + bps], 'backBuffer', [], 1);
            station_patch = Screen('MakeTexture', w, station_im );
            
            Screen('DrawTexture', w, station_patch, srcrect, screen_patch,[],[],[],[white white white white]);
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawTexture',w,maskCOS,movie_rect,screen_patch,[],[],[]);
            
            station_im_cos = Screen('GetImage',w,screen_patch,'backBuffer',[],1);
            
            for i = 1:mv_length
                im_frame = (double(station_im_cos) - background)*cmodulator(i) + background;
                if mod(i,3) == 1
                    frame(:,:,3) = im_frame;
                elseif mod(i,3) == 2
                    frame(:,:,1) = im_frame;
                elseif mod(i,3) == 0
                    frame(:,:,2) = im_frame;
                end
                if i == mv_length && mod(mv_length,3) == 1
                    frame(:,:,1) = ones(bps)*background;
                    frame(:,:,2) = ones(bps)*background;
                elseif i == mv_length && mod(mv_length,3) == 2
                    frame(:,:,2) = ones(bps)*background;
                end
                if mod(i,3) == 0 || i == mv_length
                    movie_play{p} = Screen('MakeTexture',w,frame);
                    p = p + 1;
                    frame = zeros(bps,bps,3);
                end
            end
            
            
        else
            
            bigPatch_X = ceil( bps + pixStep * mv_length);
            xy = diag([bigPatch_X, bps]) * rand( 2, round( n_total_dots * bigPatch_X/bps ) );
            
            luminance = 2*mod( randperm( size(xy,2) ), 2 ) - 1 + background;
            color_arg = [luminance; luminance; luminance; 255*ones( 1, size(xy,2) )];
            
            Screen('DrawDots', w, xy , dot_size, color_arg, [], 1);
            big_im = Screen('GetImage', w, [0, 0, bigPatch_X, bps], 'backBuffer', [], 1);
            big_im = double(big_im);
            %
            
            %create a big matrix of texture
            
            bigPatch = Screen('MakeTexture', w, (big_im - background) * amplitude + background); %make texture of big matrix
            
            
            frame=zeros(bps,bps,3);
            p=1;  %count for number of 360HZ frames
            
            if curr_stim == 1
                for i=1:mv_length
                    offset=(i-1)*pixStep;  %how many pixel moved in this frame
                    %because of left motion, source rect should move to from left
                    %to right
                    srcrect=[offset 0 bps+offset bps];
                    
                    Screen('BlendFunction', w, GL_ONE, GL_ZERO, [1,1,1,1]);
                    Screen('FillRect', w, [background, background, background, white]);
                    %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('DrawTexture',w,bigPatch,srcrect,screen_patch,[],[],[],[white white white white]);
                    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('DrawTexture',w,maskCOS,movie_rect,screen_patch,[],[],[]);
                    
                    
                    im = Screen('GetImage',w,screen_patch,'backBuffer',[],1);
                    %imageSeq(:,:,i)=im;
                    if mask~=0
                        im_mask = (mask_im - background) * amplitude + background;
                    end
                    
                    
                    % make texture out of dumped image
                    if frame_rate == 360
                        im = (double(im) - background)*cmodulator(i) + background;
                        
                        if mod(i,3) == 1
                            frame(:,:,3) = im;
                        elseif mod(i,3) == 2
                            frame(:,:,1) = im;
                        elseif mod(i,3) == 0
                            frame(:,:,2) = im;
                        end
                        if i == mv_length && mod(mv_length,3) == 1
                            frame(:,:,1) = ones(bps)*background;
                            frame(:,:,2) = ones(bps)*background;
                        elseif i == mv_length && mod(mv_length,3) == 2
                            frame(:,:,2) = ones(bps)*background;
                        end
                        if mod(i,3) == 0 || i == mv_length
                            movie_play{p} = Screen('MakeTexture',w,frame);
                            p = p + 1;
                            frame = zeros(bps,bps,3);
                        end
                    end
                end
                
            else
                
                imageSeq=zeros(bps,bps,mv_length);   %matrix contain all image
                
                for i=1:mv_length
                    offset=(i-1)*pixStep;  %how many pixel moved in this frame
                    %because of left motion, source rect should move to from left
                    %to right
                    srcrect=[offset 0 bps+offset bps];
                    
                    %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('DrawTexture',w,bigPatch,srcrect,screen_patch,[],[],[],[white white white white]);
                    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    Screen('DrawTexture',w,maskCOS,movie_rect,screen_patch,[],[],[]);
                    
                    im = Screen('GetImage',w,screen_patch,'backBuffer',[],1);
                    imageSeq(:,:,i) = double(im);
                    
                end
                
                im_station = mean(imageSeq, 3);
                
                for i = 1:mv_length
                    im_frame = (im_station - background)*cmodulator(i) + background;
                    %                        im_mask = (mask_im - background) * amplitude + background;
                    if mod(i,3) == 1
                        frame(:,:,3) = im_frame;
                    elseif mod(i,3) == 2
                        frame(:,:,1) = im_frame;
                    elseif mod(i,3) == 0
                        frame(:,:,2) = im_frame;
                    end
                    if i == mv_length && mod(mv_length,3) == 1
                        frame(:,:,1) = ones(bps)*background;
                        frame(:,:,2) = ones(bps)*background;
                    elseif i == mv_length && mod(mv_length,3) == 2
                        frame(:,:,2) = ones(bps)*background;
                    end
                    if mod(i,3) == 0 || i == mv_length
                        movie_play{p} = Screen('MakeTexture',w,frame);
                        p = p + 1;
                        frame = zeros(bps,bps,3);
                    end
                end
                
            end
        end
        
        %done with the movie
        
        %calculate the moving texture
        direction = direction_list(trial); %
        switch direction
            case 1%right
                angle = 180;
                corrKey1 = 79;  corrKey2 = 80; incorrKey1 = 81; incorrKey2 = 82;
            case 2%left
                angle = 0;
                corrKey1 = 80;	corrKey2 = 79; incorrKey1 = 81; incorrKey2 = 82;
            case 3%up
                angle = 90;
                corrKey1 = 82;	corrKey2 = 81; incorrKey1 = 79; incorrKey2 = 80;
            case 4%down
                angle = 270;
                corrKey1 = 81;	corrKey2 = 82; incorrKey1 = 80; incorrKey2 = 79;
        end
        
        
        
        %  initiate trial
        t2 = GetSecs - t1;
        
        FlushEvents('keyDown');
        WaitSecs(1-t2);
        
        Screen('FillRect',w, background);
        Screen('Flip', w);
        mm = 19;
        for i=0:4
            nn = mm-i*4;
            Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
            Screen('Flip', w);
            WaitSecs(0.05);
        end
        Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
        Screen('Flip', w);
        WaitSecs(0.36);
        
        Screen('FillRect',w, background);
        Screen('Flip', w);
        WaitSecs(0.3);
        
        
        % play the movie
        priorityLevel=MaxPriority(w);Priority(priorityLevel);
        blah=GetSecs;
        Screen('Flip',w);
        for i = 1:ceil(mv_length/3);
            Screen('DrawTexture', w, movie_play{i}, movie_rect, screen_patch, angle);
            
            % WaitSecs(1)
            Screen('Flip',w);
        end;
        
        %Screen('FillRect',w, background);
        %Screen('Flip',w);
        %WaitSecs(0.005);
        
        if mask ~= 0
            Screen('DrawTexture', w, maskPatch, srcrect, screen_patch,[],[],[],[white white white white]);
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawTexture',w,maskCOS,movie_rect,screen_patch,[],[],[]);
            Screen('Flip',w);
            
            WaitSecs(mask_duration);
        end
        
        Screen('FillRect',w, background);
        Screen('Flip',w);
        %
        
        Priority(0);
        
        % get the response && Update QUEST
        FlushEvents('keyDown');
        while (1);
            [keyIsDown,secs,keyCode] = KbCheck;
            if keyIsDown;
                if keyCode(corrKey1) || keyCode(corrKey2) || keyCode(incorrKey1) || keyCode(incorrKey2) || keyCode(41);
                    break;
                else
                    keyIsDown = 0;
                end;
            end
        end
        Screen('FillRect',w, background);
        vbl=Screen('Flip', w);
        ampS = 1;
        
        if keyCode(corrKey1) || keyCode(corrKey2)
            rs = 1;
            %Snd('Play',sin((0:1000))*ampS);
            %Snd('Wait');
        elseif keyCode(41)
            exit_index = 1;
            rs = 0;
        else
            rs=0;
        end
        %
        rs_key = find(keyCode == 1);
        data_all( trial, 2:7 ) = [current_t, direction, rs_key, rs, curr_stim, curr_stair];
        
        [new_state, new_up_down] = stair_case( log(curr_state(curr_stim_index, curr_stair).duration), curr_state(curr_stim_index, curr_stair).up_down,...
            rs, up_down_design( curr_stair, : ), step_size );
        
        curr_state(curr_stim_index, curr_stair).duration = exp( new_state );
        curr_state(curr_stim_index, curr_stair).up_down = new_up_down;
        
        trial = trial+1;
        
        FlushEvents('keyDown');
        % Close movies
        for i = 1:ceil(mv_length/3);
            % if time_gauss(i) ~= amplitude;
            Screen(movie_play{i}, 'Close');
            % end
        end
        if mask ~= 0
            Screen(maskPatch, 'Close');
        end
        clear movie motion_step;
    end
    
    if savefile~=0;
        save(filename);
        save(stair_name_save,'curr_state')
    end
    
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
catch
    ListenChar(1);
    ddd = lasterror;
    ddd.message
    ddd.stack(1,1).line
    psychrethrow(lasterror);
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
end %try..catch..
time = toc/60;


%%%%finally, do some computation of staircase,it should output thresholds
%%%%and initial values for next training blocks
