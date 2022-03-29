% this code is to measure suppression and summation curve using FAST
% 2 AFC and 4 AFC here,but we use plaids stimuli 


% We use both 4 AFC and 2 AFC here 
% suppression curve.
%04/21/14 Ruyuan Zhang change it for data collection


clc;clear all;close all;

ListenChar(2);
warning('off','MATLAB:dispatcher:InexactMatch')

% KEY Monitor Parameters
scale_factor   = 2.07;       % most important parameter - how many acrmin is one screen pixel? for SONY monitor at [1024 640] resolution, and 30.4in vieving distance, use scale_factor=2
frame_rate     = 360;        % for the actual exp: SET THIS AND MONITOR TO 120!!!!!!!
resolution     = [1280 720]; % for the actual exp: SET THIS AND MONITOR to [1024 640]
linearize 	   = 1;          % whether monitor is linearized. if=1, program will look for "MyGammaTable"
break_if_mismatch = 1;       % if 1, this will quit the program if monitor resolution does not match above settings (use 1 for the ACTUAL EXP, 0 for testing on a different monitor, eg laptop)


%% -----you may check parameter here every time you run-------------------
subject_initials    ='LS';
contrast            = 100;    %5, 100
%% -------------------

%%%%%% I SHOULD USE SMALLER TARGET SIZE NEXT TIME LIKE 0.5.... AND I CAN
%%%%%% GIVE FEEDBACK
if contrast <= 15
    sizes           = exp(linspace( log(1), log(11), 6 ));
else
    %sizes           = exp(linspace( log(0.7), log(9), 6 ));%notice the size setting here
    sizes           = exp(linspace( log(0.5), log(11), 6 ));%notice the size setting here
end

%sizes = 11;

No_contrast         = 1;
main_condts         = length(sizes);
stimulus_radius_all =sizes*60/scale_factor;

% housekeeping
Hor_eccentr         = 0;
Ver_eccentr         = 0;
V_ecc_fix           = 0;
H_ecc_fix           = 0;
V_ecc_fix = V_ecc_fix*60/scale_factor;
H_ecc_fix = H_ecc_fix*60/scale_factor;


% upper_limit   = 150;    % ms
tme = clock;
upper_bound = 0.5;
l = 7;
ww = 4;
i_result = 1;
TT1 = GetSecs;
k = 1;


%% ---stimulus parameter---
speed                 = 6;  %deg/sec
spatial_envelope      = 2;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
which_envelope 	      = 2;
background            = 126;   % background intensity, in gray scale untis
white                 = background*2;
SF                    = 1.2;          %cycles/deg
TF 					  = speed*SF;     %cycles/sec
orientation           = 0;   %deg
amplitude             = background*contrast/100;

%% data structure
total_trials          =600;
data_all              = zeros( total_trials,7); %trial number,size,duration,direcition,rs_key,rs,fast_number
data_all(:,1)         =(1:total_trials)';
fast_list             =rem(randperm(total_trials),2)+1;
%fast_list             =ones(1,total_trials)*2;
direction_list        =rem(randperm(total_trials),4)+1;
exp_order             = rem( randperm(total_trials), main_condts) + 1; %size list
%%



try
    %open Screen windows
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
    
    
    % 
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
    
  
    %set up fast
    fastSS1 = fastFull( 4,...
        'funcSS',...
        'psyWeibull',...
        { [0.000001, 0.3], [0.000001, 0.3], [0.000001, 4], [0.0002, 20] } ); %
    
    fastSS2 = fastFull( 2,...
        'funcSS',...
        'psyWeibull',...
        { [0.000001, 0.3], [0.000001, 0.3], [0.000001, 4], [0.0002, 20] } );%
 
    
    while trial <= total_trials
        
        t1 = GetSecs;

        stim_condt  = exp_order(trial);
        current_r   = sizes(stim_condt); % radius in degree
        H_ecc_stim  = Hor_eccentr*60/scale_factor;
        V_ecc_stim  = Ver_eccentr*60/scale_factor;
        
        if fast_list(trial)==1
            n_alternatives = 4;
            current_t = fastChooseY( fastSS1, current_r ); % duration in second.
        elseif fast_list(trial)==2
            n_alternatives = 2;
            current_t = fastChooseY( fastSS2, current_r ); % duration in second.
        end
        
        
        if current_r < 7
            if current_t > upper_bound
                current_t = upper_bound;
                
            end
        else
            if current_t > 0.3         % 100ms as the boundary
                current_t = 0.3;
                
            end
        end
        
        %housekeeping stuff
        TFstep = (2*pi*TF)/frame_rate;
        f=(SF*scale_factor/60)*2*pi;
        orientation=orientation*pi/180;
        a=cos(orientation)*f; b=sin(orientation)*f;
        q1 = 1; q2 = 1;
        
        
        %calculate the moving texture
        direction = direction_list(trial); %
        if n_alternatives==4
            switch direction
                case 1%right
                    angle = 135;
                    mask_angle=90;
                    corrKey = 79;	incorrKey1 = 80; incorrKey2 = 81; incorrKey3 = 82;
                case 2%left
                    angle = 315;
                    mask_angle=90;
                    corrKey = 80;	incorrKey1 = 79; incorrKey2 = 81; incorrKey3 = 82;
                case 3%down
                    angle = 225;
                    mask_angle=0;
                    corrKey = 81;	incorrKey1 = 79; incorrKey2 = 80; incorrKey3 = 82;
                case 4%up
                    angle = 45;
                    mask_angle=0;
                    corrKey = 82;	incorrKey1 = 79; incorrKey2 = 80; incorrKey3 = 81;                      
            end
        elseif n_alternatives==2
            switch direction
                case 1%right
                    angle = 135;
                    mask_angle=90;
                    corrKey1 = 79;	corrKey2 = 80; incorrKey1 = 81; incorrKey2 = 82;
                case 2%left
                    angle = 315;
                    mask_angle=90;
                    corrKey1 = 80;	corrKey2 = 79; incorrKey1 = 81; incorrKey2 = 82;
                case 3%down
                    angle = 225;
                    mask_angle=0;
                    corrKey1 = 81;	corrKey2 = 82; incorrKey1 = 79; incorrKey2 = 80;
                case 4%up
                    angle = 45;
                    mask_angle=0;
                    corrKey1 = 82;	corrKey2 = 81; incorrKey1 = 80; incorrKey2 = 79;                      
            end
        end
        
        
        %% make the spatial envelope
        stimulus_radius = stimulus_radius_all(stim_condt);
        [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
        bps = length(x);%round((stimulus_radius)*2+1);
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
        

        %create stimulus rectangles
        movie_rect = [0,0,bps,bps];
        scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
        scr_top = fix(screen_rect(4)/2)-round(bps/2);
        screen_rect_middle = movie_rect + [scr_left_middle, scr_top+V_ecc_fix, scr_left_middle, scr_top+V_ecc_fix];
        screen_patch = screen_rect_middle+[H_ecc_stim,V_ecc_stim,H_ecc_stim,V_ecc_stim];
        
        %% ------------------------make the movie----------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        %        [time_gauss,mv_length] = envelope(((time_sigmaP)/frame_rate),frame_rate,which_envelope,amplitude);
        time_sigma = 0.5*current_t;
        [time_gauss,mv_length] = envelope( time_sigma, frame_rate, which_envelope, amplitude );
        time_gauss = time_gauss/amplitude;
        
        %% -----------------------make the movie----------------------------
        %make the movie
        motion_step(1) = rand*2*pi;
        for i=2:mv_length;
            motion_step(i) = motion_step(i-1)+TFstep;
        end

        for i = 1:mv_length;
            %first grating
            grating1=round(((sin(a*x+b*y+motion_step(i)).*circle*amplitude*time_gauss(i))+background));
            %second grating,rotate 90 deg
            grating2=grating1';
            movie{i}=(grating1+grating2)/2;
            %movie{i}=grating1;
        end
        
        frame = zeros(bps,bps,3);
        for i = 1:ceil(mv_length/3)
            for j=1:3
                if ((i-1)*3+j)>mv_length
                    %frame(:,:,j) = ones(bps)*background;
                    switch j
                        case 1
                            frame(:,:,3) = ones(bps)*background;
                        case 2
                            frame(:,:,1) = ones(bps)*background;
                        case 3
                            frame(:,:,2) = ones(bps)*background;
                    end
                else
                    switch j
                        case 1
                            frame(:,:,3) = movie{(i-1)*3+j};
                        case 2
                            frame(:,:,1) = movie{(i-1)*3+j};
                        case 3
                            frame(:,:,2) = movie{(i-1)*3+j};
                    end
                end
            end
            movie_play{i} = Screen('MakeTexture',w,frame);
        end
        
        %done with the movie
        %% play the movie
       
      
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
            %Screen('DrawTexture', w, movie_play{i}, movie_rect, screen_patch, angle);
            
            % WaitSecs(1)
            Screen('Flip',w);
        end;

        Screen('FillRect',w, background);
        Screen('Flip',w);

        Priority(0);
        
        %%  get the response && Update QUEST
        FlushEvents('keyDown');
        while (1);
            [keyIsDown,secs,keyCode] = KbCheck;
            if keyIsDown;
                if keyCode(79) || keyCode(80) || keyCode(81) || keyCode(82);
                    break;
                else
                    keyIsDown = 0;
                end;
            end
        end
        Screen('FillRect',w, background);
        vbl=Screen('Flip', w);
        ampS = 1;
        if n_alternatives==2
            if keyCode(corrKey1)||keyCode(corrKey2);
                rs = 1;
                %Snd('Play',sin((0:1000))*ampS);
                %Snd('Wait');
            else
                rs = 0;
            end
        elseif n_alternatives==4
            if keyCode(corrKey);
                rs = 1;
                %Snd('Play',sin((0:1000))*ampS);
                %Snd('Wait');
            else
                rs = 0;
            end
        end
        
        
        rs_key = find(keyCode == 1);
        data_all( trial, 2:7 ) = [current_r, current_t, direction, rs, rs_key(1),fast_list(trial)];
        if fast_list(trial)==1
            [fastSS1, resample] = fastUpdate( fastSS1, [current_r, current_t, rs] );
        elseif fast_list(trial)==2
            [fastSS2, resample] = fastUpdate( fastSS2, [current_r, current_t, rs] );
        end
        
        
        FlushEvents('keyDown');
        
        
        %flush video card memory
        for i = 1:ceil(mv_length/3);
            % if time_gauss(i) ~= amplitude;
            Screen(movie_play{i}, 'Close');
            % end
        end
        clear movie moving_grattingP;
        
        %set a break
        if rem(trial,200)==0
            Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
            Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
            Screen('DrawText',w,[int2str(trial),'  trials. Press 0 key to continue'],sr_hor-300,sr_ver-80,0);Screen('Flip',w);
            KbWait;
            GetChar;
            FlushEvents('keyDown');
        end
        
    trial = trial+1;
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


data_file_name = strcat(subject_initials,'_curves',num2str(contrast),'.mat');

IsExist = exist(data_file_name,'file');

%if IsExist && block ~= 0
if IsExist ~= 0
    error('data file name exists')
end

save(data_file_name);
