
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
subject_initials    ='OK';
block               = 10;  % 0 means test,run 10 blocks total
%duration_list        = [0.005, 0.03, 0.04, 0.06, 0.08, 0.16]; %KP
%duration_list        = [0.005, 0.02, 0.03, 0.045, 0.06, 0.1]; %PF
duration_list        = [0.005, 0.02, 0.03, 0.045, 0.06, 0.1]; %LS, and OK low
%duration_list       = logspace(log10(dura_range(1)), log10(dura_range(2)), 12);%dura_range          =[0.050 0.050];%OK %seconds, change the duration range here based on individual psychmetric function
%dura_range          =[0.005 0.065];%OK %low contrast.
%dura_range          =[0.1 0.1];%for Olga.
%dura_range          =[0.005 0.2];%for KS.
%dura_range          = [0.04 0.04];%RZ,high contrast
%size                = 0.7;    %deg, radius
env_radi            = 2;    %2, 11deg, radius
contrast            = 3;   %3, 100%
speed               = 5;  %deg/sec


%%


savefile            =1;
tme                 =clock;
filename = strcat(subject_initials,'_block',num2str(block),'_',num2str(speed),'.mat');
IsExist = exist(filename,'file');
if IsExist && block ~= 0
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
theX                   = -640; 
theY                   = 360;

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
%%

%data structure
total_trials  = 30*length(duration_list);
%total_trials  =10;
data_all      =zeros(total_trials,4);% trial number,duration,direction,choice,correct
data_all(:,1) =(1:total_trials)';
t_conditions  = rem(randperm(total_trials),length(duration_list))+1;
%direction_list= rem(randperm(total_trials),4)+1;
direction_list = round((1:total_trials)*360/total_trials);
direction_list = direction_list( randperm( total_trials ) );

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
    % make arrow 
    
    load('arrow_70_140.mat');
    for i = 1:3
        arrow_matrix(:,:,i) = rot90(arrow_70_140(:,:,i), 3);
    end
    arrow_matrix( arrow_matrix == 127 ) = background;
    arrow_ext = background*ones( size( arrow_matrix ) );
    arrow_ext = uint8( arrow_ext(:,1:end-10,:) );
    arrow_matrix = [ arrow_matrix, arrow_ext ];
    arrow_height = size(arrow_matrix, 1);
    arrow_width = size(arrow_matrix, 2);
    arrow_rect = [ 0, 0, arrow_width, arrow_height ];
    arrow_left_middle = screen_rect(3)/2 - arrow_width/2;
    arrow_top_middle = screen_rect(4)/2 - arrow_height/2;
    arrow_patch = arrow_rect + [arrow_left_middle, arrow_top_middle, arrow_left_middle, arrow_top_middle];
    arrow = Screen('MakeTexture', w, arrow_matrix);
    
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
    HideCursor
    while trial <= total_trials
        
        if exit_index == 1;
            ShowCursor
            break;
        end
        
        t1 = GetSecs;
        current_t = duration_list( t_conditions(trial) );%curr_stair_t( curr_stair_index );
        
        %        [time_gauss,mv_length] = envelope(((time_sigmaP)/frame_rate),frame_rate,which_envelope,amplitude);
        time_sigma = 0.5*current_t;
        [time_gauss,mv_length] = envelope_rect( current_t, frame_rate, amplitude );
        cmodulator = time_gauss/amplitude;
        
        %------------------------make the movie----------------------------
        % make a wide image patch with dots.
        frame=zeros(bps,bps,3);
        p = 1;
        
        
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
            %imageSeq(:,:,i)=im;
            
            
            
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
        
        
        
        %done with the movie
        
        %calculate the moving texture
        direction = direction_list(trial); %
        
        
        corrKey = 82;	incorrKey1 = 80; incorrKey2 = 81; incorrKey3 = 79;
        
        
        
        
        %  initiate trial
        t2 = GetSecs - t1;
        
        FlushEvents('mouseUp','mouseDown','keyDown');
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
        FlushEvents('mouseUp','mouseDown','keyDown');
        SetMouse(theX,theY);
        
        for i = 1:ceil(mv_length/3);
            Screen('DrawTexture', w, movie_play{i}, movie_rect, screen_patch, direction);
            
 %            WaitSecs(1)
            Screen('Flip',w);
        end;
        
        Screen('FillRect',w, background);
        Screen('Flip',w);
        
        Priority(0);
        
        FlushEvents('keyDown');
        while 1
            [x,y,buttons] = GetMouse;
            if x ~= theX || y ~= theY
                break
            end
        end
        
        while 1
            
            [keyIsDown,secs,keyCode] = KbCheck;
            if keyIsDown
                exit_index = 1;
                break
            end
            
            if any(buttons) % wait for press
                break
            end
            
            arrow_direction = cart2pol( x - theX, y - theY )*180/pi;
            arrow_direction = ( arrow_direction - 180 );
            if arrow_direction < 1
                arrow_direction = arrow_direction + 360;
            end
            
            Screen('DrawTexture', w, arrow, arrow_rect, arrow_patch, arrow_direction);
            Screen('DrawDots', w, [x+1280; y], 20, 200, [], 2);
            Screen('Flip',w);
            [x,y,buttons] = GetMouse;
        end
        
        Screen('FillRect',w, background);
        vbl=Screen('Flip', w);
        ampS = 1;
        
        data_all( trial, 2:4 ) = [current_t, direction, arrow_direction ];
        
        
        trial = trial+1;
        
        FlushEvents('keyDown');
        % Close movies
        for i = 1:ceil(mv_length/3);
            Screen(movie_play{i}, 'Close');
        end
        
        clear movie motion_step;
    end
    ShowCursor
    if savefile~=0;
        save(filename);
        %     save(stair_file_name_save,'curr_stair_t')
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
