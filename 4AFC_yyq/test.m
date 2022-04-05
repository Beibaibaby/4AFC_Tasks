%2022.3.29
%不同时长

%%
clear;clc;sca;
rng('Shuffle');%根据当前时间初始化随机数生成器
addpath('utility')

%% 可调参数
screen_distance=60; %被试与屏幕的距离（需要预先测量）  厘米
screen_width=50; %屏幕水平宽度 （需要预先测量） 厘米
exp_hz=60; %实验电脑的刷新率
bgc=128; %背景颜色 rgb值

grating_type=2; %光栅交叉的方式， 1表示 水平*垂直 ，   2表示 左斜*右斜
move_speed=4; %运动速度，deg/sce
cpd=1;  %每度的周期数量%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ori=[1,2,3,4]; %运动方向：左、右、上、下
stimulus=1; % 1:高对比度大刺激；2：低对比度小刺激
duration=[20,40,60,80,100,120,140]; %运动时长

time1=0.5; %注视点时长
time2=0.5; %试次结束后的空屏

n=25;   %每个方向重复次数
m=size(ori,2)*n;  %每个时长重复次数
trial_num=m*length(duration); %试次数量

%%
if stimulus==1
    R=10; %半径 deg :  [1,10]
    contrast=0.99;  %对比度 : [0.05,0.99]
else
    R=1;
    contrast=0.05;
end
R = usrDeg2Pix(R,screen_distance,screen_width) ;
ppc = usrDeg2Pix(1,screen_distance,screen_width) / cpd; %每个周期的像素数量
move_speed_per_frame=usrDeg2Pix(move_speed,screen_distance,screen_width)/exp_hz;  %每帧的速度
if grating_type==1
    orientation1=0; %水平
    orientation2=90; %垂直
else
    orientation1=45;
    orientation2=135;
end

%% trials
durs=reshape(repmat(duration,m,1),trial_num,1);
oris=reshape(repmat(repmat(ori,1,n),1,size(duration,2)),trial_num,1);
trials=Shuffle([durs,oris],2);

%% window
sub=inputdlg({'被试编号'});
subID=sub{1};
clock0 = fix(clock);
f = sprintf('%04d%02d%02d%02d%02d%02d',clock0);
expname=[pwd,'/data/',f,'test_subID-',subID,'.mat'] ;

[w,wrect,hz,xc,yc]=init_screen(bgc,2,exp_hz);%打开窗口（颜色，skip，预设刷新率）
disp_intro(w,'instruction.png',13);

%% test_cycle
for trial=1:trial_num 
    this_ori=trials(trial,2);%本次运动方向
    time_sigma=trials(trial,1)/1000;%本次运动时长
    
    [time_Gauss,mv_length] = envelope(time_sigma/2,exp_hz,2,bgc);
    time_Gauss=time_Gauss./bgc;
    thismove_dist=move_speed_per_frame*mv_length;  %运动距离
    
    if this_ori==1 %向左运动
        gabor_size_y=R*2;  %光栅直径
        gabor_size_x=thismove_dist+gabor_size_y;
        init_dot=[xc-R,yc-R];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[-move_speed_per_frame,0,-move_speed_per_frame,0];
    elseif this_ori==2 %向右运动
        gabor_size_y=R*2;  %光栅直径
        gabor_size_x=thismove_dist+gabor_size_y;
        init_dot=[xc+R-gabor_size_x,yc-R];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[move_speed_per_frame,0,move_speed_per_frame,0];
    elseif this_ori==3 %向上运动
        gabor_size_x=R*2;
        gabor_size_y=thismove_dist+gabor_size_x;
        init_dot=[xc-R,yc-R];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[0,-move_speed_per_frame,0,-move_speed_per_frame];
    elseif this_ori==4 %向下运动
        gabor_size_x=R*2;
        gabor_size_y=thismove_dist+gabor_size_x;
        init_dot=[xc-R,yc+R-gabor_size_y];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[0,move_speed_per_frame,0,move_speed_per_frame];
    end
    
    [mask_texid,mask_img]=get_mask1(w,wrect,bgc,R); %获取本次的mask
    
    this_contrast=contrast;
    [H,L]=get_HL(this_contrast); %依据Michelson对比度公式计算对比度
    
    gabor1= make_gabor(gabor_size_x,gabor_size_y,ppc,orientation1,2*pi*rand,[H,H,H],[L,L,L]);
    gabor2= make_gabor(gabor_size_x,gabor_size_y,ppc,orientation2,2*pi*rand,[H,H,H],[L,L,L]);
    
    texid1=Screen('MakeTexture',w,gabor1);
    texid2=Screen('MakeTexture',w,gabor2);
    
    Screen('DrawDots',w,[xc,yc],10,0,[],1);%fixation(w窗口，居中，大小，颜色，反锯齿圆type)
    Screen('Flip',w);
    WaitSecs(time1);
    for flip=1:mv_length %通过连续flip实现画面运动
        Screen('DrawTexture',w,texid1,[],init_rect,[],[],time_Gauss(flip));
        Screen('DrawTexture',w,texid2,[],init_rect,[],[],time_Gauss(flip));
        Screen('DrawTexture',w,mask_texid);
        Screen('Flip',w);
        init_rect=init_rect+speed_rect;
        checkend;
    end
    Screen('Flip',w);
    resstart=GetSecs;
    while 1
        checkend;
        [kid,~,kc]=KbCheck;
        if kid
            key=find(kc);
            [inx,subres]=ismember(key,[37,39,38,40]);
            if inx==0
                continue
            end
            RT=GetSecs-resstart;
            if subres==this_ori
                judge=1;
            else
                judge=0;
            end
            break;
        end
    end
    
    Screen('Close',mask_texid);
    Screen('Close',texid1);
    Screen('Close',texid2);
    
    %第一列，运动方向， 第二列，半径， 第三列对比度，第四列运动时长，第五列无意义
    result(trial).trial=trial;
    result(trial).this_ori=this_ori; %运动方向
    result(trial).this_R=R; %半径
    result(trial).this_cont=this_contrast; %对比度
    result(trial).this_dur=time_sigma; %运动时长
    result(trial).this_mv_dur=mv_length/exp_hz; %真实运动时长
    result(trial).time_Gauss=time_Gauss; %时间轴上的透明度变化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    result(trial).subres=subres; %被试反应
    result(trial).RT=RT; %反应时
    result(trial).judge=judge; %正确性
    Screen('Flip',w);
    WaitSecs(time2);
    save(expname,'result');
    
    if mod(trial,100)==0
        draw_center_text(w,'Break. Press "SPACE" to continue.',wrect(3)/2,wrect(4)/2,0,30);
        Screen('Flip',w);
        while 1
            checkend;
            [kid,~,kc]=KbCheck;
            if kid
                key=find(kc);
                if key==32
                    break;
                end
            end
        end
    end
end
%%
save(expname);
ListenChar(0);
ShowCursor;
sca;
return


