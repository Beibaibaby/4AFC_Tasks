%2022.3.29
%��ͬʱ��

%%
clear;clc;sca;
rng('Shuffle');%���ݵ�ǰʱ���ʼ�������������
addpath('utility')

%% �ɵ�����
screen_distance=60; %��������Ļ�ľ��루��ҪԤ�Ȳ�����  ����
screen_width=50; %��Ļˮƽ��� ����ҪԤ�Ȳ����� ����
exp_hz=60; %ʵ����Ե�ˢ����
bgc=128; %������ɫ rgbֵ

grating_type=2; %��դ����ķ�ʽ�� 1��ʾ ˮƽ*��ֱ ��   2��ʾ ��б*��б
move_speed=4; %�˶��ٶȣ�deg/sce
cpd=1;  %ÿ�ȵ���������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ori=[1,2,3,4]; %�˶��������ҡ��ϡ���
stimulus=1; % 1:�߶Աȶȴ�̼���2���ͶԱȶ�С�̼�
duration=[20,40,60,80,100,120,140]; %�˶�ʱ��

time1=0.5; %ע�ӵ�ʱ��
time2=0.5; %�Դν�����Ŀ���

n=25;   %ÿ�������ظ�����
m=size(ori,2)*n;  %ÿ��ʱ���ظ�����
trial_num=m*length(duration); %�Դ�����

%%
if stimulus==1
    R=10; %�뾶 deg :  [1,10]
    contrast=0.99;  %�Աȶ� : [0.05,0.99]
else
    R=1;
    contrast=0.05;
end
R = usrDeg2Pix(R,screen_distance,screen_width) ;
ppc = usrDeg2Pix(1,screen_distance,screen_width) / cpd; %ÿ�����ڵ���������
move_speed_per_frame=usrDeg2Pix(move_speed,screen_distance,screen_width)/exp_hz;  %ÿ֡���ٶ�
if grating_type==1
    orientation1=0; %ˮƽ
    orientation2=90; %��ֱ
else
    orientation1=45;
    orientation2=135;
end

%% trials
durs=reshape(repmat(duration,m,1),trial_num,1);
oris=reshape(repmat(repmat(ori,1,n),1,size(duration,2)),trial_num,1);
trials=Shuffle([durs,oris],2);

%% window
sub=inputdlg({'���Ա��'});
subID=sub{1};
clock0 = fix(clock);
f = sprintf('%04d%02d%02d%02d%02d%02d',clock0);
expname=[pwd,'/data/',f,'test_subID-',subID,'.mat'] ;

[w,wrect,hz,xc,yc]=init_screen(bgc,2,exp_hz);%�򿪴��ڣ���ɫ��skip��Ԥ��ˢ���ʣ�
disp_intro(w,'instruction.png',13);

%% test_cycle
for trial=1:trial_num 
    this_ori=trials(trial,2);%�����˶�����
    time_sigma=trials(trial,1)/1000;%�����˶�ʱ��
    
    [time_Gauss,mv_length] = envelope(time_sigma/2,exp_hz,2,bgc);
    time_Gauss=time_Gauss./bgc;
    thismove_dist=move_speed_per_frame*mv_length;  %�˶�����
    
    if this_ori==1 %�����˶�
        gabor_size_y=R*2;  %��դֱ��
        gabor_size_x=thismove_dist+gabor_size_y;
        init_dot=[xc-R,yc-R];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[-move_speed_per_frame,0,-move_speed_per_frame,0];
    elseif this_ori==2 %�����˶�
        gabor_size_y=R*2;  %��դֱ��
        gabor_size_x=thismove_dist+gabor_size_y;
        init_dot=[xc+R-gabor_size_x,yc-R];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[move_speed_per_frame,0,move_speed_per_frame,0];
    elseif this_ori==3 %�����˶�
        gabor_size_x=R*2;
        gabor_size_y=thismove_dist+gabor_size_x;
        init_dot=[xc-R,yc-R];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[0,-move_speed_per_frame,0,-move_speed_per_frame];
    elseif this_ori==4 %�����˶�
        gabor_size_x=R*2;
        gabor_size_y=thismove_dist+gabor_size_x;
        init_dot=[xc-R,yc+R-gabor_size_y];
        init_rect=[init_dot,init_dot+[gabor_size_x,gabor_size_y]];
        speed_rect=[0,move_speed_per_frame,0,move_speed_per_frame];
    end
    
    [mask_texid,mask_img]=get_mask1(w,wrect,bgc,R); %��ȡ���ε�mask
    
    this_contrast=contrast;
    [H,L]=get_HL(this_contrast); %����Michelson�Աȶȹ�ʽ����Աȶ�
    
    gabor1= make_gabor(gabor_size_x,gabor_size_y,ppc,orientation1,2*pi*rand,[H,H,H],[L,L,L]);
    gabor2= make_gabor(gabor_size_x,gabor_size_y,ppc,orientation2,2*pi*rand,[H,H,H],[L,L,L]);
    
    texid1=Screen('MakeTexture',w,gabor1);
    texid2=Screen('MakeTexture',w,gabor2);
    
    Screen('DrawDots',w,[xc,yc],10,0,[],1);%fixation(w���ڣ����У���С����ɫ�������Բtype)
    Screen('Flip',w);
    WaitSecs(time1);
    for flip=1:mv_length %ͨ������flipʵ�ֻ����˶�
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
    
    %��һ�У��˶����� �ڶ��У��뾶�� �����жԱȶȣ��������˶�ʱ����������������
    result(trial).trial=trial;
    result(trial).this_ori=this_ori; %�˶�����
    result(trial).this_R=R; %�뾶
    result(trial).this_cont=this_contrast; %�Աȶ�
    result(trial).this_dur=time_sigma; %�˶�ʱ��
    result(trial).this_mv_dur=mv_length/exp_hz; %��ʵ�˶�ʱ��
    result(trial).time_Gauss=time_Gauss; %ʱ�����ϵ�͸���ȱ仯%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    result(trial).subres=subres; %���Է�Ӧ
    result(trial).RT=RT; %��Ӧʱ
    result(trial).judge=judge; %��ȷ��
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


