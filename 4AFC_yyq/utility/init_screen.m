function [w,wrect,hz,xc,yc]=init_screen(bgc,skip,exp_hz)

if nargin<3 %nargin�������������
    exp_hz=[];
end

Screen('CloseAll');
Screen('Preference','SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', skip); %skip=1������ͬ������������
screens=Screen('Screens');  screen_num=max(screens);  
[w,wrect]=Screen('OpenWindow',screen_num,bgc);
% [w,wrect]=PsychImaging('OpenWindow', 0, bgc, [], [], [], [], 16);
ListenChar(2) ;  %�رձ༭���Լ��̼���
HideCursor;
hz=FrameRate(w);
hz=round(hz);

if ~isempty(exp_hz)
    if hz~=exp_hz
        sca;
        error('ˢ�������趨ֵ����')
    end
end

[xc,yc]=WindowCenter(w);
AssertOpenGL; %OpenGL������ͼ�ο⣬�����ǣ�
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Priority(MaxPriority(w));
KbName('UnifyKeyNames');
Screen('TextSize',w,55);
Screen('TextFont',w,'-:lang=zh-cn');
rng('Shuffle');%���ݵ�ǰʱ���ʼ�������������