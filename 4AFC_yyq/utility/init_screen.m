function [w,wrect,hz,xc,yc]=init_screen(bgc,skip,exp_hz)

if nargin<3 %nargin：输入参数数量
    exp_hz=[];
end

Screen('CloseAll');
Screen('Preference','SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', skip); %skip=1则跳过同步性能力测试
screens=Screen('Screens');  screen_num=max(screens);  
[w,wrect]=Screen('OpenWindow',screen_num,bgc);
% [w,wrect]=PsychImaging('OpenWindow', 0, bgc, [], [], [], [], 16);
ListenChar(2) ;  %关闭编辑器对键盘监听
HideCursor;
hz=FrameRate(w);
hz=round(hz);

if ~isempty(exp_hz)
    if hz~=exp_hz
        sca;
        error('刷新率与设定值不符')
    end
end

[xc,yc]=WindowCenter(w);
AssertOpenGL; %OpenGL：开放图形库，作用是？
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Priority(MaxPriority(w));
KbName('UnifyKeyNames');
Screen('TextSize',w,55);
Screen('TextFont',w,'-:lang=zh-cn');
rng('Shuffle');%根据当前时间初始化随机数生成器