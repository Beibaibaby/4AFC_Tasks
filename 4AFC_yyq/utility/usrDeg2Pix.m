function pix = usrDeg2Pix(degree,screen_distance,screen_width) %degree为两点间视角
%根据观察距离和屏幕宽度将视角degree换算为pix（像素）
screens=Screen('Screens');  screen_num=max(screens); %选择编号�?大的窗口 
[Width]=Screen('WindowSize',screen_num);
pix=tan(degree./2./180.*pi)./(screen_width./2./screen_distance).*Width; %点距视角在全屏视角中占比*像素宽度
end