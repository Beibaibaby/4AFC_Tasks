function draw_center_text(w,str,x,y,color,strsize)

if nargin<=4
    color=[];
    strsize=[];
elseif nargin<=5
    strsize=[];
end

if isempty(color)
    color=0;
end
if isempty(strsize)
    strsize=30;
end

old_size=Screen('TextSize',w,strsize);
str=double(str);
strRect=Screen('TextBounds',w,str);
sx=x-strRect(3)/2;
sy=y-strRect(4)/2;
Screen('DrawText',w,str,sx,sy,color);
Screen('TextSize',w,old_size);
