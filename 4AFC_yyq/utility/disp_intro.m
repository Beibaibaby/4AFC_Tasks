function disp_intro(w,array_name,key)
if nargin<3
    key=32;
end
array=imread(array_name);
texid=Screen('MakeTexture',w,array);
while 1
    Screen('DrawTexture',w,texid);
    Screen('Flip',w);
    [~,~,kc]=KbCheck;
    if kc(key)
        break;
    end
    checkend
end
Screen('Close',texid);

while 1
    [~,~,kc]=KbCheck;
    if kc(key)==0
        break
    end
end