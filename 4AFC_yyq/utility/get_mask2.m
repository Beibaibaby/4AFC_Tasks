function [mask_texid,mask_img]=get_mask2(w,bgc)
%2.27 卢毓昕 绘制全屏矩形背景色的texture
Screen('Flip',w);
Screen('FillRect',w,bgc);
mask_img=Screen('GetImage',w,[],'backbuffer');
mask_texid=Screen('MakeTexture',w,mask_img);
Screen('FillRect',w,bgc);
Screen('Flip',w);

