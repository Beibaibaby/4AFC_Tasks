function [mask_texid,mask_img]=get_mask1(w,wrect,bgc,R)

Screen('Flip',w);

Screen('FillRect',w,255);

[X,Y]=meshgrid(-R:R-1,-R:R-1);

XY=((X.^2)+(Y.^2)).^(0.5);

coss=cos(XY./(R/pi*2));

coss=1-coss;
coss=coss.*255;
coss=uint8(coss);

texid=Screen('MakeTexture',w,coss);
Screen('DrawTexture',w,texid);
img=Screen('GetImage',w,wrect,'backbuffer');

mask_img=ones(size(img)).*bgc;
mask_img(:,:,4)=img(:,:,1);

mask_texid=Screen('MakeTexture',w,mask_img);

Screen('FillRect',w,bgc);
Screen('Flip',w);




