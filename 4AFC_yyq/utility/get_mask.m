function [mask_texid,mask_img]=get_mask(w,xc,yc,bgc,R,sigma)
%2.27 卢毓昕 已设置虚边缘，可考虑在主程序中输入实边缘再换算赋给R
%仅适用于视窗垂直且水平居中于mask的情况！！！否则需调整二维高斯分布函数
%双重for循环可能增大运算量，考虑前置get_mask算好后面直接调用图像？
%sigma= %单位：像素
Screen('Flip',w);
Screen('FillRect',w,bgc);
Screen('FillOval',w,0,[xc-R,yc-R,xc+R,yc+R]);
mask_img=Screen('GetImage',w,[],'backbuffer');
inx=mask_img(:,:,1)~=0;%圆形视窗外的像素赋为1，视窗内赋为0，得到1080x1920的logic数组
mask=mask_img(:,:,3);%
mask(inx)=255;%将255赋给logic=1的所有位置（视窗外）
mask_img(:,:,4)=mask;%mask_img的视窗外区域皆为不透明

[I,J]=size(mask_img(:,:,4));
mu=round(R);
x=1:(2*mu);
pd=makedist('Normal','mu',mu,'sigma',sigma);
y=pdf(pd,x);
y=255*y/y(mu);%透明度区间：0~255
% y(1)=0;

for i=1:I
    for j=1:J
        if mask_img(i,j,4)==0
            mask_img(i,j,4)=y(1+round(   ((i-round(I/2))^2+(j-round(J/2))^2) ^0.5));
            %根据半径给透明度赋值
        end
    end
end

mask_img(:,:,1:3)=bgc;

mask_texid=Screen('MakeTexture',w,mask_img);
Screen('FillRect',w,bgc);
Screen('Flip',w);

