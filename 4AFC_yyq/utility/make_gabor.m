function gabor = make_gabor(size_x,size_y,pix_per_cycle,orientation,phase,color1,color2)
%size， gabor的尺寸， 直径
%pix_per_cycle , 每个周期的像素数量
%orientation， gabor的朝向，0表示水平，该值增加，则gabor逆时针旋转，单位为°
% phase，相位
% color1  峰值颜色
% color2  谷值颜色

color1=color1/255;
color2=color2/255;

[X,Y]=meshgrid(1:size_x,1:size_y);%创建size_y*size_x的网格，Y为1:size_y列向量横向重复size_x次

a=(sin(  (sind(orientation)*X+cosd(orientation)*Y)  *  2*pi/pix_per_cycle  +phase  )+1)/2;

a1=mapminmax(a,color1(1),color2(1));%mapminmax使a归一化至峰、谷值内
a2=mapminmax(a,color1(2),color2(2));
a3=mapminmax(a,color1(3),color2(3));

gabor=cat(3,a1,a2,a3);%堆叠a1 a2 a3三个颜色矩阵
gabor=im2uint8(gabor);%将gabor数组化为8bit类图像
gabor(:,:,4)=128;%透明度50%？

