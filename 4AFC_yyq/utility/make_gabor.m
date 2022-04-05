function gabor = make_gabor(size_x,size_y,pix_per_cycle,orientation,phase,color1,color2)
%size�� gabor�ĳߴ磬 ֱ��
%pix_per_cycle , ÿ�����ڵ���������
%orientation�� gabor�ĳ���0��ʾˮƽ����ֵ���ӣ���gabor��ʱ����ת����λΪ��
% phase����λ
% color1  ��ֵ��ɫ
% color2  ��ֵ��ɫ

color1=color1/255;
color2=color2/255;

[X,Y]=meshgrid(1:size_x,1:size_y);%����size_y*size_x������YΪ1:size_y�����������ظ�size_x��

a=(sin(  (sind(orientation)*X+cosd(orientation)*Y)  *  2*pi/pix_per_cycle  +phase  )+1)/2;

a1=mapminmax(a,color1(1),color2(1));%mapminmaxʹa��һ�����塢��ֵ��
a2=mapminmax(a,color1(2),color2(2));
a3=mapminmax(a,color1(3),color2(3));

gabor=cat(3,a1,a2,a3);%�ѵ�a1 a2 a3������ɫ����
gabor=im2uint8(gabor);%��gabor���黯Ϊ8bit��ͼ��
gabor(:,:,4)=128;%͸����50%��

