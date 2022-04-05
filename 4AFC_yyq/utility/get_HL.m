function [H,L]=get_HL(this_contrast)
% 依据Michelson对比度公式， 对比度为(最亮-最暗)/(最亮+最暗)
%对比度为0时，亮度为128
%H L分别为最亮、最暗的值
L=(256-this_contrast*256)/2;
H=256-L;
L=round(L-1);
H=round(H-1);