function [H,L]=get_HL(this_contrast)
% ����Michelson�Աȶȹ�ʽ�� �Աȶ�Ϊ(����-�)/(����+�)
%�Աȶ�Ϊ0ʱ������Ϊ128
%H L�ֱ�Ϊ���������ֵ
L=(256-this_contrast*256)/2;
H=256-L;
L=round(L-1);
H=round(H-1);