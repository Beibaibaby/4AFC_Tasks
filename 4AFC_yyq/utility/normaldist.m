function y=normaldist(mu,sigma,x)
% 一维高斯分布，x为用于横轴取值的向量，在该主程序中未引用
% mu=move_time(thistrial(2))/2+fade_in;
% sigma=move_time(thistrial(2))/2;
% dt=move_time(thistrial(2))/(exp_hz*thismove_time);
% x=0:dt:(2*mu);
pd=makedist('Normal','mu',mu,'sigma',sigma);
y=pdf(pd,x);
%yt=yt/yt(round(length(xt)/2));
