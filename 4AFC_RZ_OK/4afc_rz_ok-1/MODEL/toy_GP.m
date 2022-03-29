clear
clc


x = 1:0.1:5;
x = x';
y = 2*randn( size(x) ) + 5*sin(4*x) - 1;

z = linspace( 0, 6, 100 );
z = z';



% hyp.mean = 1;
% hyp.cov = [0;1];
% hyp.lik = log(2);
% 
% hyp2 = minimize(hyp, @gp, -1000, inf, meanfunc, covfunc, likfunc, x, y);
% [m, s2] = gp( hyp2, inf, meanfunc, covfunc, likfunc, x, y, z);



[m, s2] = gaussian_process( x, y, z,...
    2.5, 1.3, 0.46, 1 );
hold off
% f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)]; 
% hold on
% fill([z; flipdim(z,1)], f, [7 7 7]/8)
plot( z, m )
hold on
plot( x, y, '*' )
plot( z, m+diag(sqrt(s2)) )
plot( z, m-diag(sqrt(s2)) )

hold off

