load SA_block1_5
data = data_all;
load SA_block2_5
data = [data; data_all];
load SA_block3_5
data = [data; data_all];
load SA_block4_5
data = [data; data_all];
load SA_block5_5
data = [data; data_all];
load SA_block6_5
data = [data; data_all];
% load KP_block7_5
% data = [data; data_all];
% load KP_block8_5
% data = [data; data_all];
% load KP_block9_5
% data = [data; data_all];
% load KP_block10_5
% data = [data; data_all];
% load OK_block11_5
% data = [data; data_all];
% load OK_block12_5
% data = [data; data_all];

doi = data( data(:,2) == 0.03, : );
%doi = data( 1:480, : );

diff_all = doi(:,3) - doi(:,4);

diff_all = diff_all + 90;

diff_all( diff_all < - 180 ) = 360 + diff_all( diff_all < - 180 );
diff_all( diff_all >= 180 ) = diff_all( diff_all >= 180 ) - 360;
