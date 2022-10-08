clearvars;
addpath(genpath('matlabGiftiCifti\'));
FileList = dir('ExtractedData');

% positive cohenced
v = gifti('HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s2_MSMAll_L_Language.func.gii');
t = v.cdata > 0;
t = t .* v.cdata;
C = v.cdata > mean(nonzeros(t));

g1 = gifti('../HCP_S1200_GroupAvg_v1/HCP_S1200_GroupAvg_v1/180areas.L.label.gii');
g2 = gifti('../HCP_S1200_GroupAvg_v1/HCP_S1200_GroupAvg_v1/180areas.R.label.gii');

% PSL: 25
% 31pd: 161
% STSVA: 176
% STSVP: 130
% STSDP: 129
% STSDA: 128
% 55b: 12
% SFL: 26
% 44: 74
% 45: 75
% PGi: 150
% STGa: 123
% 47L: 76

A = [ ['47L' 45 'IFJa' 'IFSp' 44] ['SFL' 'SCEF' '8BM'] ['55b' '8AV' 4 'FEF'] ['STSda' 'STSva' 'TE1a'] ['STSvp' 'STSdp' 'TPOJ1' 'PHT'] ['PGi' 'STV'] ['PSL'] ['MBelt' 'TA2' 'A5' 'A4' 'PBelt' 'LBelt' 'A1' 'RI' 52] ['STGa' 'TGd']];
B_Broca = [76 75 79 81 74];
B_SFL = [26 43 63];
B_55b = [12 67 8 10];
B_STSa = [128 176 132];
B_STSp = [130 129 139 137];
B_PGi = [150 28];
B_PSL = [25];
B_Auditory = [173 107 125 175 124 174 24 104 103];
B_STGa = [123 131];

A = {'Broca' 'SFL' '55b' 'STSa' 'STSp' 'PGi' 'PSL' 'Auditory' 'STGa'};
values = [-49 21 31 41 -9 -29 -19 -39 11];

res_L = C;

% Broca
tmp = (g1.cdata == (B_Broca + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (-50);

% SFL
tmp = (g1.cdata == (B_SFL + 180)) .* C;
res_L = res_L + sum(tmp, 2) * 20;

% 55b
tmp = (g1.cdata == (B_55b + 180)) .* C;
res_L = res_L + sum(tmp, 2) * 30;

% STSa
tmp = (g1.cdata == (B_STSa + 180)) .* C;
res_L = res_L + sum(tmp, 2) * 40;

% STSp
tmp = (g1.cdata == (B_STSp + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (-10);

% PGi
tmp = (g1.cdata == (B_PGi + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (-30);

% PSL
tmp = (g1.cdata == (B_PSL + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (-20);

% Auditory
tmp = (g1.cdata == (B_Auditory + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (-40);

% STGa
tmp = (g1.cdata == (B_STGa + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (10);

g3 = gifti('HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s2_MSMAll_L_Language.func.gii');
s = g3;

tmp = res_L;
tmp(res_L == 1) = 0;
tmp(res_L == -49) = 11;
tmp(res_L == 21) = 8;
tmp(res_L == 31) = 9;
tmp(res_L == 41) = 16;
tmp(res_L == -9) = 2;
tmp(res_L == -29) = 5;
tmp(res_L == -19) = 6;
tmp(res_L == -39) = 0;
tmp(res_L == 11) = 7;
tmp(res_L == 0) = 0;

s.cdata = tmp;
% save(s, 'ROI_Parcellation_L.func.gii');

g3 = gifti('ExtractedData/990366_Lan.R.func.gii');
s = g3;
s.cdata = tmp;
% save(s, 'ROI_Parcellation_R.func.gii')

s = g3;
s.cdata = (res_L == 1);
% save(s, 'ROI_Parcellation_missed_vertices.func.gii');

s = gifti('ExtractedData/990366_Lan.L.func.gii');
s.cdata = C;
save(s, 'SFig2.func.gii');
