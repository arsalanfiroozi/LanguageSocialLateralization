%% 
clearvars
addpath(genpath('matlabGiftiCifti\'));

FileList_Language = dir('ExtractedData');
clear files_language;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    files_language((i+1)/4) = str2num(a(1:6));
end

clear checklist3;
clear SL SR;
L = 0;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    LH_Soc = gifti(['ExtractedData\' FileList_Language(i+2).name]);
    RH_Soc = gifti(['ExtractedData\' FileList_Language(i+3).name]);
    L = L + 1;
    checklist3(L) = str2num(a(1:6));
    SL(:,L) = LH_Soc.cdata;
    SR(:,L) = RH_Soc.cdata;
    fprintf('%d ',i);
end

clear checklist2;
clear LL LR;
L = 0;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    LH_Lan = gifti(['ExtractedData\' FileList_Language(i).name]);
    RH_Lan = gifti(['ExtractedData\' FileList_Language(i+1).name]);
    L = L + 1;
    checklist2(L) = str2num(a(1:6));
    LL(:,L) = LH_Lan.cdata;
    LR(:,L) = RH_Lan.cdata;
    fprintf('%d ',i);
end

LL = LL(:,[1:902 904:end]);
LR = LR(:,[1:902 904:end]);
SL = SL(:,[1:902 904:end]);
SR = SR(:,[1:902 904:end]);
%%
LD = LL - LR;
SD = SL - SR;
LD = mean(LD, 2);
SD = mean(SD, 2);
LDB = LD > mean(LD(LD>0));
LDS = LD < mean(LD(LD<0));
SDB = SD > mean(SD(SD>0));
SDS = SD < mean(SD(SD<0));
t = LDB .* SDB; % R
t = t + 2 * LDS .* SDB; % B
t = t + 3 * LDB .* SDS; % G
t = t + 4 * LDS .* SDS; % Y
LH_Lan.cdata = t;
save(LH_Lan, 'SFig5.func.gii');
% .\wb_command.exe -metric-label-import Rev03.func.gii SFig5.txt Rev03.label.gii