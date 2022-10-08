clearvars;
addpath(genpath('matlabGiftiCifti\'));
FileList = dir('ExtractedData');

for i = 1:1045
    a = FileList((i*4)-1).name;
    LH = gifti(['ExtractedData\' FileList((i*4)-1).name]);
    RH = gifti(['ExtractedData\' FileList((i*4)+0).name]);
    LanLH(:,i) = LH.cdata;
    LanRH(:,i) = RH.cdata;
    checklist1(i) = str2num(a(1:6));
    fprintf('%d ',i);
end

LanLH = LanLH(:,[1:902 904:end]);
LanRH = LanRH(:,[1:902 904:end]);
checklist1 = checklist1([1:902 904:end]);

for i = 1:1045
    a = FileList((i*4)+1).name;
    LH = gifti(['ExtractedData\' FileList((i*4)+1).name]);
    RH = gifti(['ExtractedData\' FileList((i*4)+2).name]);
    SocLH(:,i) = LH.cdata;
    SocRH(:,i) = RH.cdata;
    checklist2(i) = str2num(a(1:6));
    fprintf('%d ',i);
end

SocLH = SocLH(:,[1:902 904:end]);
SocRH = SocRH(:,[1:902 904:end]);
checklist2 = checklist2([1:902 904:end]);

Id = checklist1;

%%
addpath(genpath('matlabGiftiCifti\'));
FileList = dir('ExtractedData');

LanLH = mean(LanLH,2);
LanRH = mean(LanRH,2);
SocLH = mean(SocLH,2);
SocRH = mean(SocRH,2);
LH = gifti(['ExtractedData\' FileList((1*4)-1).name]);
RH = gifti(['ExtractedData\' FileList((1*4)+0).name]);

LH.cdata = LanLH;
RH.cdata = LanRH;
save(LH,'1044_LH_LanMean.func.gii');
save(RH,'1044_RH_LanMean.func.gii');

LH.cdata = SocLH;
RH.cdata = SocRH;
save(LH,'1044_LH_SocMean.func.gii');
save(RH,'1044_RH_SocMean.func.gii');

LH.cdata = SocLH - SocRH;
save(LH,'1044_Diff_SocMean.func.gii');
LH.cdata = LanLH - LanRH;
save(LH,'1044_Diff_LanMean.func.gii');
