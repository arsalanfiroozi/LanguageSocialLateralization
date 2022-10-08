clearvars
addpath(genpath('matlabGiftiCifti\'));
FileList = dir('ExtractedData');
%% Language
for i = 1:1045
    LH = gifti(['ExtractedData\' FileList((i*4)-1).name]);
    RH = gifti(['ExtractedData\' FileList((i*4)+0).name]);
    DiffL(:,i) = LH.cdata - RH.cdata;
    LanLH(:,i) = LH.cdata;
    LanRH(:,i) = RH.cdata;
    fprintf('%d ',i);
end

DiffL = DiffL(:,[1:902 904:end]);
LanLH = LanLH(:,[1:902 904:end]);
LanRH = LanRH(:,[1:902 904:end]);

%% Social
for i = 1:1045
    LH = gifti(['ExtractedData\' FileList((i*4)+1).name]);
    RH = gifti(['ExtractedData\' FileList((i*4)+2).name]);
    DiffS(:,i) = LH.cdata - RH.cdata;
    SocLH(:,i) = LH.cdata;
    SocRH(:,i) = RH.cdata;
    fprintf('%d ',i);
end

DiffS = DiffS(:,[1:902 904:end]);
SocLH = SocLH(:,[1:902 904:end]);
SocRH = SocRH(:,[1:902 904:end]);
%% Lan
map = zeros(32492,1);
for i=1:32492
    [h,ttest_p] = ttest(DiffL(i,:),0);
    [~,~,map(i)] = fdr(ttest_p);
end
map = -log10(map);
map = map .* sign(mean(DiffL,2));

g = gifti('Fig1f.func.gii');
g.cdata = map;
save(g, 'SFig4_Lan.func.gii');
%% Soc
map = zeros(32492,1);
for i=1:32492
    [h,ttest_p] = ttest(DiffS(i,:),0);
    [~,~,map(i)] = fdr(ttest_p);
end
map = -log10(map);
map = map .* sign(mean(DiffS,2));

g = gifti('Fig1f.func.gii');
g.cdata = map;
save(g, 'SFig4_Soc.func.gii');