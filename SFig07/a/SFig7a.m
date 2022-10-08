%% Lan
clearvars
ROI_Parcellation
addpath(genpath('matlabGiftiCifti\'));

FileList_Social = dir('ExtractedData');
clear files_social;
for i = 3:4:size(FileList_Social,1)
    a = FileList_Social(i).name;
    files_social((i+1)/4) = str2num(a(1:6));
end

FileList_Language = dir('ExtractedData');
clear files_language;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    files_language((i+1)/4) = str2num(a(1:6));
end

[val,pos]=intersect(files_social, files_language);

clear checklist1;
clear DiffS;
L = 0;
for i = 3:4:size(FileList_Social,1)
    a = FileList_Social(i).name;
    if(sum(ismember(files_language,str2num(a(1:6)))))
        LH_Soc = gifti(['ExtractedData\' FileList_Language(i+2).name]);
        RH_Soc = gifti(['ExtractedData\' FileList_Language(i+3).name]);
        L = L + 1;
        checklist1(L) = str2num(a(1:6));
        DiffS(:,L) = LH_Soc.cdata - RH_Soc.cdata;
    end
end

DiffS = DiffS(:,[1:902 904:end]);

clear checklist2;
clear DiffL;
L = 0;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    if(sum(ismember(files_social,str2num(a(1:6)))))
        LH_Lan = gifti(['ExtractedData\' FileList_Language(i).name]);
        RH_Lan = gifti(['ExtractedData\' FileList_Language(i+1).name]);
        L = L + 1;
        checklist2(L) = str2num(a(1:6));
        DiffL(:,L) = LH_Lan.cdata - RH_Lan.cdata;
    end
end

DiffL = DiffL(:,[1:902 904:end]);
%%
A = {'All vertices except PGi and Auditory'};
B = [-49];

res_LB = res_L;
res_L(res_L == -39) = 0;
res_L(res_L == -29) = 0;
res_L(res_L == 1) = 0;
res_L(res_L ~= 0) = -49;

for j = 1:size(B,2)
    roi = B(j);

    ROI_Language = (res_L == roi);
    ROI_Curvature = (res_L == roi);

    ROI_DiffL = transpose(DiffL) * ROI_Language;
    ROI_DiffC = transpose(DiffS) * ROI_Curvature;

    ROI_DiffC = ROI_DiffC ./ sum(ROI_Curvature);
    ROI_DiffL = ROI_DiffL ./ sum(ROI_Language);

    t = -30;

    tail_x = ROI_DiffL(ROI_DiffL <= t);
    tail_y = ROI_DiffC(ROI_DiffL <= t);
    body_x = ROI_DiffL(ROI_DiffL > t);
    body_y = ROI_DiffC(ROI_DiffL > t);
end
%%
Sel = ROI_DiffL <= t;
DiffS_Sel = DiffS(:,Sel);
p_mapL = zeros(32492,1);
p_mapR = zeros(32492,1);
for i=1:32492
    [~, p_mapR(i)] = ttest(DiffS_Sel(i,:),0,'Tail','left'); % Greater than 0
    [~, p_mapL(i)] = ttest(DiffS_Sel(i,:),0,'Tail','right'); % Lower than 0
end
mean_DiffS_Sel = mean(DiffS_Sel,2);

res_L = res_LB;
res_L(res_L == 1) = 0;
res_L(res_L == -39) = 0;
res_L(res_L ~= 0) = 1;
sel_pL = p_mapL(res_L==1);
sel_pR = p_mapR(res_L==1);
f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
p = pie([sum(sel_pR<0.05) sum(sel_pL<0.05) sum(sel_pL>0.05 & sel_pR>0.05)],{'RH lateralization' 'LH lateralization' 'No lateralization'});
p(1).FaceColor = 'green';
p(3).FaceColor = 'red';
p(5).FaceColor = [0.5 0.5 0.5];
% export_fig('Rev06_PieChart.png','-r600');