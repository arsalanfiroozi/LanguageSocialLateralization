%% Asymmetry with behavior
clearvars
ROI_Parcellation

load('SFig9_TimeAccuracy.mat');
load('SFig9_subj_lan_measures.mat');
load('SFig9_subj_soc_measures.mat');
load('SFig9_id_subj.mat');
addpath(genpath('matlabGiftiCifti\'));

FileList_Language = dir('ExtractedData');
clear files_language;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    files_language((i+1)/4) = str2num(a(1:6));
end

clear checklist3;
clear DiffS;
clear C;
L = 0;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    LH_Soc = gifti(['ExtractedData\' FileList_Language(i+2).name]);
    RH_Soc = gifti(['ExtractedData\' FileList_Language(i+3).name]);
    L = L + 1;
    checklist3(L) = str2num(a(1:6));
    DiffS(:,L) = LH_Soc.cdata - RH_Soc.cdata;
    C(L,:) = T(ismember(T(:,1),str2num(a(1:6))),:);
    fprintf('%d ',i);
end

clear checklist2;
clear DiffL;
clear Subj_Lan;
clear Subj_Soc;
L = 0;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name;
    LH_Lan = gifti(['ExtractedData\' FileList_Language(i).name]);
    RH_Lan = gifti(['ExtractedData\' FileList_Language(i+1).name]);
    L = L + 1;
    checklist2(L) = str2num(a(1:6));
    DiffL(:,L) = LH_Lan.cdata - RH_Lan.cdata;
    Subj_Lan(L,:) = subj_language_measures(ismember(id,str2num(a(1:6))),:);
    Subj_Soc(L,:) = subj_social_measures(ismember(id,str2num(a(1:6))),:);
    fprintf('%d ',i);
end

DiffL = DiffL(:,[1:902 904:end]);
DiffS = DiffS(:,[1:902 904:end]);
Subj_Lan = Subj_Lan([1:902 904:end],:);
Subj_Soc = Subj_Soc([1:902 904:end],:);
% DiffS = DiffS(:, (Subj_Soc(:,1) ~= 0));
Subj_Soc = Subj_Soc((Subj_Soc(:,1) ~= 0), :);
DiffL = DiffL(:, (Subj_Lan(:,1) ~= 0));
Subj_Lan = Subj_Lan((Subj_Lan(:,1) ~= 0), :);
C = C([1:902 904:end], :);
%% Lan
A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
B = [-9 41 11 -29 -19 -49 31 21];

Y1 = Subj_Lan(:,2);
Y2 = Subj_Lan(:,4);
Y1 = Y1 - mean(Y1);
Y2 = Y2 - mean(Y2);
Y1 = Y1 ./ std(Y1);
Y2 = Y2 ./ std(Y2);
Y = (Y1 + Y2) ./ 2;
% Y = Y2;

G1 = Y <= (mean(Y) - std(Y));
G2 = (Y > (mean(Y) - std(Y))) & (Y < (mean(Y) + std(Y)));
G3 = Y >= (mean(Y) + std(Y));
data_anova = [];
g1_anova = [];
g2_anova = [];
Bar = zeros(8,3);
Errors = zeros(8,3);
Anova_y = [];
Anova_g1 = {};
Anova_P = zeros(8,3);
ttest_p2 = zeros(8,1);
for i=1:size(B,2)
    ROI = res_L == B(i);
    ROI_DiffL = (DiffL' * ROI) ./ sum(ROI);
    Errors(i,1) = std(ROI_DiffL(G1)) / sqrt(length(ROI_DiffL(G1)));
    Errors(i,2) = std(ROI_DiffL(G2)) / sqrt(length(ROI_DiffL(G2)));
    Errors(i,3) = std(ROI_DiffL(G3)) / sqrt(length(ROI_DiffL(G3)));
    Bar(i,1) = mean(ROI_DiffL(G1));
    data_anova = [data_anova;ROI_DiffL(G1)];
    g1_anova = [g1_anova;ones(size(ROI_DiffL(G1)))*1];
    Bar(i,2) = mean(ROI_DiffL(G2));
    Bar(i,3) = mean(ROI_DiffL(G3));
    data_anova = [data_anova;ROI_DiffL(G3)];
    g1_anova = [g1_anova;ones(size(ROI_DiffL(G3)))*2];

    g2_anova = [g2_anova;ones(size(ROI_DiffL(G3 | G1)))*i];

    [~,ttest_p2(i),~,~] = ttest2(ROI_DiffL(G1),ROI_DiffL(G3),'Tail','left');
end

[p,tbl,stats] = anovan(data_anova,{g1_anova, g2_anova},'model','interaction','varnames',{'g1','g2'});
[results,m,h,gnames] = multcompare(stats,'Dimension',[1 2]);    

addpath(genpath('matlabGiftiCifti\'));
figure;
hold on;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
% b = bar(Bar(:,end:-1:1));
b = bar(Bar(:,[3 1]));
tmp = autumn;
tmp = tmp(175:-1:1,:);
b(1).FaceColor = tmp(end,:);
% b(2).FaceColor = tmp(end/2-7,:);
b(2).FaceColor = tmp(1,:);
set(gca,'XTickLabel', A);
% h = legend({'High' 'Medium' 'Low'},'Location','northwest');
h = legend({'High' 'Low'},'Location','northwest');
legend boxoff;
title(h, {'Average score in language tests'});
set(gca,'Box','off');
add_errorbar(Errors(:,[3 1]), Bar(:,[3 1]));
xlabel('Language areas','FontName','arial','FontWeight','bold','FontSize',12);
ylabel({'Language Task' 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
sigasterisk(2,1,4,4,'*',Bar(:,[3 1]),Errors(:,[3 1]));
sigasterisk(2,1,5,5,'*',Bar(:,[3 1]),Errors(:,[3 1]));
ylim([0 80]);
% export_fig(['fig6_barplot/Language_ZScore.png'],'-r600');

%% Soc
Names = {'Social_Task_Perc_Random'	'Social_Task_Perc_TOM'	'Social_Task_Perc_Unsure'	'Social_Task_Perc_NLR'	'Social_Task_Median_RT_Random'	'Social_Task_Median_RT_TOM'	'Social_Task_Median_RT_Unsure'	'Social_Task_Random_Perc_Random'	'Social_Task_Random_Median_RT_Random'	'Social_Task_Random_Perc_TOM'	'Social_Task_Random_Median_RT_TOM'	'Social_Task_Random_Perc_Unsure'	'Social_Task_Random_Median_RT_Unsure'	'Social_Task_Random_Perc_NLR'	'Social_Task_TOM_Perc_Random'	'Social_Task_TOM_Median_RT_Random'	'Social_Task_TOM_Perc_TOM'	'Social_Task_TOM_Median_RT_TOM'	'Social_Task_TOM_Perc_Unsure'	'Social_Task_TOM_Median_RT_Unsure'	'Social_Task_TOM_Perc_NLR'};
A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
B = [-9 41 11 -29 -19 -49 31 21];

Bars = zeros(8, 3);
Errors = zeros(8, 3);
Anova_y = [];
Anova_g1 = {};
Anova_P = zeros(8,3);
ttest_p = zeros(8,1);
for i=1:length(B)
    f = (C(:,18) == -1) | (C(:,9) == -1);
    ROI = res_L == B(i);
    ROI_DiffS = (DiffS' * ROI) ./ sum(ROI);
    ROI_DiffS = ROI_DiffS(~f);
    Y = (C(:,18) + C(:,9)) ./ 2;
    Y = Y(~f);
    G1 = Y <= (median(Y) - mad(Y,1));
    G2 = (Y > (median(Y) - mad(Y,1))) & (Y < (median(Y) + mad(Y,1)));
    G3 = Y >= (median(Y) + mad(Y,1));
    Bars(i,1) = mean(ROI_DiffS(G1));
    Bars(i,2) = mean(ROI_DiffS(G2));
    Bars(i,3) = mean(ROI_DiffS(G3));
    Errors(i,1) = std(ROI_DiffS(G1)) / sqrt(sum(G1));
    Errors(i,2) = std(ROI_DiffS(G2)) / sqrt(sum(G2));
    Errors(i,3) = std(ROI_DiffS(G3)) / sqrt(sum(G3));
    
    [~,ttest_p(i),~,~] = ttest2(ROI_DiffS(G1),ROI_DiffS(G3),'Tail','right');
%     ttest_p(i) = ranksum(ROI_DiffS(G1),ROI_DiffS(G3),'Tail','right');
    
    Gs = {'Low' 'Medium' 'High'};
    Anova_y = [ROI_DiffS(G1)' ROI_DiffS(G2)' ROI_DiffS(G3)'];
    Anova_g1 = {Gs{ones(1, sum(G1))} Gs{ones(1, sum(G2))*2} Gs{ones(1, sum(G3))*3}};
    [p,tbl,stats] = anovan(Anova_y,{Anova_g1});
    [results,m,h,gnames] = multcompare(stats);
    Anova_P(i,:) = results(:,6)';
end

[~,~,ttest_p_corr] = fdr(ttest_p);

figure;
hold on;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
% b = bar(Bars(:,end:-1:1));
b = bar(Bars(:,[3 1]));
tmp = winter;
b(1).FaceColor = tmp(end,:);
b(2).FaceColor = tmp(1,:);
set(gca,'XTickLabel', A);
% h = legend({'High' 'Medium' 'Low'},'Location','southeast');
h = legend({'High' 'Low'},'Location','southeast');
legend boxoff;
title(h, 'Accuracy in social task');
set(gca,'Box','off');
add_errorbar(Errors(:,[3 1]), Bars(:,[3 1]));
xlabel('Language areas','FontName','arial','FontWeight','bold','FontSize',12);
ylabel({'Social Task' 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
sigasterisk(2,1,1,1,'*',Bars(:,[3 1]),Errors(:,[3 1]),35 * 2 / 125);
sigasterisk(2,1,2,2,'*',Bars(:,[3 1]),Errors(:,[3 1]),35 * 2 / 125);

height_v = ylim;
height = height_v(2) - height_v(1);
threshold1 = height / 125;
threshold2 = height * 2 / 125;
threshold3 = height * 1.5 / 125;

ylim([-18 4]);
% export_fig('fig6_barplot/Acc.png','-r600');
