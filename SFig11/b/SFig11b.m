%% Lan
clearvars
load('SFig11b.mat')
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

clear checklist1;
clear LHLan RHLan LHMath RHMath LHSvsM RHSvsM;
L = 0;
for i = 3:4:size(FileList_Language,1)
    a = FileList_Language(i).name
    if(sum(ismember(files_social,str2num(a(1:6)))))
        LH_Lan = gifti(['ExtractedData\' FileList_Language(i).name]);
        RH_Lan = gifti(['ExtractedData\' FileList_Language(i+1).name]);
        LH_Math = gifti(['ExtractedData2\' a(1:6) '_MB.L.func.gii']);
        RH_Math = gifti(['ExtractedData2\' a(1:6) '_MB.R.func.gii']);
        LH_SvsM = gifti(['ExtractedData2\' a(1:6) '_SM.L.func.gii']);
        RH_SvsM = gifti(['ExtractedData2\' a(1:6) '_SM.R.func.gii']);
        L = L + 1;
        checklist1(L) = str2num(a(1:6));
        LHLan(:,L) = LH_Lan.cdata;
        RHLan(:,L) = RH_Lan.cdata;
        LHMath(:,L) = LH_Math.cdata;
        RHMath(:,L) = RH_Math.cdata;
        LHSvsM(:,L) = LH_SvsM.cdata;
        RHSvsM(:,L) = RH_SvsM.cdata;
    end
end

LHmeanLanSocVSbase = (LHLan + LHMath) / 2;
RHmeanLanSocVSbase = (RHLan + RHMath) / 2;

LHStoryVSmean = LHLan - LHmeanLanSocVSbase;
RHStoryVSmean = RHLan - RHmeanLanSocVSbase;

LHStoryVSmean = LHStoryVSmean(:,[1:902 904:end]);
LHSvsM = LHSvsM(:,[1:902 904:end]);
RHStoryVSmean = RHStoryVSmean(:,[1:902 904:end]);
RHSvsM = RHSvsM(:,[1:902 904:end]);

LHLan = LHLan(:,[1:902 904:end]);
LHMath = LHMath(:,[1:902 904:end]);
%%
gL = gifti('../HCP_S1200_GroupAvg_v1/HCP_S1200_GroupAvg_v1/180areas.L.label.gii');
gR = gifti('../HCP_S1200_GroupAvg_v1/HCP_S1200_GroupAvg_v1/180areas.R.label.gii');
gL = gL.cdata;
gR = gR.cdata;

bars = zeros(8,2);
err = zeros(8,2);
pvals = zeros(8,1);

for i=1:8
    bars(i,1) = mean(mean(LHLan(gL==B(i)+180,:),1));
    bars(i,2) = mean(mean(LHMath(gR==B(i),:),1));
    [~,pvals(i)] = ttest(mean(LHLan(gL==B(i)+180,:),1),mean(LHMath(gR==B(i),:),1));
    err(i,1) = std( mean(LHLan(gL==B(i)+180,:),1))/sqrt(size(LHLan,2));
    err(i,2) = std(mean(LHMath(gR==B(i),:),1))/sqrt(size(LHMath,2));
end

[~,~,padj] = fdr(pvals);

f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
b = bar(1:size(A,2), bars);
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 1 0];
legend({'Story' 'Math'},'FontName','arial','FontSize',10,'Location','southeast');
legend boxoff;
set(gca,'Box','off');
ylabel('Activation relative to baseline','FontName','arial','FontWeight','bold','FontSize',12);
%xlabel('ROIs outside language-auditory network','FontName','arial','FontWeight','bold','FontSize',12);

hold on
ngroups = size(A,2);
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
e = zeros(2,1);
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    e = errorbar(x, bars(:,i), err(:,i), '.', 'HandleVisibility','off');
    e.Color = 'black';
end
set(gca,'XTickLabel', A);
% export_fig('LH_Math_Story_vs_Old_Baseline_2.png','-r600');