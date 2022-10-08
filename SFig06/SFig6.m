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

N = 100;
%%
v = gifti('HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s2_MSMAll_L_Language.func.gii');
t = v.cdata > 0;
t = t .* v.cdata;
C = v.cdata > mean(nonzeros(t));
B_Auditory = [173 107 125 175 124 174 24 104 103];
g1 = gifti('Fig4a_180areas.L.label.gii');
res_L = C;
tmp = (g1.cdata == (B_Auditory + 180)) .* C;
res_L = res_L + sum(tmp, 2) * (-40);
res_L(res_L == -39) = 0;
C = res_L == 1;
LL_Sel = LL(C, :);
LL_Sel = reshape(LL_Sel, [size(LL_Sel,1)*size(LL_Sel,2),1]);
% LanRange = linspace(mean(LL_Sel)-2*std(LL_Sel),mean(LL_Sel)+2*std(LL_Sel),N);
LanRange = linspace(0,mean(LL_Sel)+std(LL_Sel),N);
dp = [];
for i=1:size(LL,2)
    for j=1:N
        sel_vertices = (LL(:,i) > LanRange(j)) .* C;
        sel_vertices = sel_vertices == 1;
        LanDiff = LL(:,i) - LR(:,i);
        SocDiff = SL(:,i) - SR(:,i);
        dp(i,j,:) = [mean(LanDiff(sel_vertices)) mean(SocDiff(sel_vertices))];
    end
end

dpm = squeeze(nanmean(dp, 1));

m = csvread('SFig6.txt');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'Box','on');
set(gca,'FontName','arial','FontSize',10);
hold on
c = floor(linspace(1,size(m,1),N));
scatter(dpm(:,1),dpm(:,2),25,m(c,:),'filled');
box off
ylabel({'Social Task', 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
xlabel({'LH activation - RH activation', 'Language Task'},'FontName','arial','FontWeight','bold','FontSize',12);
axis square
% ylim([-7.5 -6.5])
% export_fig('Rev01_0.png', '-r600');
%%
v = gifti('HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s2_MSMAll_L_Language.func.gii');
t = v.cdata > 0;
t = t .* v.cdata;
C = v.cdata > mean(nonzeros(t));
B_Auditory = [173 107 125 175 124 174 24 104 103];
g1 = gifti('Fig4a_180areas.L.label.gii');
% LanRange = linspace(mean(v.cdata)-2*std(v.cdata),mean(v.cdata)+2*std(v.cdata),N);
tm = nonzeros(v.cdata .* (v.cdata > 0));
LanRange = linspace(0,mean(tm)+std(tm),N);
dp = [];
for i=1:size(LL,2)
    i
    for j=1:N
        sel_vertices = v.cdata > LanRange(j);
        tmp = (g1.cdata == (B_Auditory + 180)) .* C;
        sel_vertices = sel_vertices + sum(tmp, 2) * (-40);
        sel_vertices(sel_vertices == -39) = 0;
        sel_vertices = sel_vertices == 1;
        
        LanDiff = LL(:,i) - LR(:,i);
        SocDiff = SL(:,i) - SR(:,i);
        dp(i,j,:) = [mean(LanDiff(sel_vertices)) mean(SocDiff(sel_vertices))];
    end
end

dpm = squeeze(mean(dp, 1));

m = csvread('SFig6.txt');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'Box','on');
set(gca,'FontName','arial','FontSize',10);
hold on
c = floor(linspace(1,size(m,1),N));
scatter(dpm(:,1),dpm(:,2),25,m(c,:),'filled');
box off
ylabel({'Social Task', 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
xlabel({'LH activation - RH activation', 'Language Task'},'FontName','arial','FontWeight','bold','FontSize',12);
axis square;
% export_fig('Rev01_1.png', '-r600');
