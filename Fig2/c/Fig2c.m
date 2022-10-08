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

res_L(res_L == -39) = 0;
res_L(res_L == -29) = 0;
res_L(res_L == 1) = 0;
res_L(res_L ~= 0) = -49;

for j = 1:size(B,2)
    roi = B(j);

    ROI_Language = (res_L == roi);
    ROI_Social = (res_L == roi);

    ROI_DiffL = transpose(DiffL) * ROI_Language;
    ROI_DiffS = transpose(DiffS) * ROI_Social;

    ROI_DiffS = ROI_DiffS ./ sum(ROI_Social);
    ROI_DiffL = ROI_DiffL ./ sum(ROI_Language);

    t = -30;

    tail_x = ROI_DiffL(ROI_DiffL <= t);
    tail_y = ROI_DiffS(ROI_DiffL <= t);
    body_x = ROI_DiffL(ROI_DiffL > t);
    body_y = ROI_DiffS(ROI_DiffL > t);

    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    scatter(body_x, body_y, 10,'o', 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    hold on;
    scatter(tail_x, tail_y, 10,'o', 'filled', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
    plot( t * ones(301,1) , [-150:150], [':', 'k']);

    axis square;
    ylabel({'Social Task' 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
    xlabel({'LH activation - RH activation' 'Language Task'},'FontName','arial','FontWeight','bold','FontSize',12);
    %export_fig('fig2_scatterplot_less_than_-30.png','-r600');
end

q = ROI_DiffL <= t;
%     q(903) = 0;
selected = DiffS(:,q);
vertices = mean(selected, 2);
g1 = gifti('Fig1e.func.gii');
g1.cdata = vertices;
save(g1, 'fig2_selected_subjects_Social_less_-30.func.gii');

q = ROI_DiffL <= t;
%     q(903) = 0;
selected = DiffL(:,q);
vertices = mean(selected, 2);
g1 = gifti('Fig1e.func.gii');
g1.cdata = vertices;
save(g1, 'fig2_selected_subjects_Lan_less_-30.func.gii');