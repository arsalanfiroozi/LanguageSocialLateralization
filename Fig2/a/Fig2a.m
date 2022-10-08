clearvars
ROI_Parcellation

t = -30;
load('ROI_DiffL.mat');

addpath(genpath('matlabGiftiCifti\'));
FileList = dir('ExtractedData');

% A = {'Broca' 'SFL' '55b' 'STSa' 'STSp' 'PGi' 'PSL' 'STGa'};
% values = [-49 21 31 41 -9 -29 -19 11];

A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
values = [-9 41 11 -29 -19 -49 31 21];

%% Language
for i = 1:1045
    LH = gifti(['ExtractedData\' FileList((i*4)-1).name]);
    RH = gifti(['ExtractedData\' FileList((i*4)+0).name]);
    DiffL(:,i) = LH.cdata - RH.cdata;
    fprintf('%d ',i);
end

DiffL = DiffL(:,[1:902 904:end]);
% DiffL = DiffL(:,ROI_DiffL <= t);

%% Social
for i = 1:1045
    LH = gifti(['ExtractedData\' FileList((i*4)+1).name]);
    RH = gifti(['ExtractedData\' FileList((i*4)+2).name]);
    DiffS(:,i) = LH.cdata - RH.cdata;
    fprintf('%d ',i);
end

DiffS = DiffS(:,[1:902 904:end]);
% DiffS = DiffS(:,ROI_DiffL <= t);

%% Bar

bar_data = zeros(8, 2);
bar_error = zeros(8, 2);
p_data = zeros(8, 2);

for i = 1:size(A, 2)
    ROI = values(i);
    selected_vertices = res_L == ROI;
    
    Lan_vertices = DiffL(logical(selected_vertices), :);
    Soc_vertices = DiffS(logical(selected_vertices), :);
    Lan_bar = mean(Lan_vertices(:));
    Soc_bar = mean(Soc_vertices(:));
    bar_data(i, :) = [Lan_bar Soc_bar];
    
    Lan_Data = mean(Lan_vertices, 1);
    Soc_Data = mean(Soc_vertices, 1);
    std_Lan = std(Lan_Data);
    std_Soc = std(Soc_Data);
    err_Lan = std_Lan / sqrt(length(Lan_Data));
    err_Soc = std_Soc / sqrt(length(Soc_Data));
    bar_error(i, :) = [err_Lan err_Soc];
    
    [~, p] = ttest(Lan_Data);
    p_data(i,1) = p;
    [~, p] = ttest(Soc_Data);
    p_data(i,2) = p;
end

p_data_fdr = zeros(8,2);
[~,~,p_data_fdr(:,1)] = fdr(p_data(:,1));
[~,~,p_data_fdr(:,2)] = fdr(p_data(:,2));

f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
b = bar(1:8, bar_data);
b(1).FaceColor = [1 1 0];
b(2).FaceColor = [0 1 1];
% legend({'Language Task' 'Social Task'},'Location','southeast','FontName','arial','FontSize',10);
legend({'Language Task' 'Social Task'},'FontName','arial','FontSize',10);
legend boxoff;
set(gca,'Box','off');
ylabel('LH activation - RH activation','FontName','arial','FontWeight','bold','FontSize',12);
xlabel('Language areas','FontName','arial','FontWeight','bold','FontSize',12);

hold on
ngroups = 8;
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
e = zeros(2,1);
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    e = errorbar(x, bar_data(:,i), bar_error(:,i), '.', 'HandleVisibility','off');
    e.Color = 'black';
end
set(gca,'XTickLabel', A);

% export_fig('sfig2_barplot.png','-r600'); % dpi: 600
