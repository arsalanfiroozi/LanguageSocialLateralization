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
% Lan_mean = mean(DiffL,2);
% Soc_mean = mean(DiffS,2);
g1 = gifti('Fig4a_180areas.L.label.gii');

RActivity_Lan = zeros(180,1044);
RActivity_Soc = zeros(180,1044);
for i=1:180
    RActivity_Lan(i,:) = mean(DiffL(g1.cdata == (i+180),:));
    RActivity_Soc(i,:) = mean(DiffS(g1.cdata == (i+180),:));
end
Error_Lan = std(RActivity_Lan') ./ sqrt(1044);
Error_Soc = std(RActivity_Soc') ./ sqrt(1044);
Names = g1.labels.name;
Names = Names(2:181);
Lan_Diff = mean(RActivity_Lan,2);
Soc_Diff = mean(RActivity_Soc,2);

[Lan_Sorted, index_Lan] = sort(abs(mean(RActivity_Lan,2)));
Lan_Label = Names(index_Lan);
Lan_Error = Error_Lan(index_Lan);
[Soc_Sorted, index_Soc] = sort(abs(mean(RActivity_Soc,2)));
Soc_Label = Names(index_Soc);
Soc_Error = Error_Soc(index_Soc);
x=1:180;
[i,Soc_x]=knee_pt(Soc_Sorted,x); 
x=1:180;
[i,Lan_x]=knee_pt(Lan_Sorted,x); 

g2 = gifti('Fig4a_CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR_L.func.gii');
Cole_Data = g2.cdata;

g2 = gifti('Fig1e.func.gii');
C = zeros(size(g2.cdata));

selected_index_Soc = 1:180;
selected_index_Soc = selected_index_Soc(abs(mean(RActivity_Soc,2)) >= Soc_Sorted(Soc_x));
selected_index_Lan = 1:180;
selected_index_Lan = selected_index_Lan(abs(mean(RActivity_Lan,2)) >= Lan_Sorted(Lan_x));
value_Lan = Lan_Sorted .* (Lan_Sorted >= Lan_Sorted(Lan_x));

for i=1:length(selected_index_Soc)
    C(g1.cdata == (selected_index_Soc(i) + 180)) = Soc_Diff(selected_index_Soc(i));
end

for i=1:length(selected_index_Lan)
    C(g1.cdata == (selected_index_Lan(i) + 180)) = Lan_Diff(selected_index_Lan(i));
end

g3 = gifti('Fig1e.func.gii');
g3.cdata = C;
% save(g3,'All_ROIs.func.gii');

C = zeros(size(g2.cdata));
Filtered_Regions = [3 6 8];
Filtered_ROIs = [76 150 130 176];
T = C;

for i=1:length(selected_index_Lan)
    if((sum(selected_index_Lan(i) == Filtered_ROIs) == 0) && (sum(mode(Cole_Data(g1.cdata == (selected_index_Lan(i) + 180))) == Filtered_Regions) == 0))
        C(g1.cdata == (selected_index_Lan(i) + 180)) = i;
        T(g1.cdata == (selected_index_Lan(i) + 180)) = selected_index_Lan(i);
    end
end

g3 = gifti('Fig1e.func.gii');
g3.cdata = C;
% save(g3,'Selected_ROIs.func.gii');

res_L = T;
A = Names(nonzeros(unique(T)));
A = erase(A,'_ROI');
A = erase(A,'R_');
B = nonzeros(unique(T));
values = B;

vals = 1:8;
for i=1:8
    vals(i) = Lan_Diff(B(i));
end
[~,I] = sort(vals);
B = B(I);
A = A(I);
values = B;

B = B(8:-1:1);
A = A(8:-1:1);
values = B;
%% Bar
bar_data = zeros(size(A,2), 2);
bar_data_LH = zeros(size(A,2), 2);
bar_data_RH = zeros(size(A,2), 2);
bar_error = zeros(size(A,2), 2);
bar_error_LH = zeros(size(A,2), 2);
bar_error_RH = zeros(size(A,2), 2);

pvalsLanSoc = zeros(8,1);
pvalsLan = zeros(8,1);
pvalsSoc = zeros(8,1);

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
    [~,pvalsLan(i)] = ttest(Lan_Data);
    [~,pvalsSoc(i)] = ttest(Soc_Data);
    
    Lan_vertices = LanLH(logical(selected_vertices), :);
    Soc_vertices = SocLH(logical(selected_vertices), :);
    Lan_bar = mean(Lan_vertices(:));
    Soc_bar = mean(Soc_vertices(:));
    bar_data_LH(i, :) = [Lan_bar Soc_bar];
    
    Lan_Data = mean(Lan_vertices, 1);
    Soc_Data = mean(Soc_vertices, 1);
    std_Lan = std(Lan_Data);
    std_Soc = std(Soc_Data);
    err_Lan = std_Lan / sqrt(length(Lan_Data));
    err_Soc = std_Soc / sqrt(length(Soc_Data));
    bar_error_LH(i, :) = [err_Lan err_Soc];
    
    Lan_vertices = LanRH(logical(selected_vertices), :);
    Soc_vertices = SocRH(logical(selected_vertices), :);
    Lan_bar = mean(Lan_vertices(:));
    Soc_bar = mean(Soc_vertices(:));
    bar_data_RH(i, :) = [Lan_bar Soc_bar];
    
    Lan_Data = mean(Lan_vertices, 1);
    Soc_Data = mean(Soc_vertices, 1);
    
    std_Lan = std(Lan_Data);
    std_Soc = std(Soc_Data);
    err_Lan = std_Lan / sqrt(length(Lan_Data));
    err_Soc = std_Soc / sqrt(length(Soc_Data));
    bar_error_RH(i, :) = [err_Lan err_Soc];
end

[~,~,pvalsLan] = fdr(pvalsLan);
[~,~,pvalsSoc] = fdr(pvalsSoc);

f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
b = bar(1:size(A,2), bar_data);
b(1).FaceColor = [1 1 0];
b(2).FaceColor = [0 1 1];
legend({'Language Task' 'Social Task'},'FontName','arial','FontSize',10);
legend boxoff;
set(gca,'Box','off');
ylabel('LH activation - RH activation','FontName','arial','FontWeight','bold','FontSize',12);
xlabel('ROIs outside language-auditory network','FontName','arial','FontWeight','bold','FontSize',12);

hold on
ngroups = size(A,2);
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
% export_fig('fig2_barplot.png','-r600'); % dpi: 600