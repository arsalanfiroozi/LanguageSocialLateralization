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
%%
res_L(res_L ~= 0) = 1;
A = {'All ROIs outside language-auditory network'};
B = 1;
values= 1;
for i=1:size(A,2)
    selected_vertices = res_L == values(i);
    DiffL_mean = transpose(DiffL) * selected_vertices;
    DiffS_mean = transpose(DiffS) * selected_vertices;
    DiffL_mean = DiffL_mean / sum(selected_vertices);
    DiffS_mean = DiffS_mean / sum(selected_vertices);
    figure;
    hold on;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    scatter(DiffL_mean, DiffS_mean, 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0 1]);
    axis square;
    XLim = xlim;
    XLim = ceil(max([abs(XLim(1)) abs(XLim(2))]));
    xlim([-XLim XLim]);
    xticks(-XLim:XLim/2:XLim);
    YLim = ylim;
    YLim = max([abs(YLim(1)) abs(YLim(2))]);
    ylim([-YLim YLim]);
    yticks(-YLim:YLim/2:YLim);
    h = lsline;
    plot(zeros(2 * YLim + 1,1), [-YLim:YLim], [':', 'k']);
    plot([-XLim:XLim], zeros(2 * XLim + 1,1), [':', 'k']);
    h.Color = 'r';
    h.LineWidth = 1;
    [r,p] = corrcoef(DiffL_mean, DiffS_mean);
    if(mod(floor(abs(r(1,2)*100)),100) == 0)
        R = sprintf('%1.2E', r(1,2));
    else
        R = sprintf('%1.2f', r(1,2));
    end
    if(p(1,2) > 0.05)
        P = '> 0.05';
    else
        P = '< 0.05';
    end
    if(p(1,2) < 0.05)
        P = [P ' *'];
        if(p(1,2) < 0.005)
            P = [P '*'];
            if(p(1,2) < 0.0005)
                P = [P '*'];
            end
        end
    end
    title({['\bf ' A{i}] ['\rm r = ' R ' , p ' P]},'FontName','arial','FontSize',12);
    set(gca,'Box','off');
    ylabel({'Social Task' 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
    xlabel({'LH activation - RH activation' 'Language Task'},'FontName','arial','FontWeight','bold','FontSize',12);
end

