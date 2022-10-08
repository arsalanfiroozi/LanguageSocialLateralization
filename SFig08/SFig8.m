%% Lan
clearvars;
ROI_Parcellation

addpath(genpath('matlabGiftiCifti\'));
load('SFig8_Handedness.mat');
load('SFig8_id.mat')

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

clear checklist0;
clear handedness;
L = 0;
for i = 1:length(id)
    if(sum(ismember(files_language,id(i,1))))
        L = L + 1;
        checklist0(L) = id(i,1);
        handedness(L) = Data_Handedness(i);
    end
end
handedness = handedness([1:902 904:end]);

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
DiffLB = DiffL;
DiffSB = DiffS;
handednessB = handedness;

handedness = handedness > 0;
DiffL = DiffL(:,handedness);
DiffS = DiffS(:,handedness);

if(~sum(checklist1 - checklist2))
    A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
    values = [-9 41 11 -29 -19 -49 31 21];
    
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
        xlim([-200 200]);
        ylim([-100 100]);
        h = lsline;
        h.Color = 'r';
        h.LineWidth = 1;
        [r,p] = corrcoef(DiffL_mean, DiffS_mean);
        if(mod(floor(abs(r(1,2)*100)),100) == 0)
            R = sprintf('%1.2E', r(1,2));
        else
            R = sprintf('%1.2f', r(1,2));
        end
        if(i == 6)
            R = '0.0016';
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
        %title({['\fontsize{14}', AirfoilName];'\fontsize{10}\color{blue} ';'Coefficient of Lift';'vs Angle of Attack'},'fontweight','b');
        set(gca,'Box','off');
        ylabel({'Social Task' 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
        xlabel({'LH activation - RH activation' 'Language Task'},'FontName','arial','FontWeight','bold','FontSize',12);
        plot(zeros(201,1), [-100:100], [':', 'k']);
        plot([-200:200], zeros(401,1), [':', 'k']);
        xlim([-200 200]);
        ylim([-100 100]);
        yticks([-100 -50 0 50 100]);
        yticklabels({'-100','-50','0','50','100'});
        axis square;
%         export_fig(['fig3_scatter_right_handed/' A{i} '.png'], '-r600');
    end
end

%%

handedness = handednessB < 0;
DiffL = DiffLB(:,handedness);
DiffS = DiffSB(:,handedness);

if(~sum(checklist1 - checklist2))
    A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
    values = [-9 41 11 -29 -19 -49 31 21];
    
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
        xlim([-200 200]);
        ylim([-100 100]);
        h = lsline;
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
        %title({['\fontsize{14}', AirfoilName];'\fontsize{10}\color{blue} ';'Coefficient of Lift';'vs Angle of Attack'},'fontweight','b');
        set(gca,'Box','off');
        ylabel({'Social Task' 'LH activation - RH activation'},'FontName','arial','FontWeight','bold','FontSize',12);
        xlabel({'LH activation - RH activation' 'Language Task'},'FontName','arial','FontWeight','bold','FontSize',12);
        plot(zeros(201,1), [-100:100], [':', 'k']);
        plot([-200:200], zeros(401,1), [':', 'k']);
        xlim([-200 200]);
        ylim([-100 100]);
        yticks([-100 -50 0 50 100]);
        yticklabels({'-100','-50','0','50','100'});
        axis square;
%         export_fig(['fig3_scatter_left_handed/' A{i} '.png'], '-r600');
    end
end

