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

g3 = gifti('ExtractedData/990366_Lan.L.func.gii');

if(~sum(checklist1 - checklist2))
    percent = 100:-1:1;
    A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
    values = [-9 41 11 -29 -19 -49 31 21];
    DiffL_mean = mean(DiffL,2);
    DiffS_mean = mean(DiffS,2);
    dice_value=zeros(size(percent,2),size(A,2));
    jaccard_value=zeros(size(percent,2),size(A,2));
    
    for i=1:size(percent,2)
        for j=1:size(A,2)
            selected_vertices_roi = res_L == values(j);
            selected_vertices_L = selected_vertices_roi .* DiffL_mean;
            selected_vertices_S = selected_vertices_roi .* DiffS_mean;
            val_L = selected_vertices_L(selected_vertices_L > 0);
            val_S = selected_vertices_S(selected_vertices_S < 0);
            t_L = quantile(val_L, 1-percent(i)/100);
            t_S = quantile(val_S, percent(i)/100);
            data_L = selected_vertices_L >= t_L;
            data_S = selected_vertices_S <= t_S;
            dice_value(i, j) = dice(data_L,data_S);
            jaccard_value(i, j) = jaccard(data_L,data_S);
            if( mod(percent(i) ,10) == 0)
                g3.cdata = data_L;
%                 save(g3,['fig2_percentile/' num2str(percent(i))  '_Lan_' A{j} '.func.gii']);
                g3.cdata = data_S;
%                 save(g3,['fig2_percentile/' num2str(percent(i))  '_Soc_' A{j} '.func.gii']);
            end
        end
    end
    
%     A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
%     values = [-9 41 11 -29 -19 -49 31 21];
%     tmp(res_L == 1) = 0;
%     tmp(res_L == -49) = 11;
%     tmp(res_L == 21) = 8;
%     tmp(res_L == 31) = 9;
%     tmp(res_L == 41) = 16;
%     tmp(res_L == -9) = 2;
%     tmp(res_L == -29) = 5;
%     tmp(res_L == -19) = 6;
%     tmp(res_L == -39) = 0;
%     tmp(res_L == 11) = 12;

%   STGa ==> 7
    C = [[255 0 0];[246 132 9];[255 211 0];[0 149 246];[0 255 0];[246 158 220];[211 18 255];[106 18 141]]/255;
    
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'Box','on');
    set(gca,'FontName','arial','FontSize',10); % Check this
    hold on;
    for i=1:size(A,2)
        plot(percent, jaccard_value(:,i),'Color', C(i,:),'LineWidth',1);
    end
    legend(A,'Location','northwest','FontName','arial','FontSize',10);
    legend boxoff;
    axis square;
    ylabel({'Spatial overlap' 'between language and social asymmetry maps'},'FontName','arial','FontWeight','bold','FontSize',12);
    xlabel({'Percentage of vertices' 'with the highest asymmetry values'},'FontName','arial','FontWeight','bold','FontSize',12);
    plot(30*ones(201,1), 0:1/200:1, [':', 'k'], 'HandleVisibility','off');
    xticks([0 20 40 60 80 100]);
    xticklabels({'0','20','40','60','80','100'});
    yticks([0 0.2 0.4 0.6 0.8 1]);
    yticklabels({'0','0.2','0.4','0.6','0.8','1'});
    export_fig('fig2_Jaccardlineplot.png','-r600');
    
end
