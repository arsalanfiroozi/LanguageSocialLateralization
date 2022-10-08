%%
clearvars;
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
% % 
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
%% Lan
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


figure;
bar(abs(mean(RActivity_Lan,2)));
add_errorbar(Error_Lan', abs(mean(RActivity_Lan,2)));
xticks(1:180);
set(gca,'XTickLabel', Names);
title('Language');

figure;
bar(abs(mean(RActivity_Soc,2)));
add_errorbar(Error_Soc', abs(mean(RActivity_Soc,2)));
xticks(1:180);
set(gca,'XTickLabel', Names);
title('Social');

Tmp = Lan_Label;
Tmp(:) = {''};
figure;
hold on;
set(gcf,'Color',[1 1 1]);
set(gca,'TickLength',[0.005, 0.01])
set(gca,'TickDir','out');
set(gca,'FontName','arial','FontSize',10); % Check this
h = bar(flip(Lan_Sorted),'FaceColor',[.5 .5 .5]);
YLim = ylim;
plot(ones(YLim(2) - YLim(1) + 1,1)*(360-2*Lan_x+3)/2, [YLim(1):YLim(2)], [':', 'k']);
bar(flip(Lan_Sorted .* (Lan_Sorted >= Lan_Sorted(Lan_x))),'FaceColor',[1 0 0]);
add_errorbar(flip(Lan_Error)', flip(Lan_Sorted));
xticks(1:180);
set(gca,'XTickLabel', flip(Tmp));
set(gca,'xtick',[]);
% title('Language');
ylabel({'Language task' '|LH activation - RH activation|'},'FontWeight','bold');
% export_fig('Lan_180barplot.png','-r600');