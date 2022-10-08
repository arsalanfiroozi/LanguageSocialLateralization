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

% [val,pos]=intersect(files_social, files_language);


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
%%
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
g2 = gifti('Fig4a_CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR_L.func.gii');
A = g2.cdata;
C = zeros(size(g2.cdata));

Cole_Data = g2.cdata;
Ls = mode(Cole_Data(g1.cdata == (129 + 180))); % 3a

selected_index_Lan = 1:180;
selected_index_Lan = selected_index_Lan(abs(mean(RActivity_Lan,2)) >= Lan_Sorted(Lan_x));

for i=1:length(B)
    C(g1.cdata == (B(i) + 180)) = Lan_Diff(B(i));
end

g3 = gifti('Fig1e.func.gii');
g3.cdata = C;
save(g3,'All_ROIs_Lan_2.func.gii');
%%
g2 = gifti('Fig4a_CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR_L.func.gii');
A = g2.cdata;
C = zeros(size(g2.cdata));

Cole_Data = g2.cdata;
Ls = mode(Cole_Data(g1.cdata == (129 + 180))); % 3a

selected_index_Soc = 1:180;
selected_index_Soc = selected_index_Soc(abs(mean(RActivity_Soc,2)) >= Soc_Sorted(Soc_x));
% T = (1:180);
% selected_index_Lan = selected_index_Lan((T == 161) | (T == 67) | (T == 70) | (T == 87) | (T == 171));

for i=1:length(B)
    C(g1.cdata == (B(i) + 180)) = Soc_Diff(B(i));
end

g3 = gifti('Fig1e.func.gii');
g3.cdata = C;
save(g3,'All_ROIs_Soc_2.func.gii');