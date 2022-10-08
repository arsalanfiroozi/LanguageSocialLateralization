%% Lan
clearvars
ROI_Parcellation

addpath(genpath('matlabGiftiCifti\'));

g = gifti('SFig15.R.func.gii');
g = g.cdata;

Names = {'Social Cognition' 'Incentive Processing' 'Emotion Processing' 'Relational Processing' 'Working Memory' 'Face Localizer' 'Place Localizer' 'Body Localizer' 'Tool Localizer'};
m = {g(:,74), g(:,36), g(:,83), g(:,78), g(:,11), g(:,20), g(:,21), g(:,19), g(:,22)};

A = {'STSp' 'STSa' 'STGa' 'PGi' 'PSL' 'Broca' '55b' 'SFL'};
values = [-9 41 11 -29 -19 -49 31 21];

res_L(sum(res_L==values,2)==1) = 1;

data = zeros(9,1);
err = zeros(9,1);
for i=1:9
    t = m{i};
    data(i) = mean(t(res_L==1));
    err(i) = std(t(res_L==1))/sqrt(sum(res_L==1));
end

f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
bar(1:9, data)
hold on
e = errorbar(1:9, data, err);
e.LineStyle = 'none';
e.Color = [0 0 0];
box off
xticklabels(Names)
ylabel("Cohen's d");
ylim([-0.3 0.7])
% export_fig('Rev07_2_barplot.png','-r600');