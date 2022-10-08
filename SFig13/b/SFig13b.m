%%
clearvars
addpath(genpath('matlabGiftiCifti\'));
B_Size = 20;
load('SFig13b_Files.mat');
Names_Files = Files(:,1);
id_Files = Files(:,2);
Mental_Files = Names_Files(strcmp(id_Files,'Mental'));
Random_Files = Names_Files(strcmp(id_Files,'Random'));

%% Social
D_Social = zeros(2,2);
for j=1:length(Mental_Files)/2
    videoFReader = VideoReader(['SOCIAL Stimuli/' Mental_Files{j}]);
    numFrames = floor(videoFReader.FrameRate * 18);
%     j
%     videoFReader.width
    B_Size_h = 768/2;
    B_Size_w = 1024/2;
    C_Contrast = zeros(768/B_Size_h,1024/B_Size_w);
    for i= 1:numFrames
        Frame = ones(768,1024) * 128;
        Frame((-videoFReader.height/2 + 768/2):(+videoFReader.height/2 + 768/2 - 1),(-videoFReader.width/2 + 1024/2):(+videoFReader.width/2 + 1024/2 - 1)) = rgb2gray(readFrame(videoFReader));
    
        fun = @(block_struct) std2(block_struct.data);
        B = blockproc(Frame,[B_Size_h B_Size_w],fun);
        C_Contrast = C_Contrast + B;
    end
    C_Contrast = C_Contrast ./ numFrames;
    D_Social = D_Social + C_Contrast;
    size(C_Contrast)
end
D_Social = D_Social ./ 10;

%% Random
D_Random = zeros(2,2);
for j=1:length(Random_Files)/2
    videoFReader = VideoReader(['SOCIAL Stimuli/' Random_Files{j}]);
    numFrames = floor(videoFReader.FrameRate * 18);
%     j
%     videoFReader.width
    B_Size_h = 768/2;
    B_Size_w = 1024/2;
    C_Contrast = zeros(768/B_Size_h,1024/B_Size_w);
    for i= 1:numFrames
        Frame = ones(768,1024) * 128;
        Frame((-videoFReader.height/2 + 768/2):(+videoFReader.height/2 + 768/2 - 1),(-videoFReader.width/2 + 1024/2):(+videoFReader.width/2 + 1024/2 - 1)) = rgb2gray(readFrame(videoFReader));
    
        fun = @(block_struct) std2(block_struct.data);
        B = blockproc(Frame,[B_Size_h B_Size_w],fun);
        C_Contrast = C_Contrast + B;
    end
    C_Contrast = C_Contrast ./ numFrames;
    D_Random = D_Random + C_Contrast;
    size(C_Contrast)
end
D_Random = D_Random ./ 10;

%%
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
imagesc(D_Social,[min([min(D_Social(:)) min(D_Random)]) max([max(D_Social(:)) max(D_Random)])]);
title('Social');
set(gca,'YTick',[]);
set(gca,'XTick',[]);
colorbar;
% export_fig('Fig_Social.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
imagesc(D_Random);
title('Random');
set(gca,'YTick',[]);
set(gca,'XTick',[]);
colorbar;
% export_fig('Fig_Random.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
imagesc(D_Social - D_Random);
title('Difference');
set(gca,'YTick',[]);
set(gca,'XTick',[]);
colorbar;
% export_fig('Fig_Difference.png','-r600');
