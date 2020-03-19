function Seg_final = RcSeg(points, ft_threshold, nThreads)
%% recursive segmentation



%%
normals = pcnormals(pointCloud(points(:,1:3)),10);
% ft_threshold = 0.125;
k = 10;
Label = GraphRG(points, normals(:,3), k, ft_threshold);

t = accumarray(Label,[1:length(Label)]',[],@(x) {x});
Seg = cellfun(@(x) points(x,:),t,'UniformOutput',0);
lb = cellfun(@(x) size(x,1),Seg);

Seg_reserve = Seg(lb<=10);
Seg_processing = Seg(lb>10);
clear Seg t lb Label normals

p = gcp('nocreate');
if isempty(p)
    parpool('local',nThreads);
end
%% recursive segmentation
for ii = 1:inf
    
    TMP = cell(length(Seg_processing),1);
    indx = false(length(Seg_processing),1);
    parfor i = 1:length(Seg_processing)
        points_n = Seg_processing{i};
        normals = pcnormals(pointCloud(points_n(:,1:3)),10);
        Label = GraphRG(points_n, normals(:,3), k, ft_threshold);
        
        if max(Label)==1
            indx(i) = true;
        else
            t = accumarray(Label,[1:length(Label)]',[],@(x) {x});
            Seg = cellfun(@(x) points_n(x,:),t,'UniformOutput',0);
            
            TMP{i} = Seg;
        end
    end
    
    Seg = cat(1, TMP{:});
    
    unchange = Seg_processing(indx);
    
    if sum(indx) == length(Seg_processing) || ii == 10 % very unlikely the iteration will exceed 10
        Seg_final = [Seg;Seg_reserve;unchange];
        break
    else
        lb = cellfun(@(x) size(x,1),Seg);
        
        Seg_processing = Seg(lb>10);
        
        Seg_reserve = [Seg_reserve;unchange;Seg(lb<=10)];
    end
end

clear Seg_reserve unchange Seg_processing Seg
%%
% cmap = hsv(length(Seg_final));
% pz = randperm(size(cmap,1),size(cmap,1));
% cmap = cmap(pz,:);
% col = cell(length(Seg_final),1);
% for i =1:length(Seg_final)
%     ww = Seg_final{i};
%     col(i) = {repmat(cmap(i,:),size(ww,1),1)};
% end
% figure;pcshow(cell2mat(Seg_final),cell2mat(col));grid off;
% clear cmap pz col
%%
Inputs.PatchDiam1 = 0.1; % Patch size of the first uniform-size cover
Inputs.PatchDiam2Min = 0.03; % Minimum patch size of the cover sets in the second cover
Inputs.PatchDiam2Max = 0.08; % Maximum cover set size in the stem's base in the second cover
Inputs.lcyl = 3; % Relative (length/radius) length of the cylinders
Inputs.FilRad = 3; % Relative radius for outlier point filtering
Inputs.BallRad1 = Inputs.PatchDiam1+0.02; % Ball radius in the first uniform-size cover generation
Inputs.BallRad2 = Inputs.PatchDiam2Max+0.01; % Maximum ball radius in the second cover generation
Inputs.nmin1 = 3; % Minimum number of points in BallRad1-balls, generally good value is 3
Inputs.nmin2 = 1; % Minimum number of points in BallRad2-balls, generally good value is 1
Inputs.OnlyTree = 1; % If 1, point cloud contains points only from the tree

%% slipt segments
% addpath('./TreeQSM_src')
% linearity = Cal_Linearity(Seg_final);
Seg_post = cell(length(Seg_final),1);
parfor i = 1:length(Seg_final)
    
    pts = Seg_final{i};
    if size(pts,1)>100 && (max(pts(:,3)) - min(pts(:,3)))>1
        % Generate cover sets
        cover1 = cover_sets(pts,Inputs);
        if length(cover1.ball) > 2
            Base = find(pts(cover1.center,3) == min(pts(cover1.center,3)));
            Forb = false(length(cover1.center),1);
            segment1 = segments(cover1,Base,Forb);
            Seg_tmp = cell(length(segment1.segments),1);
            for j = 1:length(segment1.segments)
                ids = cell2mat(cover1.ball(cell2mat(segment1.segments{j})));
                Seg_tmp{j} = pts(ids,:);
            end
            Seg_post{i} = Seg_tmp;
        else
            Seg_post{i} = {pts};
        end
    else
        Seg_post{i} = {pts};
    end
end
Seg_final= cat(1, Seg_post{:});
clear Seg_post


end