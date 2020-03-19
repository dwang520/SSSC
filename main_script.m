
%% recursive segmentation
ft_threshold = 0.125;
nThreads = 10;

Seg_final = RcSeg(points, ft_threshold, nThreads);

%% visualize segmentation
cmap = hsv(length(Seg_final));
pz = randperm(size(cmap,1),size(cmap,1));
cmap = cmap(pz,:);
col = cell(length(Seg_final),1);
for i =1:length(Seg_final)
    ww = Seg_final{i};
    col(i) = {repmat(cmap(i,:),size(ww,1),1)};
end
figure;pcshow(cell2mat(Seg_final),cell2mat(col));grid off;
clear cmap pz col

%% calculate segment attribtues
[PtsAttri, SegAtrri] = Seg_attributes(Seg_final);

poolobj = gcp('nocreate');
delete(poolobj);
%% construct edges for supergraph
Edge = ConstructEdges(20, PtsAttri, SegAtrri);

%% construct a weighted graph
weight1 = sqrt(sum((SegAtrri.Lp(Edge(:,1),:) - SegAtrri.Lp(Edge(:,2),:)).^2,2));
% construct a graph based on edge list
adj = sparse(Edge(:,1),Edge(:,2),weight1,length(Seg_final),length(Seg_final));
Gp = graph(adj);
clear adj weight1 weight2 weight

%% updae segment Pli with local graph optimization
[SegAtrri, PtsAttri] = GraphOpt(Edge, Seg_final, SegAtrri, PtsAttri);
clear Edge

%% find roots
SegAtrri = FindRoots_1010(SegAtrri, Gp, 1, 0.8, 100, 1, 0.5);

% visual inspection
Roots = Seg_final(SegAtrri.Root_id);
cmap = hsv(length(Roots));
pz = randperm(size(cmap,1),size(cmap,1));
cmap = cmap(pz,:);
col = cell(length(Roots),1);
for i =1:length(Roots)
    ww = Roots{i};
    col(i) = {repmat(cmap(i,:),size(ww,1),1)};
end
figure;pcshow(cell2mat(Roots),cell2mat(col));grid off;
clear cmap pz col

ia = SegAtrri.Lp(:,3) <= 5;
dump = cell2mat(Seg_final(ia));
dump(ismember(dump,cell2mat(Roots),'rows'),:) = [];

hold on
pcshow(dump,[.3 .3 .3])
hold off

%% apply shortest path
DBH = ones(length(SegAtrri.Root_id),1);
% shortest path to each root
d = distances(Gp,SegAtrri.Root_id,'Method','positive');
% weight by DBH
w = DBH.^(2/3)';
d = d'./w;
d(d == inf) = nan;
[dump,SegAtrri.treeid] = min(d,[],2,'omitnan');
% [dump,SegAtrri.treeid] = min(d,[],2);
SegAtrri.treeid(isnan(dump)) = nan;

clear d

figure
pcshow(SegAtrri.C,SegAtrri.treeid)
colormap(colorcube(max(SegAtrri.treeid)))

%% Path optimization
SegAtrri = PathOpt(Gp, SegAtrri);

%% segment level tree id to point wise
Tid = cell(size(Seg_final,1),1);% point tree id
Wood = cell(size(Seg_final,1),1);% point wood/leaf
for i = 1:size(Seg_final,1)
    dump = Seg_final{i};
    Tid{i} = SegAtrri.treeid(i)*ones(size(dump,1),1);
    Wood{i} = SegAtrri.Final_Pli(i)*ones(size(dump,1),1);
end
PtsAttri.Treeid = cell2mat(Tid);
PtsAttri.Pbfinal = cell2mat(Wood);
clear Tid Wood dump
%% post processing
[PtsAttri, ~] = PostProssing(Seg_final, SegAtrri, PtsAttri, 0.13, 0.1);


%% results

results = [PtsAttri.P, PtsAttri.Treeid, double(PtsAttri.Pbfinal>=0.5)];

%%
figure
subplot(1,2,1)
pcshow(PtsAttri.P, PtsAttri.Treeid)
colormap(lines(max(PtsAttri.P)))
title('Tree Isolation','Color', 'k','FontWeight','Normal')
set(gcf,'color','w');
set(gca,'color','w');
set(gca, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'ZColor', [0.15 0.15 0.15])
axis off

subplot(1,2,2)
tmp = PtsAttri.Pbfinal>=0.5;
wood = PtsAttri.P(tmp,:);
leaf = PtsAttri.P(~tmp,:);

pcshow(wood, repmat([0.4471, 0.3216, 0.1647],size(wood,1),1),'markersize',1);
hold on
pcshow(leaf, repmat([0.2667, 0.5686, 0.1961],size(leaf,1),1),'markersize',1);
hold off
title('Leaf-wood Separation','Color', 'k','FontWeight','Normal')
set(gcf,'color','w');
set(gca,'color','w');
set(gca, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'ZColor', [0.15 0.15 0.15])
axis off