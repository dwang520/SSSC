function Edge = ConstructEdges(k, PtsAttri, SegAtrri)
%% construct a superpoint graph
% this graph is firstly constructed on points P, and then transfered to
% superpoints
% k: number of neighbors for knn graph
% P: points
% L: per point segment label
% SegAtrri: segment attributes
%% construct a connected knn graph for points
% point wise
% k  = 20;
P = PtsAttri.P;
L = PtsAttri.L;
C = SegAtrri.C;
%%
n_point = size(P,1);
%---compute full adjacency graph-------------------------------------------
[neighbors,~] = knnsearch(P,P,'K', k+1);
% remove point itself
neighbors = neighbors(:,2:end);
% convert to nx2 matrix
source = reshape(repmat(1:n_point, [k 1]), [1 (k * n_point)])';
target = reshape(neighbors', [1 (k * n_point)])';
clear neighbors distance prune2 pruned

%% -----transform point wise to segment wise------------------------------
S_source = L(source);
S_target = L(target);

% remove those with very long horizontal distances
dump = C(S_source,1:2) - C(S_target,1:2);
dis = sqrt(sum(dump.^2,2));
disthres = prctile(dis,99.9);
S_source(dis > disthres) = [];
S_target(dis > disthres) = [];
clear dump dis
% remove redundant edges
selfedge = S_source==S_target;
S_source = S_source(~selfedge);
S_target = S_target(~selfedge);

Edge = [S_source,S_target];
Edge = sort(Edge,2);
Edge = unique(Edge, 'rows');
% symmetric
Edge2 = [Edge(:,2),Edge(:,1)];
Edge = [Edge;Edge2];

clear Edge2 source target S_source S_target selfedge

end