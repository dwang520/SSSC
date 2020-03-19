function SegAtrri = FindRoots_1010(SegAtrri, Gp, Lp_thres, direction_thres, S_thres, H_thres, Merge_thres)
%% find root segments


Updated_Pli = SegAtrri.Updated_Pli;
Lp = SegAtrri.Lp;
direction = SegAtrri.direction;
S = SegAtrri.S;
H = SegAtrri.H;
Hp = SegAtrri.Hp;

%% find roots
% wood segs
ia = find(Updated_Pli >= 0.5);
% % candidate root segs
ia(Lp(ia,3)>Lp_thres | abs(direction(ia,3))<direction_thres | S(ia)<S_thres | H(ia)<H_thres ) = [];

%% graph merge
% travel steps bwteen any of these segs
d = distances(Gp,ia,ia,'Method','unweighted');
d(d ==0 ) = inf;

[mind,ib] = min(d,[],2);
% merge those with very short steps
ic = find(mind <= 1);
if isempty(ic)~=1
    % connected segs
    Cs = [ic,ib(ic)];
    Cs = [Cs;Cs(:,[2,1])];
    
    adj2 = sparse(Cs(:,1),Cs(:,2),ones(size(Cs,1),1),length(ia),length(ia));
    G2 = graph(adj2);
    bins = conncomp(G2);
    % merge conncted segs
    t = accumarray(bins',ia,[],@(x) {x});
    
    Rt2 = nan(length(t),1);
    for i = 1:length(t)
        css = t{i};
        if length(css)~=1
            Rt2(i) = css(find(Lp(css,3) == min(Lp(css,3)),1));
        else
            Rt2(i) = css;
        end
    end
    
else
    Rt2 = ia;
end

%% furhter merge segments
% sort
[~, ic] = sort(Hp(Rt2,3),'ascend');
Rt2 = Rt2(ic);

merge = cell(length(Rt2),1);
mm = 1;
while isempty(Rt2) ~= 1
    
    % highest point of current segment
    Tp = Hp(Rt2(1),:);
    Mg = Rt2(1);
    while isempty(Tp) ~= 1
        % horizontal distance between lowest points to current highest point
        Hdis = sqrt(sum((Lp(Rt2,1:2) - Tp(1:2)).^2,2));
        %
        Pclust = find(Hdis<=Merge_thres & Lp(Rt2,3)>=Tp(3)*0.8);
        Pclust(ismember(Rt2(Pclust),Mg)) = [];
        
        % find nearest one
        nrs = Pclust(find(Hdis(Pclust) == min(Hdis(Pclust)),1));
        % merge this one
        Mg = [Mg;Rt2(nrs)];
        % update highest point
        Tp = Hp(Rt2(nrs),:);
    end
    Mg = unique(Mg);
    
    merge{mm} = Mg;
    mm = mm+1;
    Rt2(ismember(Rt2,Mg)) = [];
end

merge(cellfun(@isempty,merge)) = [];

%% Final root segments
Root_id = nan(length(merge),1);
for i = 1:length(merge)
    tmp = merge{i};
    if length(tmp)>1
        
        Root_id(i) = tmp(find(Lp(tmp,3) == min(Lp(tmp,3)),1));
    else
        Root_id(i) = tmp;
    end
end

%%
SegAtrri.Root_id = Root_id;
end
