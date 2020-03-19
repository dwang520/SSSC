function SegAtrri = PathOpt(Gp, SegAtrri)


treeids = unique(SegAtrri.treeid(~isnan(SegAtrri.treeid)));

visited_paths = cell(length(treeids),1);
for i = 1:length(treeids)
    ia = find(SegAtrri.treeid == treeids(i));
    Pia = SegAtrri.Updated_Pli(ia);
    Wia = ia(Pia>=0.8); % wood segments connected to Root_id(treeids(i))
    
    if isempty(Wia)~=1
        TR = shortestpathtree(Gp,Wia,SegAtrri.Root_id(treeids(i)),'OutputForm','cell');
        % redistribution
        path_u = cell(length(TR),1);
        for j = 1:length(TR)
            dump = TR{j};
            wsn = sum(SegAtrri.Updated_Pli(dump)>=0.5);
            dump = flip(dump);
            path_u{j} = dump(1:wsn);
        end
        visited_paths{i} = cell2mat(path_u')';
    end
end
visited_paths = cell2mat(visited_paths);

tbl = tabulate(visited_paths);

Seg_Visit_Freq = zeros(size(SegAtrri.C,1),1);
if isempty(tbl) ~= 1
    Seg_Visit_Freq(tbl(:,1)) = tbl(:,2);
end
%% apply graph path optimization
ia = find(Seg_Visit_Freq>0);
ib = ia(SegAtrri.Updated_Pli(ia)<0.5);
% ****** final segment level wood/leaf probability ********
Final_Pli = SegAtrri.Updated_Pli;
Final_Pli(ib) = 0.51;

%%
SegAtrri.Final_Pli = Final_Pli;

end