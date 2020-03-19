function [SegAtrri, PtsAttri] = GraphOpt(Edge, Seg_final, SegAtrri, PtsAttri)
%% update segment level Pli with graph optimization
% Updated_Pli: updated Pli for segments
% Pboriginal: original point level pli 
% PbUpdate: updated point level pli 

%%
S = SegAtrri.S;
Pli = SegAtrri.Pli;


[t_label,~,ic] = unique(Edge(:,1));
t = accumarray(ic,Edge(:,2),[],@(x) {x});

Initial_Pli = Pli;
aa = 1;
ite = 0;
while aa>0.05 && ite<3
    
    Npli = nan(length(t_label),1);
    for i = 1:length(t_label)
        current_Seg = t_label(i);
        neighbor_Seg = t{i};
        
        all = [current_Seg;neighbor_Seg];
        Npli(i) = sum(S(all).*Pli(all))/sum(S(all));
    end
    
    Updated_Pli = Pli;
    Updated_Pli(t_label) = Npli;
    
    aa = mean(abs(Updated_Pli - Pli));
    
    Pli = Updated_Pli;
end

% upwrap segments to points, so we can operate on point level
Pboriginal = cell(length(Seg_final),1);
PbUpdate = cell(length(Seg_final),1);
for i = 1:length(Seg_final)
    dump = Seg_final{i};
    Pboriginal{i} = repmat(Initial_Pli(i),size(dump,1),1);
    PbUpdate{i} = repmat(Updated_Pli(i),size(dump,1),1);
end
Pboriginal = cell2mat(Pboriginal);
PbUpdate = cell2mat(PbUpdate);

% figure
% subplot(1,2,1)
% pcshow(P, Pboriginal)
% axis off
% subplot(1,2,2)
% pcshow(P, PbUpdate)
% axis off
% 
% clear Pli ic t t_label

SegAtrri.Updated_Pli = Updated_Pli;
PtsAttri.Pboriginal = Pboriginal;
PtsAttri.PbUpdate = PbUpdate;
end