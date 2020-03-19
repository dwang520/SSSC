function [PtsAttri, SegAtrri] = PostProssing(Seg_final, SegAtrri, PtsAttri, lthres, dthres)
% new version
%% remove low vegetations

ia = find(PtsAttri.P(:,3)<lthres);
rootpoints = cell2mat(Seg_final(SegAtrri.Root_id));
rootpoints(rootpoints(:,3)>lthres,:) = [];

[~,D] = knnsearch(rootpoints(:,1:2),PtsAttri.P(ia,1:2));

% vegetation points
vegep = ia(D>dthres); 

PtsAttri.Treeid(vegep) = 0;

%%
%% deal with unconnected points (to nearest valid points)
ia = isnan(PtsAttri.Treeid);
unlabeledP = PtsAttri.P(ia,:);
restP = PtsAttri.P(~ia,:);
restTid = PtsAttri.Treeid(~ia);

idx = knnsearch(restP, unlabeledP);
PtsAttri.Treeid(ia) = restTid(idx);

%% old version
% % figure
% % pcshow(PtsAttri.P,PtsAttri.Treeid)
% % colormap(gca,[colorcube(max(PtsAttri.Treeid))])
% % grid off
% % figure
% % pcshow(PtsAttri.P,double(PtsAttri.Pbfinal>=0.5))
% % colormap(gca,'parula')
% % grid off
% 
% %% deal with unconnected points (to nearest valid points)
% ia = isnan(PtsAttri.Treeid);
% unlabeledP = PtsAttri.P(ia,:);
% restP = PtsAttri.P(~ia,:);
% restTid = PtsAttri.Treeid(~ia);
% 
% idx = knnsearch(restP, unlabeledP);
% PtsAttri.Treeid(ia) = restTid(idx);
% 
% 
% % ia = isnan(PtsAttri.Treeid);
% % unlabeledP = PtsAttri.P(ia,:);
% % 
% % idx = knnsearch(SegAtrri.Lp(SegAtrri.Root_id,1:2), unlabeledP(:,1:2));
% % PtsAttri.Treeid(ia) = idx;
% % clear idx ia unlabeledP
% 
% 
% 
% %% remove low vegetations
% 
% ia = find(PtsAttri.P(:,3)<lthres);
% rootpoints = cell2mat(Seg_final(SegAtrri.Root_id));
% rootpoints(rootpoints(:,3)>lthres,:) = [];
% 
% [~,D] = knnsearch(rootpoints(:,1:2),PtsAttri.P(ia,1:2));
% 
% % vegetation points
% vegep = ia(D>dthres); 
% 
% PtsAttri.Treeid(vegep) = 0;
% 
% % figure
% % pcshow(PtsAttri.P,PtsAttri.Treeid)
% % colormap(colorcube(max(PtsAttri.Treeid)))
% % grid off

end