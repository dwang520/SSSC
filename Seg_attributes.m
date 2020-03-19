function [PtsAttri, SegAtrri] = Seg_attributes(Seg_final)
%% calculate segments attributes
% SegAtrri.C: segment center point - superpoints
% SegAtrri.S: segment size
% SegAtrri.Lp: segment lowest point
% SegAtrri.Hp: segment highest point
% SegAtrri.: segment length
% SegAtrri.direction:segment orientation
% SegAtrri.Pli: segment wood probability
%% also export point attributes
% PtsAttri.P: unwrapped points
% PtsAttri.L: per point segment label
%% ------------------------------
% calculate linearity of each segment
linearity = Cal_Linearity(Seg_final);
% calculate size of each segment
sl = cellfun(@(x)size(x,1),Seg_final);

% test a range of thresholds
Lthres_list = 0.70:0.02:0.95;
Sthres_list = 10:2:50;
% all combinations of two thresholds
allc = combvec(Lthres_list,Sthres_list)';
%
Freq = zeros(length(Seg_final),1);
for i = 1:size(allc,1)
    % find wood segments based on two thresholds
    ia = linearity >= allc(i,1) & sl >= allc(i,2);
    % count the frequency of being identified as wood
    Freq = Freq + ia;
end

% Probability that a segment is wood !!!
Pli = Freq/size(allc,1);

clear linearity sl Lthres_list Sthres_list allc Freq ia
%% segment features
C = nan(size(Seg_final,1),3); % segment center point - superpoints
L = cell(size(Seg_final,1),1);% point segment label
S = nan(size(Seg_final,1),1); % segment size
Lp = nan(size(Seg_final,1),3); % segment lowest point
Hp = nan(size(Seg_final,1),3); % segment highest point
H = nan(size(Seg_final,1),1); % segment length 
direction = zeros(length(Seg_final),3); %segment orientation
parfor i = 1:size(Seg_final,1)
    dump = Seg_final{i};
    ia = size(dump,1);
    S(i) = ia;
    H(i) = max(dump(:,3)) - min(dump(:,3));
    L{i} = i*ones(ia,1);
    
    dis = sqrt(sum((dump - mean(dump,1)).^2,2));
    C(i,:) = dump(find(dis == min(dis),1),:);
    
    Lp(i,:) = dump(find(dump(:,3) == min(dump(:,3)),1),:);
    Hp(i,:) = dump(find(dump(:,3) == max(dump(:,3)),1),:);
    
    if size(dump,1)>3
        coeff = pca(dump);
        direction(i,:) = coeff(:, 1);
    end
end

L = cell2mat(L);
P = cell2mat(Seg_final); % point cloud

% poolobj = gcp('nocreate');
% delete(poolobj);

%% write to a structure 
SegAtrri.C = C;
SegAtrri.S = S;
SegAtrri.Lp = Lp;
SegAtrri.H = H;
SegAtrri.Hp = Hp;
SegAtrri.direction = direction;
SegAtrri.Pli = Pli;

PtsAttri.P = P;
PtsAttri.L = L;
end