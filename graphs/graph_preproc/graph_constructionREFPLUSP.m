function CM = graph_constructionREFPLUSP(A, method, referenceGroup, similarityMeasure)
% Jackknife Bias Estimation for network construction (Das et al., 2018)
%
% INPUT: matrix format, dimensions: n people x m features 
%ƒ
% Steps: 
%   1. compute reference network from reference group provided
%   (referenceGroup)
%   2. perform jackknife bias estimation for each individual, i.e. 
%   reestimate network for reference group + individual and then calculate 
%   the difference = individual network
%   3. bring data in the right format (flat connectivity matrix, upper
%       triangle); dimensions: n people x (m/2-m) edges 
%  
%   --> The different variable types of a pair of variables have to be   
%   taken in consideration as they determine which metric will be used to 
%   quantify the relationship between the respective pair of variables
%       - continuous (gaussean), continuous (gaussean): Pearson correlation
%       - continuous (gaussean/ not gaussean), continuous (not gaussean): 
%       Spearman correlation
%       - continuous (gaussean/ not gaussean), categorical (ordinal):
%       Spearman correlation
%       - categorical (ordinal), categorical (ordinal): Spearman
%       correlation
%       - continuous, categorical (2 levels): point-biseral 
%       correlation (equivalent to Pearson correlation)
%       - continuous, categorical (>2 levels): eta correlation
%       - categorical, categorical: Chi-Square/ Cramer's V (?)
%
%
% OUTPUT: matrix format, dimensions: n people x (m/2-m) edges 
global NM
% % read in reference group data 
if size(referenceGroup,1) >1
    RG = referenceGroup; 
else
    RG = readtable(referenceGroup);
    RG = table2array(RG);
end

% if size(RG,2) == size(A,2) % check whether ID and case label are columns 

% else
%     % remove case and label variable 
%     label_var = NM.datadescriptor{1,1}.input_settings.label_edit;
%     case_var = NM.datadescriptor{1,1}.input_settings.case_edit;
%     RG = table2array(removevars(RG, {label_var, case_var}));
% end
% RG = referenceGroup;
% read in variable types vector
%varTypes = readtable(variableTypesVec);

% compute reference group correlation matrix (for now only continuous 
% variables) 
% if ~ismember('c',varTypes)
%     
% end
switch similarityMeasure
    case 'Pearson correlation'
        RCM = corr(RG, 'Type', 'Pearson');
        for i = 1:size(A,1)
            CURSAMPLE = [RG;A(i,:)];
            CURCM = corr(CURSAMPLE, 'Type', 'Pearson');
            DIFCM = RCM - CURCM;
            UPPTRI = triu(DIFCM,1);
            tri  = triu(true(size(DIFCM)),1);
            if ~exist('R', 'var')
                R = UPPTRI(tri)';
            else
                R = [R;UPPTRI(tri)'];
            end
        end
    case 'Spearman`s rho'
        RCM = corr(RG, 'Type', 'Spearman');
        for i = 1:size(A,1)
            CURSAMPLE = [RG;A(i,:)];
            CURCM = corr(CURSAMPLE, 'Type', 'Spearman');
            DIFCM = RCM - CURCM;
            UPPTRI = triu(DIFCM,1);
            tri  = triu(true(size(DIFCM)),1);
            if ~exist('R', 'var')
                R = UPPTRI(tri)';
            else
                R = [R;UPPTRI(tri)'];
            end
        end
    case 'Kendall`s Tau'
        RCM = corr(RG, 'Type', 'Kendall');
        for i = 1:size(A,1)
            CURSAMPLE = [RG;A(i,:)];
            CURCM = corr(CURSAMPLE, 'Type', 'Kendall');
            DIFCM = RCM - CURCM;
            UPPTRI = triu(DIFCM,1);
            tri  = triu(true(size(DIFCM)),1);
            if ~exist('R', 'var')
                R = UPPTRI(tri)';
            else
                R = [R;UPPTRI(tri)'];
            end
        end
    case 'Mutual information'
        R = pyrunfile('cv_computeMInetworks.py', 'netwX', ...
            ref = RG, ...
            X = A);
end
CM = double(R);
end

function MIM = mi_mat(DATA)
MIM = zeros([size(DATA,2), size(DATA,2)]);
bins = 10; 
    for i = 1:size(DATA,2)
        for j = 1:size(DATA,2)
            if i<=j
                a = discretize(DATA(:,i));
                b = discretize(DATA(:,j));
                cur_mi = mi(a,b);
                
%                 cur_mi = pyrunfile('py_calc_MI.py', 'py_mi', ...
%                     veci = DATA(:,i), vecj = DATA(:,j), ...
%                     nbins = int64(bins));
                MIM(i,j) = cur_mi;
                MIM(j,i) = cur_mi;
            end
        end
    end
end
