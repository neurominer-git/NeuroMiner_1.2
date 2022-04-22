function CM = graph_constructionJBE(A, method, referenceGroup, variableTypesVec)
% Jackknife Bias Estimation for network construction (Das et al., 2018)
%
% INPUT: matrix format, dimensions: n people x m features 
%
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
% read in reference group data 
RG = readtable(referenceGroup);

% remove case and label variable 
label_var = NM.datadescriptor{1,1}.input_settings.label_edit;
case_var = NM.datadescriptor{1,1}.input_settings.case_edit;
RG = table2array(removevars(RG, {label_var, case_var}));

% read in variable types vector
varTypes = readtable(variableTypesVec);

% compute reference group correlation matrix (for now only continuous 
% variables) 
% if ~ismember('c',varTypes)
%     
% end
RCM = corrcoef(RG);
% compute individual networks 
%R = zeros(size(A));
for i = 1:size(A,1)
    CURSAMPLE = [RG;A(i,:)];
    CURCM = corrcoef(CURSAMPLE);
    DIFCM = RCM - CURCM;
    UPPTRI = triu(DIFCM,1);
    tri  = triu(true(size(DIFCM)),1);
    if ~exist('R', 'var')
        R = UPPTRI(tri)';
    else 
        R = [R;UPPTRI(tri)'];
    end
end
CM = R;
end