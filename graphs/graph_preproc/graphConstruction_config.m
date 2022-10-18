function [GRAPHCONSTRUCTION, PX, act ] = graphConstruction_config(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)
global NM 
% what method should be used for the construction of the networks?
ConstructionMethod = 'Normative network + 1';
% ParcellationAtlas = 'Hammers.nii';
% VariableTypesVec = '';
ReferenceGroup = '';
SimilarityMeasure = 'Mutual information';

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    if isempty(GRAPHCONSTRUCTION), [GRAPHCONSTRUCTION, PX] = graphConstruction_config(GRAPHCONSTRUCTION, PX, parentstr, true); end
    if isfield(GRAPHCONSTRUCTION,'method'), ConstructionMethod = GRAPHCONSTRUCTION.method; end
    if isfield(GRAPHCONSTRUCTION,'parcellation'), ParcellationAtlas = GRAPHCONSTRUCTION.parcellation; end
    %if isfield(GRAPHCONSTRUCTION, 'variableTypes'), VariableTypesVec = GRAPHCONSTRUCTION.variableTypes; end
    if isfield(GRAPHCONSTRUCTION, 'refGroup'), ReferenceGroup = GRAPHCONSTRUCTION.refGroup; end
    if isfield(GRAPHCONSTRUCTION,'simMeasure'), SimilarityMeasure = GRAPHCONSTRUCTION.simMeasure; end
    % BETTER, as some of the options are only relevant for some
    % construction methods, to open an additional menu after method is
    % selected
    if ~ischar(ConstructionMethod)
        GCSTR_METHOD = 'undefined';
    else
        GCSTR_METHOD = ConstructionMethod;
    end

    if strcmp(GCSTR_METHOD,'KL divergence')
        if ~exist('ParcellationAtlas', 'var') || isempty(ParcellationAtlas)
            GCSTR_PARC = '(undefined)';
        else
            GCSTR_PARC = ParcellationAtlas;
        end
        menustr = ['Select graph construction method [ ' GCSTR_METHOD ' ]|', ...
            'Define parcellation atlas [ ' GCSTR_PARC ' ]'];
        menuact = [1,4];
    end

    if strcmp(GCSTR_METHOD, 'Normative network + 1')
        if isempty(SimilarityMeasure)
            GCSTR_SIMMEASURE = '(undefined)';
        else
            GCSTR_SIMMEASURE = SimilarityMeasure;
        end

        if isempty(ReferenceGroup)
            GCSTR_REFG = '(undefined)';
        else
            if size(ReferenceGroup,1) >1
                GCSTR_REFG = 'from Matlab workspace';
            else
                GCSTR_REFG = ReferenceGroup;
            end
        end
        
%         if isempty(VariableTypesVec)
%             GCSTR_VARTVEC = '(undefined)';
%         else
%             GCSTR_VARTVEC = sprintf('vector length: %d', length(VariableTypesVec));
%         end

        menustr = ['Select graph construction method [ ' GCSTR_METHOD ' ]|', ...
            'Define similarity measure [' GCSTR_SIMMEASURE ']|', ...
            'Define reference group [' GCSTR_REFG ']'];
            % 'Define variable types (vector) [ ' GCSTR_VARTVEC ']'] , ...
             % only important if we want to implement a network estimation similar to mgm 
        menuact = 1:3;
    end

    %[menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);

    nk_PrintLogo
    mestr = 'Graph construction'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr);
    act = nk_input(mestr,0,'mq', menustr, menuact);

    switch act
        case 1
            [GRAPHCONSTRUCTION,PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr);
            %ConstructionMethod = 'KL divergence';
        case 2
             simMeasureNo = nk_input('Define measure to quantify similarity between variables', 0, 'mq', ...
                ['Mutual information |', ...
                'Pearson correlation |', ...
                'Spearman`s rho |',...
                'Kendall`s tau'], 1:4, 0);
             switch simMeasureNo 
                 case 1
                     GRAPHCONSTRUCTION.simMeasure = 'Mutual information';
                 case 2
                     GRAPHCONSTRUCTION.simMeasure = 'Pearson correlation';
                 case 3
                     GRAPHCONSTRUCTION.simMeasure = 'Spearman`s rho';
                 case 4 
                     GRAPHCONSTRUCTION.simMeasure = 'Kendall`s tau';
             end
  
        case 3
            readRefG = nk_input('Read in reference group data externally', 0, 'mq', ...
                ['From MATLAB workspace |' ...
                'From file'], [0,1], 0);

            switch readRefG
                case 0
                    GRAPHCONSTRUCTION.refGroup = nk_input('Define reference group dataset (as named in MATLAB workspace)',0,'e',[]);
          
                case 1
                    GRAPHCONSTRUCTION.refGroup = nk_FileSelector(1,0,'Select file with reference group data (in the same format as main data, incl. group and ID column)','.*\.txt$|.*\.csv');
            end
        case 4
                       GRAPHCONSTRUCTION.parcellation = nk_FileSelector(1,0,'Select parcellation atlas', '.*\.nii$|.*\.img$',pwd);
          
    end


else
    [GRAPHCONSTRUCTION,PX] = return_graphconstr(GRAPHCONSTRUCTION,PX, parentstr, true);
    act = 0;
end

if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end
end


function [GRAPHCONSTRUCTION, PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)
global NM
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    menustr = ['Normative network + 1 (e.g., Drenthen et al., 2018)                (Normative network + 1)|', ... 
        'KL divergence of ROI voxel density distributions (Kong et al., 2014)               (KL divergence)' ];
    menuact = {'Normative network + 1', ...
        'KL divergence'};
        
if isfield(GRAPHCONSTRUCTION,'method'), def = find(strcmp(menuact,GRAPHCONSTRUCTION.method)); else, def = 1; end %what does this line do?

nk_PrintLogo
mestr = 'Graph construction methods'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
if isempty(def), def=1; end

act = nk_input(mestr,0, 'mq', menustr, menuact, def);
if ~strcmp(act,'BACK'), GRAPHCONSTRUCTION.method = act; end
if ~exist('PX','var'), PX = []; end

switch act{1}
    case 'KL divergence'
        GRAPHCONSTRUCTION.method = 'KL divergence';
        PX = nk_AddParam([], [], [], PX,'reset');
        if isfield(GRAPHCONSTRUCTION,'KL divergence')
            ParcellationAtlas = GRAPHCONSTRUCTION.parcellation;
            
        else
            ParcellationAtlas = '(undefined)';
        end
        %GRAPHCONSTRUCTION.parcellation = %nk_input('Parcellation atlas (incl. path to file)', 0, 's');
    case 'Normative network + 1'
        GRAPHCONSTRUCTION.method = 'Normative network + 1';
        PX = nk_AddParam([], [], [], PX,'reset');
        
        GRAPHCONSTRUCTION.simMeasure = 'Mutual information';

%         readVarT = nk_input('Read in variable types externally', 0, 'mq', ...
%             ['From MATLAB workspace |' ...
%             'From file |'], [0,1], 0);
% 
%         switch readVarT
%             case 0
%                 GRAPHCONSTRUCTION.variableTypes = nk_input('Define variable types vector',0,'e',[],[1 (size(NM.Y{1},2)-2)]);
%             case 1
%                 GRAPHCONSTRUCTION.variableTypes = nk_input('Define path to file with variable types vector',0,'s');
%         end
end
else
    GRAPHCONSTRUCTION.method = 'Normative network + 1'; 
    GRAPHCONSTRUCTION.refGroup = 'undefined'; 
    GRAPHCONSTRUCTION.simMeasure = 'Mutual information';
    PX = nk_AddParam([], [], [], PX,'reset');
end
% if ~strcmp(MethodOption,'BACK'), GRAPHCONSTRUCTION.method = MethodOption; end
% if ~exist('PX','var'), PX = []; end
% 
% else
%     GRAPHCONSTRUCTION.method = 'KL divergence';
% 
end

% function [GRAPHCONSTRUCTION, PX] = return_graphSimMeasure(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)
% global NM
% 
% if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
% 
% if ~defaultsfl
%     menustr = ['Mututal information                (Mutual information)|', ...
%         'Pearson correlation (check assumptions)   (Pearson correlation)|', ...
%         'Spearson correlation                      (Spearson correlation)']; % 'KL divergence (Kong et al., 2014)               (KL divergence)|' ...
%     menuact = {'Mutual information', ...
%         'Pearson correlation', ...
%         'Spearson correlation'}; %'KL divergence', ...
% 
%     if isfield(GRAPHCONSTRUCTION,'simMeasure'), def = find(strcmp(menuact,GRAPHCONSTRUCTION.simMeasure)); else, def = 1; end %what does this line do?
% 
%     nk_PrintLogo
%     mestr = 'Measure to quantify similarity between variables'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr);
%     if isempty(def), def=1; end
% 
%     act = nk_input(mestr,0, 'mq', menustr, menuact, def);
%     if ~strcmp(act,'BACK'), GRAPHCONSTRUCTION.simMeasure= act; end
%     if ~exist('PX','var'), PX = []; end
%     switch act{1}
%         case 'Mutual information'
%             GRAPHCONSTRUCTION.simMeasure = 'Mutual information';
%             PX = nk_AddParam([], [], [], PX,'reset');
%         
%         case 'Pearson correlation'
%             GRAPHCONSTRUCTION.simMeasure = 'Pearson correlation';
%             PX = nk_AddParam([], [], [], PX,'reset');
%         case ' Spearman correlation'
%             GRAPHCONSTRUCTION.simMeasure = 'Spearman correlation';
%             PX = nk_AddParam([], [], [], PX,'reset');
%     end
% else
%     GRAPHCONSTRUCTION.simMeasure = 'Mutual information'; PX = nk_AddParam([], [], [], PX,'reset');
% end
% end


