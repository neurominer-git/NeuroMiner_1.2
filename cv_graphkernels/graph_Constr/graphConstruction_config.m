function [GRAPHCONSTRUCTION, PX, act ] = graphConstruction_config(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)
global NM
% what method should be used for the construction of the networks?
ConstructionMethod = 'KL divergence';
ParcellationAtlas = 'Hammers.nii';
VariableTypesVec = '';
ReferenceGroup = '';

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    if isempty(GRAPHCONSTRUCTION), [GRAPHCONSTRUCTION, PX] = graphConstruction_config(GRAPHCONSTRUCTION, PX, parentstr, true); end
    if isfield(GRAPHCONSTRUCTION,'method'), ConstructionMethod = GRAPHCONSTRUCTION.method; end
    if isfield(GRAPHCONSTRUCTION,'parcellation'), ParcellationAtlas = GRAPHCONSTRUCTION.parcellation; end
    if isfield(GRAPHCONSTRUCTION, 'variableTypes'), VariableTypesVec = GRAPHCONSTRUCTION.variableTypes; end
    if isfield(GRAPHCONSTRUCTION, 'refGroup'), VariableTypesVec = GRAPHCONSTRUCTION.refGroup; end
    % BETTER, as some of the options are only relevant for some
    % construction methods, to open an additional menu after method is
    % selected
    if ~ischar(ConstructionMethod)
        GCSTR_METHOD = 'undefined';
    else
        GCSTR_METHOD = ConstructionMethod;
    end

    if strcmp(GCSTR_METHOD,'KL divergence')
        if isempty(ParcellationAtlas)
            GCSTR_PARC = '(undefined)';
        else
            GCSTR_PARC = ParcellationAtlas;
        end
        menustr = ['Select graph construction method [ ' GCSTR_METHOD ' ]|', ...
            'Define parcellation atlas [ ' GCSTR_PARC ' ]'];
        menuact = 1:2;
    end

    if strcmp(GCSTR_METHOD, 'Group deviation')
        if isempty(VariableTypesVec)
            GCSTR_VARTVEC = '(undefined)';
        else
            GCSTR_VARTVEC = sprintf('vector length: %d', length(VariableTypesVec));
        end

        if isempty(ReferenceGroup)
            GCSTR_REFG = '(undefined)';
        else
            GCSTR_REFG = sprintf('reference group: %s', ReferenceGroup);
        end
        menustr = ['Select graph construction method [ ' GCSTR_METHOD ' ]|', ...
            'Define variable types (vector) [ ' GCSTR_VARTVEC ']|', ...
            'Define reference group [' GCSTR_REFG ']'];
        menuact = 1:3;
    end

    %[menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, CALIBUSE);


    nk_PrintLogo
    mestr = 'Graph construction'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);

    switch act
        case 1
            [GRAPHCONSTRUCTION,PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr);
            %ConstructionMethod = 'KL divergence';
        case 2
            switch GCSTR_METHOD
                case 'KL divergence'
                    readParcA = nk_input('Read in parcellation atlas externally', 0, 'mq', ...
                        ['From MATLAB workspace |' ...
                        'From file |'], [0,1], 0);

                    switch readParcA
                        case 0
                            GRAPHCONSTRUCTION.parcellation = nk_input('Parcellation atlas variable name in Matlab workspace',0,'e');
                        case 1
                            GRAPHCONSTRUCTION.parcellation = nk_FileSelector(1,0,'Select parcellation atlas', '.*\.nii$|.*\.img$',pwd);
                    end
                               
                case 'Group deviation'
                    readVarT = nk_input('Read in variable types externally', 0, 'mq', ...
                        ['From MATLAB workspace |' ...
                        'From file |'], [0,1], 0);

                    switch readVarT
                        case 0
                            GRAPHCONSTRUCTION.variableTypes = nk_input('Define variable types vector',0,'e',[],[1 (size(NM.Y{1},2)-2)]);
                        case 1
                            GRAPHCONSTRUCTION.variableTypes = nk_FileSelector(1,0,'Select file containing variable types vector','.*\.txt$|.*\.csv');
                    end
                   
                    %GRAPHCONSTRUCTION.refGroup = nk_input('Reference group name (as defined in dataset)', 0, 's');
        
            end
        case 3
            readRefG = nk_input('Read in reference group data externally', 0, 'mq', ...
                ['From MATLAB workspace |' ...
                'From file |'], [0,1], 0);

            switch readRefG
                case 0
                    GRAPHCONSTRUCTION.refGroup = nk_input('Define reference group dataset (as named in MATLAB workspace)',0,'e',[],[1 (size(NM.Y{1},2)-2)]);
                case 1
                    GRAPHCONSTRUCTION.refGroup = nk_FileSelector(1,0,'Select file with reference group data (in the same format as main data, incl. group and ID column)','.*\.txt$|.*\.csv');
            end

             
            %ParcellationAtlas = 'Hammers.nii';
    end
else
    [GRAPHCONSTRUCTION,PX] = return_graphconstr(GRAPHCONSTRUCTION,PX, parentstr, true);
    act = 0;
end
%GRAPHCONSTRUCTION.method = ConstructionMethod;
%GRAPHCONSTRUCTION.parcellation = ParcellationAtlas;

if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end
end

%     
%     switch act
%         case 1
%             [GRAPHCONSTRUCTION, PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr);
%             %GRAPHCONSTRUCTION.parcellation = cv_ExtGraphConstr_config(GRAPHCONSTRUCTION.method, [], [], true);
%         case 2 % for compatibility
%             %t_act = 1; while t_act > 0, [ParcellationAtlas, t_act] = cv_ExtGraphConstr_config(ConstructionMethod, ParcellationAtlas, 0, navistr);end
%             GRAPHCONSTRUCTION.parcellation = nk_input('Define parcellation atlas file path',0, 's', ParcellationAtlas);
%                 ;
%         case 1000
%             CALIBUSE = nk_AskCalibUse_config(mestr, CALIBUSE);
%             
%     end
% else
%     [GRAPHCONSTRUCTION, PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr, true);
%     %GRAPHCONSTRUCTION.parcellation = cv_ExtGraphConstr_config(GRAPHCONSTRUCTION.method, [], [], true);
%     act = 0;
% end

% 
function [GRAPHCONSTRUCTION, PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)
global NM
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    menustr = ['KL divergence (Kong et al., 2014)               (KL divergence)|' ...
        'Deviation from group (Das et al., 2018)                (Group deviation)|'];
    menuact = {'KL divergence', ...
        'Group deviation'};

    
if isfield(GRAPHCONSTRUCTION,'method'), def = find(strcmp(menuact,GRAPHCONSTRUCTION.method)); else, def = 1; end %what does this line do?

nk_PrintLogo
mestr = 'Graph construction methods'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
if isempty(def), def=1; end

act = nk_input(mestr,0, 'mq', menustr, menuact, def);
if ~strcmp(act,'BACK'), GRAPHCONSTRUCTION.method = act; end
if ~exist('PX','var'), PX = []; end

switch act{1}
    case 'KL divergence'
        GRAPHCONSTRUCTION.method = 'KL divergence';
        PX = nk_AddParam([], [], [], PX,'reset');
%         if isfield(GRAPHCONSTRUCTION,'KL divergence')
%             ParcellationAtlas = GRAPHCONSTRUCTION.parcellation;
%             
%         else
%             ParcellationAtlas = 'Hammers.nii';
%         end
        GRAPHCONSTRUCTION.parcellation = nk_input('Parcellation atlas (incl. path to file)', 0, 's');
    case 'Group deviation'
        GRAPHCONSTRUCTION.method = 'Group deviation';
        PX = nk_AddParam([], [], [], PX,'reset');
        readVarT = nk_input('Read in variable types externally', 0, 'mq', ...
            ['From MATLAB workspace |' ...
            'From file |'], [0,1], 0);

        switch readVarT
            case 0
                GRAPHCONSTRUCTION.variableTypes = nk_input('Define variable types vector',0,'e',[],[1 (size(NM.Y{1},2)-2)]);
            case 1
                GRAPHCONSTRUCTION.variableTypes = nk_input('Define path to file with variable types vector',0,'s');
        end
end
else
    GRAPHCONSTRUCTION.method = 'Group deviation'; GRAPHCONSTRUCTION.variableTypes = 9999; PX = nk_AddParam([], [], [], PX,'reset');
end
% if ~strcmp(MethodOption,'BACK'), GRAPHCONSTRUCTION.method = MethodOption; end
% if ~exist('PX','var'), PX = []; end
% 
% else
%     GRAPHCONSTRUCTION.method = 'KL divergence';
% 



end


