function [GRAPHCONSTRUCTION, PX, act ] = graphConstruction_config(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)

% what method should be used for the construction of the networks?
ConstructionMethod = 'KL divergence';
ParcellationAtlas = 'Hammers.nii';

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    if isempty(GRAPHCONSTRUCTION), [GRAPHCONSTRUCTION, PX] = graphConstruction_config(GRAPHCONSTRUCTION, PX, parentstr, true); end
    if isfield(GRAPHCONSTRUCTION,'method'), ConstructionMethod = GRAPHCONSTRUCTION.method; end
    if isfield(GRAPHCONSTRUCTION,'parcellation'), ParcellationAtlas = GRAPHCONSTRUCTION.parcellation; end
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
    end

    
    menustr = ['Select graph construction method [ ' GCSTR_METHOD ' ]|', ...
                   'Define dimensionality of mapping results [ ' GCSTR_PARC ' ]'];
    menuact = 1:2;

    nk_PrintLogo
    mestr = 'Graph construction'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);

    switch act
        case 1
            ConstructionMethod = 'KL divergence';
        case 2 
            ParcellationAtlas = 'Hammers.nii';
    end
else
    act = 0;
end
GRAPHCONSTRUCTION.method = ConstructionMethod;
GRAPHCONSTRUCTION.parcellation = ParcellationAtlas;

if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
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
% function [GRAPHCONSTRUCTION, PX] = return_graphconstr(GRAPHCONSTRUCTION, PX, parentstr, defaultsfl)
% global EXPERT NM
% if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end
% 
% if ~defaultsfl
%     menustr = ['KL divergence (Kong et al., 2014)               (KL divergence)|'];
%     menuact = {'KL divergence'};
%     
% if isfield(GRAPHCONSTRUCTION,'method'), def = find(strcmp(menuact,GRAPHCONSTRUCTION.method)); else, def = 1; end
% 
% nk_PrintLogo
% mestr = 'Graph construction methods'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
% if isempty(def), def=1; end
% 
% act = nk_input(mestr,0, 'mq', menustr, menuact, def);
% switch act
%     case 1
%         ConstructionMethod = nk_input(sprintf('Define %s decomposition',RedMode),0,'m', ...
%                              ['KL divergence|' ...
%                              'alternative'],1:2, PercMode); 
%         MethodOption = 'KL divergence';
%     case 2
%         MethodOption = 'BACK';
% end
% if ~strcmp(MethodOption,'BACK'), GRAPHCONSTRUCTION.method = MethodOption; end
% if ~exist('PX','var'), PX = []; end
% 
% else
%     GRAPHCONSTRUCTION.method = 'KL divergence';
% 



end


