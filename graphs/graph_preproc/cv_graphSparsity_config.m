function [ GRAPHSPARSITY, PX, act ] = cv_graphSparsity_config(GRAPHSPARSITY, PX, parentstr, defaultsfl)

SparsityPerc = []; 

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
   
    if isfield(GRAPHSPARSITY,'perc') && ~isempty(GRAPHSPARSITY.perc), SparsityPerc = GRAPHSPARSITY.perc; end

    if ~isempty(SparsityPerc)
        PercDef = 1;
        SPARSITYPERCSTR = ['yes, ' nk_ConcatParamstr(SparsityPerc)];
    else
        PercDef = 2;
        SPARSITYPERCSTR = 'no';
    end
        
    menustr = ['Define sparsity threshold(s) [' SPARSITYPERCSTR ']'];
    menuact = 1;
    
    nk_PrintLogo
    mestr = 'Graph sparsity threshold setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            if PercDef == 1, PercDef = 2; elseif PercDef == 2, PercDef = 1; end
            if PercDef == 1 
                SparsityPerc = nk_input('Define percentage(s) to determine edge weight cutoff',0, 'e', SparsityPerc);
                PX = nk_AddParam(SparsityPerc, 'SparsityPerc', 1, []);
            else
                SparsityPerc = [];
            end
    end
else
    act = 0;
end

GRAPHSPARSITY.perc = SparsityPerc;
% Generate parameter array for preprocessing pipeline runner
if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end

