function [ParcellationAtlas, act] = cv_ExtGraphConstr_config(ConstructionMethod, ParcellationAtlas, defaultsfl, parentstr)
global NM
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if (~exist('ParcellationAtlas','var') || isempty(ParcellationAtlas)); 
         ParcellationAtlas = cv_ExtGraphConstr_config(ConstructionMethod, [], [], 1);
    end
    mn_str = []; ParcellationAtlasStr = {'Hammers','Schaefer100'}; mn_act=[];
    switch ConstructionMethod
        case {'KL divergence'}
            mn_str = [mn_str sprintf('Define atlas for %s [ %s ]|',ConstructionMethod, ParcellationAtlasStr{ParcellationAtlas})]; mn_act = 1;
        case {'alternative'}
            ParcellationAtlas = 1; L = NM.label; L(isnan(L))=[]; act=0; return
        otherwise
            ParcellationAtlas = 1;
    end
    
    switch act
        
        case 1
        
            ParcellationAtlas = nk_input(sprintf('Define %s decomposition',ConstructionMethod),0,'m', ...
                            ['Absolute number range [ 1 ... n ] of eigenvectors|' ...
                             'Percentage range [ 0 ... 1 ] of max dimensionality|' ...
                             'Energy range [ 0 ... 1 ]of maximum decomposition'],1:3, ParcellationAtlas);   
            switch ParcellationAtlas 
                case {2,3}
                    dims = 0.8;
                case 1
                    dims = floor(size(NM.Y{NM.TrainParam.FUSION.M(1)},2)/10);
            end
        case 2
            switch ParcellationAtlas
                case 1
                    inpstr = 'Dimensionalities to project data on (e.g: 1 5 10 or Start:Step:Stop)';
                case 2
                    inpstr = 'Percentages of max dimensionality as defined by training sample size (e.g: 25 50 75 or Start:Step:Stop)';
                case 3
                    inpstr = 'PCA ratios of complete decomposition (e.g: 0.25 0.5 0.75 or Start:Step:Stop)';
            end
            dims = nk_input(inpstr, 0, 'e', defdims); 
    end
     if numel(mn_act)<2, act = 0; end
else
    switch ConstructionMethod
        case {'PCA', 't-SNE'}
            dims = 0.8; ParcellationAtlas = 3;
        case 'SparsePCA'
            dims = floor(numel(NM.cases)/10); ParcellationAtlas = 1;
        case 'PLS'
            dims = 1;
        case {'LDA','GDA'}
            L = NM.label; L(isnan(L))=[]; dims = numel(unique(L)); 
        otherwise
            dims = floor(size(NM.Y{NM.TrainParam.FUSION.M(1)},2)/10); ParcellationAtlas = 1;
    end
    act = 0;
end

end