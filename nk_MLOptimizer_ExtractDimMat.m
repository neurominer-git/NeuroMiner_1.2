function mapYi = nk_MLOptimizer_ExtractDimMat(mapY, dim_index, cPs)
global PREPROC STACKING MULTILABEL RAND

[m,n,o]  = size(mapY.Tr);
mapYi    = mapY;
mapYi.Tr = cell(m,n,o);
mapYi.CV = cell(m,n,o);
mapYi.Ts = cell(m,n,o);
IMPUTE.flag = 0;
if iscell(PREPROC), iPREPROC = PREPROC{1}; else, iPREPROC = PREPROC; end
BINMOD = iPREPROC.BINMOD;
if isfield(RAND,'Decompose') && RAND.Decompose == 2
    BINMOD = 0;
end
if isfield(iPREPROC,'LABELMOD') && isfield(iPREPROC.LABELMOD,'LABELIMPUTE') && ~strcmp(iPREPROC.LABELMOD.LABELIMPUTE.method,'none') 
    IMPUTE = iPREPROC.LABELMOD.LABELIMPUTE; IMPUTE.flag = true; 
end
nclass = numel(mapY.TrL{1,1,1});
% Check if each label in MULTILABEL mode has it ownn data container
label_dim=1;
if MULTILABEL.curdim>1
    if size(mapY.Tr{1,1,1}{1},2) > 1
        label_dim = MULTILABEL.curdim;
    end
end
for i = 1:m
     for j = 1:n
         for k = 1:o
             if BINMOD || STACKING.flag == 1
                 mapYi.Tr{i,j,k} = cell(1,nclass);
                 mapYi.CV{i,j,k} = cell(1,nclass);
                 mapYi.Ts{i,j,k} = cell(1,nclass);
                 for l = 1:nclass
                    if iscell(mapY.Tr{i,j,k}{l})
                        mapYi.Tr{i,j,k}{l} = mapY.Tr{i,j,k}{l}{dim_index, label_dim};
                        mapYi.CV{i,j,k}{l} = mapY.CV{i,j,k}{l}{dim_index, label_dim};
                        mapYi.Ts{i,j,k}{l} = mapY.Ts{i,j,k}{l}{dim_index, label_dim};
                        if isfield(mapYi,'VI'), mapYi.VI{i,j,k}{l} = mapY.VI{i,j,k}{l}{dim_index, label_dim}; end
                    else
                        mapYi.Tr{i,j,k}{l} = mapY.Tr{i,j,k}{l};
                        mapYi.CV{i,j,k}{l} = mapY.CV{i,j,k}{l};
                        mapYi.Ts{i,j,k}{l} = mapY.Ts{i,j,k}{l};
                        if isfield(mapYi,'VI'), mapYi.VI{i,j,k}{l} = mapY.VI{i,j,k}{l}; end
                    end
                    if size(mapY.TrL,3)>1
                        [ mapYi.TrL{i,j,k}{l}, mapYi.Tr{i,j,k}{l}, mapYi.TrInd{i,j,k}{l} ] = nk_LabelImputer( mapY.TrL{i,j,k}{l}, mapYi.Tr{i,j,k}{l}, mapYi.TrInd{i,j,k}{l}, cPs{l}, IMPUTE);
                    else
                        [ mapYi.TrL{i,j}{l}, mapYi.Tr{i,j,k}{l}, mapYi.TrInd{i,j}{l} ] = nk_LabelImputer( mapY.TrL{i,j}{l}, mapYi.Tr{i,j,k}{l}, mapYi.TrInd{i,j}{l}, cPs{l}, IMPUTE);
                    end
                 end
             else
                if iscell(mapY.Tr{i,j,k})
                    mapYi.Tr{i,j,k} = mapY.Tr{i,j,k}{dim_index, label_dim};
                    mapYi.CV{i,j,k} = mapY.CV{i,j,k}{dim_index, label_dim};
                    mapYi.Ts{i,j,k} = mapY.Ts{i,j,k}{dim_index, label_dim};
                    if isfield(mapYi,'VI'), mapYi.VI{i,j,k} = mapY.VI{i,j,k}{dim_index, label_dim}; end
                else
                    mapYi.Tr{i,j,k} = mapY.Tr{i,j,k};
                    mapYi.CV{i,j,k} = mapY.CV{i,j,k};
                    mapYi.Ts{i,j,k} = mapY.Ts{i,j,k};
                    if isfield(mapYi,'VI'), mapYi.VI{i,j,k} = mapY.VI{i,j,k}; end
                    
                end
                for l = 1:nclass
                    if size(mapY.TrL,3)>1
                        [ mapYi.TrL{i,j,k}{l}, mapYi.Tr{i,j,k}, mapYi.TrInd{i,j,k}{l} ] = nk_LabelImputer( mapY.TrL{i,j,k}{l}, mapYi.Tr{i,j,k}, mapYi.TrInd{i,j,k}{l}, cPs{l}, IMPUTE);
                    else
                        [ mapYi.TrL{i,j}{l}, mapYi.Tr{i,j,k}, mapYi.TrInd{i,j}{l} ] = nk_LabelImputer( mapY.TrL{i,j}{l}, mapYi.Tr{i,j}, mapYi.TrInd{i,j}{l}, cPs{l}, IMPUTE);
                    end
                end
             end
         end
    end
end