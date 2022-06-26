function [NM,MLI] = make_default_MLI_struct(defaultfl, NM, nAnalysis)
% might be useful for older NM analyses; as the function adds it at the
% right location; perhaps we can somehow integrate it in the interface of
% NM too

if ~exist('defaultfl','var') || isempty(defaultfl) || defaultfl
    MLI.method  = 'medianflip';
    MLI.upper_thresh = 95;
    MLI.lower_thresh = 5;
    MLI.nperms  = 1000;
    MLI.max_iter = 1000;
    MLI.n_visited = 100;
    MLI.frac    = .1;
    MLI.usemap  = 0;
    MLI.mapfeat = 'cvr';
    MLI.cutoff = [-2 2];
    MLI.cutoffmode = 'absolute';
    MLI.cutoffoperator = 1;
    MLI.znormdata = 1;
    if exist('NM', 'var') && exist(nAnalysis', 'var')
        NM.analysis{nAnalysis}.params.TrainParam.MLI = MLI; 
    end
else 
    sprintf('Define through NM or set defaultfl=1'); 
end