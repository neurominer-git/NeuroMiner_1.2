% =========================================================================
% FORMAT param = nk_MLI_config(param)
% =========================================================================
% NeuroMiner menu configurator for the interpretation of model predictions 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2022
function [param, act] = nk_MLI_config(param, defaultsfl)

if ~exist("defaultsfl","var") || isempty(defaultsfl)
    defaultsfl=false;
end

method  = 'posneg';
upper_thresh = 95;
lower_thresh = 5;
nperms  = 1000;
max_iter = 1000;
n_visited = 100;
frac    = .1;
usemap  = 0;
mapfeat = 'cvr';
cutoff = [-2 2];
cutoffmode = 'absolute';
cutoffoperator = 1;
znormdata = 1;

if ~defaultsfl
    
    if isfield(param,"MLI")
        method = param.MLI.method;
        nperms = param.MLI.nperms;
        frac = param.MLI.frac;
        upper_thresh = param.MLI.upper_thresh;
        lower_thresh = param.MLI.lower_thresh;
        max_iter = param.MLI.max_iter;
        n_visited = param.MLI.n_visited;
        usemap = param.MLI.MAP.flag;
        mapfeat = param.MLI.MAP.map;
        cutoff = param.MLI.MAP.cutoff;
        cutoffmode = param.MLI.MAP.percentmode;
        cutoffoperator = param.MLI.MAP.operator;
        znormdata = param.MLI.znormdata;
    end

    switch method
        case 'posneg'
            MethodStr = sprintf('Upper/lower percentiles (%g%%/%g%%)', upper_thresh, lower_thresh);
            OcclusionUpperThreshStr = sprintf('Define upper percentile [ %g ]|', upper_thresh);
            OcclusionLowerThreshStr = sprintf('Define lower percentile [ %g ]|', lower_thresh);
            DEFMETHOD = 1;
            mnuact = [1 2 3];
        case 'median'
            MethodStr = 'Median';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 2;
            mnuact = 1;
        case 'medianflip'
            MethodStr = 'Median span flipped';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 3;
            mnuact = 1;
        case 'random'
            MethodStr = 'Random value';
            OcclusionUpperThreshStr = '';
            OcclusionLowerThreshStr = '';
            DEFMETHOD = 4;
            mnuact = 1;
    end
    OcclusionMethodStr = ['Define occlusion method [ ' MethodStr ' ]|']; 
    
    if isinf(nperms)
        IterStr = 'Automated stopping';
    else
        IterStr = sprintf('%g iterations', nperms);
    end
    OcclusionIterStr = ['Define no. of iterations [ ' IterStr ' ]|']; 
    mnuact = [mnuact 4];
    
    if isinf(nperms)
        OcclusionVisitedStr = ['Define minimum number of feature visits [ ' num2str(n_visited) ' ]|'];
        mnuact = [ mnuact 5 ];
    else
        OcclusionVisitedStr = '';
    end
    OcclusionFracStr = ['Define fraction of features to be visited per iteration [ ' num2str(frac) ' ]|'];
    mnuact = [ mnuact 6 ];
    
    if ~usemap 
        MapFlagStr = 'no map specified'; 
    else
        MapFlagStr = 'yes';
    end
    OcclusionMapFlagStr = ['Use map from model visualization to operate in pre-determined feature space [ ' MapFlagStr ' ]|'];
    mnuact = [ mnuact 7 ];
    if usemap
        switch mapfeat
            case 'cvr'
               MapFeatStr = 'Cross-validation ratio map'; DEFMAP = 1; 
            case 'p_sgn'
               MapFeatStr = 'Sign-based consistency map (p value, uncorrected)'; DEFMAP = 2;
            case 'p_FDR_sgn'
               MapFeatStr = 'Sign-based consistency map (p value, FDR-corrected)'; DEFMAP = 3;   
        end
        OcclusionMapFeatStr = ['Define which map to use [ ' MapFeatStr ' ]|'];
        mnuact = [ mnuact 8 ];
        OcclusionMapCutoffStr = ['Define cutoff value (scalar or 2-value vector) [ ' num2str(cutoff) ' ]|'];
        mnuact = [ mnuact 9 ];
        switch cutoffmode
            case 'absolute'
                DEFCUTOFFMODE = 1;
            case 'percentile'
                DEFCUTOFFMODE = 2;
        end
        OcclusionMapCutoffModeStr = ['Define cutoff method [ ' cutoffmode ' ]|'];
        mnuact = [ mnuact 10 ];
        if numel(cutoff)>1
            MapCutoffOperatorStr = {'<,>', '<=,>=', '>,<', '>=,<='};
        else
            MapCutoffOperatorStr = {'<', '<=', '>', '>=', '==', '~='};
        end
        OcclusionMapCutoffOperator = ['Define cutoff operator [ ' MapCutoffOperatorStr{cutoffoperator} ' ]|'];
        mnuact = [ mnuact 11 ];
    else
        OcclusionMapFeatStr = '';
        OcclusionMapCutoffStr = '';
        OcclusionMapCutoffModeStr = '';
        OcclusionMapCutoffOperator = '';
    end
    
    ZnormDataOpts = {'None','Mean centering','Z-normalization'};
    OcclusionZnormData = ['Produce z-normalized prediction change estimates [ ' ZnormDataOpts{znormdata} ' ]' ];
    mnuact = [ mnuact 12 ];

    mnustr = [OcclusionMethodStr ...
              OcclusionUpperThreshStr ...
              OcclusionLowerThreshStr ...
              OcclusionIterStr ...
              OcclusionVisitedStr ...
              OcclusionFracStr ...
              OcclusionMapFlagStr ...
              OcclusionMapFeatStr ...
              OcclusionMapCutoffStr ...
              OcclusionMapCutoffModeStr ...
              OcclusionMapCutoffOperator ...
              OcclusionZnormData];
    
    nk_PrintLogo
    act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> Model interpreter settings',0,'mq', mnustr, mnuact);
    
    switch act
        case 1
            method = char(nk_input('Define occlusion method', 0, 'm', ...
                ['Replace feature with upper and lower percentile cutoffs|' ...
                 'Replace feature with median|' ...
                 'Replace feature by adding/subtracting the median quantile|' ...
                 'Replace feature by randomly picking a value from the feature''s training sample distribution'], ...
                 {'posneg','median','medianflip','random'}, DEFMETHOD));
        case 2
            upper_thresh = nk_input('Define upper percentile', 0, 'e', upper_thresh);
        case 3
            lower_thresh = nk_input('Define upper percentile', 0, 'e', lower_thresh);
        case 4
            nperms = nk_input('Define no. of iterations [inf for automated mode]', 0, 'e', nperms);
        case 5
            n_visited = nk_input('Define minimum number of visits per feature for automated mode', 0, 'i', n_visited);
        case 6
            frac = nk_input('Define fraction of data space to be visited per iteration', 0, 'e', frac);
        case 7
            if ~usemap, usemap=1; else, usemap = 0; end
        case 8
            mapfeat = char(nk_input('Define occlusion method', 0, 'm', ...
                ['Cross-validation ratio map|' ...
                 'Sign-based consistency map (p value, uncorrected)|' ...
                 'Sign-based consistency map (p value, FDR-corrected)'], ...
                 {'cvr','p_sgn','p_FDR_sgn'}, DEFMAP));
        case 9
            switch mapfeat
                case 'cvr'
                   MapFeatStr = 'Cross-validation ratio map';
                case 'p_sgn'
                   MapFeatStr = 'Sign-based consistency map (p value, uncorrected)';
                case 'p_FDR_sgn'
                   MapFeatStr = 'Sign-based consistency map (p value, FDR-corrected)';
            end
            cutoff = nk_input(['Define cutoff(s) to be applied to ' MapFeatStr ], 0, 'e', cutoff);
        case 10
            cutoffmode = char(nk_input('Define map thresholding mode', 0, 'm', ...
                ['Absolute mode|' ...
                 'Percentile mode' ], ...
                 {'absolute','percentile'}, DEFCUTOFFMODE));
        case 11
            if numel(cutoff) > 1

                cutoffoperator = nk_input('Define map thresholding operator', 0, 'm', ...
                    ['<,>|' ...
                     '<=,>=|' ...
                     '>,<|' ...
                     '>=,<='], ...
                     1:4, cutoffoperator);
            else
                cutoffoperator = nk_input('Define map thresholding operator', 0, 'm', ...
                    ['<|' ...
                     '<=|' ...
                     '>|' ...
                     '>=|' ...
                     '==|', ...
                     '~='], ... 
                     1:6, cutoffoperator);
            end
        case 12
            znormdata = nk_input('Define centering procedure', 0, 'm', 'None|Mean centering|Z-normalization', 1:3, znormdata);
    end
end
param.MLI.method = method;
param.MLI.nperms = nperms;
param.MLI.frac = frac;
param.MLI.upper_thresh = upper_thresh;
param.MLI.lower_thresh = lower_thresh;
param.MLI.max_iter = max_iter;
param.MLI.n_visited = n_visited;
param.MLI.znormdata = znormdata;
param.MLI.MAP.flag = usemap;
param.MLI.MAP.map = mapfeat;
param.MLI.MAP.cutoff = cutoff;
param.MLI.MAP.percentmode = cutoffmode;
param.MLI.MAP.operator = cutoffoperator;

