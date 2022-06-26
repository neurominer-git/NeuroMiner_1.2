function [TrainParam, act] = nk_Grid_config(TrainParam, SVM, varind, defaultsfl, parentstr)
global NM

OptimFlag = 1;
if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl=0; end
[~, CompStr] = nk_ReturnEvalOperator(SVM.GridParam);

nP = 1;
if isfield(TrainParam,'PREPROC')
    if ~isempty(varind)
        PX_preML = nk_ReturnParamChain(TrainParam.PREPROC(varind), true);
    else
        PX_preML = nk_ReturnParamChain(TrainParam.PREPROC, true); 
    end
    nP = prod(PX_preML.steps);
end

switch OptimFlag
    
    case 1
        
        % Define defaults
        switch NM.modeflag
            case 'classification'
                Cdefs                   = 2.^(-6:1:4);
            case 'regression'
                Cdefs                   = 2.^(-6:1:0);
        end
        Gdefs                           = 2.^(-15:2:3);
        Epsdefs                         = [ 0.05 0.075 0.1 0.125 0.15 0.2];
        Nudefs                          = [ 0.2 0.5 0.7];
        Toldefs                         = [1e-6 1e-5 1e-4 1e-3 1e-2];
        Tdefs                           = [.2 .5 .7];
        PolyCoefdefs                    = 0;
        PolyDegrdefs                    = 3;
        WLiters                         = [2 4 6 8 10];
        Neurondefs                      = [25 50 75 100];
        Leafdefs                        = logspace(1,2,10);
        % Random Forest
        Treedefs                        = [25 50 75 100 150 200];
        NumDdefs                        = -1; % = "sqrt";
        Critdefs                        = 1; % = "gini";
        MaxDdefs                        = 0; % = "None";
        MinSSplitdefs                   = 2;
        MinSLeafdefs                    = 1;
        MinWeightFractLeafdefs          = 0.0;
        MaxLeafNodesdefs                = 0; % = "None";
        MinImpDecrdefs                  = 0.0;
        Bootstrapdefs                   = 1;
        OobScoredefs                    = 0;
        %NJobsdefs                       = "None";
        %RandomStatedefs                 = "None";
        %WarmStartdefs                   = "False";
        ClassWeightdefs                 = 0; % = "None";
        CcpAlphadefs                    = 0.0;
        MaxSampdefs                     = 0; % = "None";

        Kdefs                           = 7;
        Weightdefs                      = 1;
        CoxCutoffsdefs                  = [20:10:80];
        SEQOPTstepsdefs                 = 10;
        SEQOPTlimsLdefs                 = 50;
        SEQOPTlimsUdefs                 = 50;
        OptRegul.flag                   = 0;
        GridMaxType                     = 1;
        OptRegul.type                   = 1;
        OptRegul.lambda                 = 1;
        OptRegul.big_gamma              = .5;
        OptRegul.RegulTypeComplexity    = nk_RegulFunc_config(1);
        OptRegul.RegulTypeDiversity     = nk_RegulFunc_config(1);
        NodeSelect.mode                 = 1;
        if isfield(NM.TrainParam.SVM.kernel, 'customfunc_nargin')
            if NM.TrainParam.SVM.kernel.customfunc_nargin >0
                for n = 1:NM.TrainParam.SVM.kernel.customfunc_nargin
                    argName = sprintf('customkernel_arg%d', n); 
                    eval(sprintf("%s = 0;", argName)); 
                end
            end
        end
        

        
        switch CompStr
            case 'above'
                 NodeSelect.perc        = 95;
                 NodeSelect.percvec     = 75:5:95;
            case 'below'
                 NodeSelect.perc        = 5;
                 NodeSelect.percvec     = 5:5:25;
        end
        act                             = 0;
        PX                              = [];
        GRD                             = [];

        if ~defaultsfl 
            
            if isfield(TrainParam,'GRD'), GRD = TrainParam.GRD; end
            
            %% Get current values 
            %==============================================================
            if isfield(GRD,'Cparams'),                  Cdefs   = GRD.Cparams; end
            if isfield(GRD,'Epsparams'),                Epsdefs = GRD.Epsparams; end
            if isfield(GRD,'Nuparams'),                 Nudefs  = GRD.Nuparams; end
            if isfield(GRD,'Gparams'),                  Gdefs = GRD.Gparams; end
            if isfield(GRD,'PolyCoefparams'),           PolyCoefdefs = GRD.PolyCoefparams; end
            if isfield(GRD,'PolyDegrparams'),           PolyDegrdefs = GRD.PolyDegrparams; end
            if isfield(GRD,'Neuronparams'),             Neurondefs = GRD.Neuronparams; end
            if isfield(GRD,'Kparams'),                  Kdefs = GRD.Kparams; end
            if isfield(GRD,'Tolparams'),                Toldefs = GRD.Tolparams; end
            if isfield(GRD,'Leafparams'),               Leafdefs = GRD.Leafparams; end
            % random forest
            if isfield(GRD,'Treeparams'),               Treedefs = GRD.Treeparams; end
            if isfield(GRD,'NumDparams'),               NumDdefs = GRD.NumDparams; end
            if isfield(GRD,'RFCritdefsparams'),               Critdefs = GRD.RFCritdefsparams; end
            if isfield(GRD,'MaxDparams'),               MaxDdefs = GRD.MaxDparams; end
            if isfield(GRD,'RFMinSSplitparams'),               MinSSplitdefs = GRD.RFMinSSplitparams; end
            if isfield(GRD,'RFMinSLeafparam'),               MinSLeafdefs = GRD.RFMinSLeafparam; end
            if isfield(GRD,'RFMinWFLeafparams'),               MinWeightFractLeafdefs = GRD.RFMinWFLeafparams; end
            if isfield(GRD,'RFMaxLeafNparams'),               MaxLeafNodesdefs = GRD.RFMaxLeafNparams; end
            if isfield(GRD,'RFMinImpDecrparams'),               MinImpDecrdefs = GRD.RFMinImpDecrparams; end
            if isfield(GRD,'RFBootstrapparams'),               Bootstrapdefs = GRD.RFBootstrapparams; end
            if isfield(GRD,'RFOobScoreparams'),               OobScoredefs = GRD.RFOobScoreparams; end
            if isfield(GRD,'RFClassWeightparams'),               ClassWeightdefs = GRD.RFClassWeightparams; end
            if isfield(GRD,'RFCcpAlphaparams'),               CcpAlphadefs = GRD.RFCcpAlphaparams; end
            if isfield(GRD,'MaxSampparams'),               MaxSampdefs = GRD.MaxSampparams; end


            if isfield(GRD,'Weightparams'),             Weightdefs = GRD.Weightparams; end
            if isfield(GRD,'CutOffparams'),             SEQOPTstepsdefs = GRD.CutOffparams; end
            if isfield(GRD,'LimsLparams'),              SEQOPTlimsLdefs = GRD.LimsLparams; end
            if isfield(GRD,'LimsUparams'),              SEQOPTlimsUdefs = GRD.LimsUparams; end
            if isfield(GRD,'CoxCutoffparams'),          CoxCutoffsdefs = GRD.CoxCutoffparams; end
            if isfield(GRD,'OptRegul'),                 OptRegul = GRD.OptRegul; end
            if isfield(GRD,'NodeSelect'),               NodeSelect = GRD.NodeSelect; end
            if isfield(GRD,'WLiters'),                  WLiters = GRD.WLiters; end
            if isfield(NM.TrainParam.SVM.kernel,'customfunc_nargin')
                if NM.TrainParam.SVM.kernel.customfunc_nargin >0
                    for n = 1:NM.TrainParam.SVM.kernel.customfunc_nargin
                        argName = sprintf('customkernel_arg%d', n); 
                        if isfield(GRD, argName) 
                            eval(sprintf("%s = GRD.%s", argName, argName)); 
                        end
                    end
                end
            end
            %if isfield(GRD, 'CustomKernel'),            CustomKernel = GRD.CustomKernel; end
            
            %============================================================== 
            menustr = []; menuact = []; n_pars = [];
            
            %% Slack setup
            switch SVM.prog
                
                case {'LIBSVM','LIBLIN','SVMLIT','SVMPRF', 'MEXELM'}
                    ctype = nk_GetLIBSVMClassType(SVM);
                    switch ctype
                        case 1
                            Rparstr = 'Nu-SVC parameter(s)'; [Nustr, n_pars(end+1)] = nk_ConcatParamstr(Nudefs); 
                            PX = nk_AddParam(Nudefs, ['ML-' Rparstr], 2, PX);
                            menustr = sprintf('Define %s [ %s ]', Rparstr, Nustr);                                  menuact = [ menuact 3 ]; 
                        otherwise
                            Cparstr = 'Slack/Regularization parameter(s)'; [Cstr, n_pars(end+1)] = nk_ConcatParamstr(Cdefs); 
                            PX = nk_AddParam(Cdefs, ['ML-' Cparstr], 2, PX);
                            menustr = sprintf('Define %s [ %s ]', Cparstr, Cstr);                                   menuact = [ menuact 1 ]; 
                    end
                    rtype = nk_GetLIBSVMRegrType(SVM);
                    switch rtype
                        case 1
                            Rparstr = 'Eps-SVR parameter(s)'; [Epsstr, n_pars(end+1)] = nk_ConcatParamstr(Epsdefs);
                            PX = nk_AddParam(Epsdefs, ['ML-' Rparstr], 2, PX);
                            menustr = sprintf('%s|Define %s [ %s ]', menustr, Rparstr, Epsstr);                     menuact = [ menuact 2 ];
                        case 2
                            Rparstr = 'Nu-SVR parameter(s)'; [Nustr, n_pars(end+1)] = nk_ConcatParamstr(Nudefs);
                            PX = nk_AddParam(Nudefs, ['ML-' Rparstr], 2, PX);
                            menustr = sprintf('%s|Define %s [ %s ]', menustr, Rparstr, Nustr);                      menuact = [ menuact 3 ];
                    end
                    
                case 'IMRELF'
                    Cparstr = 'Lambda parameter(s)'; [Cstr, n_pars(end+1)] = nk_ConcatParamstr( Cdefs );
                    PX = nk_AddParam(Cdefs, ['ML-' Cparstr], 2, PX);
                    menustr = sprintf('Define %s [ %s ]',Cparstr, Cstr);                                            menuact = [ menuact 1 ];
                    
                case 'matLRN'
                    if ~isfield(GRD,'matLearn') || ~isfield(GRD.matLearn,'Params')
                        GRD.matLearn = nk_matLearn_config(SVM.matLRN,SVM.matLRN.learner.framework,3);
                    end
                    if ~isempty(GRD.matLearn.Params)
                        for i=1:numel(GRD.matLearn.Params)
                            parstr = GRD.matLearn.Params(i).name; [~, n_pars(end+1)] = nk_ConcatParamstr(  GRD.matLearn.Params(i).range );
                            PX = nk_AddParam(GRD.matLearn.Params(i).range, ['ML-' parstr], 2, PX);
                        end
                        str = numel(GRD.matLearn.Params);
                        menustr = sprintf('%s|Define matLearn parameters [ %g main parameters defined ]', menustr, str); 
                        menuact = [ menuact 20 ];
                    else
                        return
                    end
                case {'GLMNET','GRDBST','ROBSVM','ELASVM'}
                    if isfield(GRD,SVM.prog) && ~isempty(GRD.(SVM.prog).Params)
                        for i=1:numel(GRD.(SVM.prog).Params)
                            parstr = GRD.(SVM.prog).Params(i).name; [~, n_pars(end+1)] = nk_ConcatParamstr(  GRD.(SVM.prog).Params(i).range );
                            PX = nk_AddParam(GRD.(SVM.prog).Params(i).range, ['ML-' parstr], 2, PX);
                        end
                        str = sprintf('%g parameters defined', numel(GRD.(SVM.prog).Params));
                    else
                        str = 'undefined'; 
                    end
                    menustr = sprintf('%s|Define %s parameters [ %s ]|', menustr, SVM.prog, str); 
                    menuact = [ menuact 21 ];
                
                    
            end
            
            %% Kernel setup
            switch SVM.prog
                
                case 'IMRELF'
                     Gparstr = 'Sigma parameter(s)'; [Gstr, n_pars(end+1)] = nk_ConcatParamstr(Gdefs);
                     PX = nk_AddParam(Gdefs, ['ML-' Gparstr], 2, PX);
                     menustr = sprintf('%s|Define %s [ %s ]', menustr, Gparstr, Gstr);                              menuact = [ menuact 4 ];                     
                case 'matLRN'
                    
                otherwise
                    
                     switch SVM.kernel.kernstr
                        case {' -t 1', ...
                            ' -t 2', ...
                            ' --t 2', ...
                            ' -t 3', ...
                            ' -t 6', ...
                            ' -t 7', ...
                            'poly', ...
                            'polynomial', ...
                            'Polynomial', ...
                            'polyN', ...
                            'hpolyN', ...
                            'rbf', ...
                            'expo', ...
                            'laplace', ...
                            'cauchy', ...
                            'cubic', ...
                            'tps', ...
                            'r', ...
                            'gauss', ...
                            'gaussian'}
                            Gparstr = 'Kernel parameter(s)'; [Gstr, n_pars(end+1)] = nk_ConcatParamstr(Gdefs);
                            PX = nk_AddParam(Gdefs, ['ML-' Gparstr], 2, PX);
                            menustr = sprintf('%s|Define %s [ %s ]', menustr, Gparstr, Gstr);                           menuact = [ menuact 4 ];
                            % Additional polynomial kernel params
                            switch SVM.kernel.kernstr
                                
                                case {' -t 1', 'poly', 'polynomial', 'Polynomial', 'hpolyN', 'polyN'}
                                    Pcparstr = 'Polynomial coefficients'; [Pcstr, n_pars(end+1)] = nk_ConcatParamstr(PolyCoefdefs);
                                    PX = nk_AddParam(PolyCoefdefs, ['ML-' Pcparstr], 2, PX);
                                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Pcparstr, Pcstr);                 menuact = [ menuact 5 ];
                                    Pdparstr = 'Polynomial degree'; [Pdstr, n_pars(end+1)] = nk_ConcatParamstr(PolyDegrdefs);
                                    PX = nk_AddParam(PolyDegrdefs, ['ML-' Pdparstr], 2, PX);
                                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Pdparstr, Pdstr);                 menuact = [ menuact 6 ];
                                    
                                case ' -t 3'
                                    Pcparstr = 'Sigmoid coefficients'; [Pcstr, n_pars(end+1)] = nk_ConcatParamstr(PolyCoefdefs);
                                    PX = nk_AddParam(PolyCoefdefs, ['ML-' Pcparstr], 2, PX);
                                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Pcparstr, Pcstr);                 menuact = [ menuact 5 ];
                            end
                         case ' -t 4' % precomputed kernels
                             if NM.TrainParam.SVM.kernel.kerndef <8
                                WLparstr = 'WL iterations'; [WLstr, n_pars(end+1)] = nk_ConcatParamstr(WLiters);
                                PX = nk_AddParam(WLiters, ['ML-' WLparstr], 2, PX);
                                menustr = sprintf('%s|Define %s [ %s ]', menustr, WLparstr, WLstr);                        menuact = [ menuact 100];
                             elseif NM.TrainParam.SVM.kernel.kerndef == 8
                                 if NM.TrainParam.SVM.kernel.customfunc_nargin > 0
                                    for n = 1:NM.TrainParam.SVM.kernel.customfunc_nargin
                                        argParstr = sprintf('CFuncparstr_%d', n); 
                                        argStr = sprintf('CFuncstr_%d', n);
                                        arg = sprintf('customkernel_arg%d', n);
                                        
                                        eval(sprintf("%s = 'Custom kernel argument %d'", argParstr, n));
                                        eval(sprintf('[%s, n_pars(end+1)] = nk_ConcatParamstr(%s)', argStr, arg));
                                        eval(sprintf("PX = nk_AddParam(%s, ['ML-' %s], 2, PX)", arg, argParstr)); 
                                        eval(sprintf('X = %s; Y = %s', argParstr, argStr)); 
                                        eval(sprintf("menustr = '%s|Define %s [ %s ]'", menustr, X, Y));
                                        
%                                       
                                        menuact = [ menuact 101];
                                    end
                                  
                                 end
                                                                        
                           
                             end
%                              switch SVM.kernel.kerndef
%                                  case 8 % custom kernel function
%                                      CKFnamestr = 'Custom kernel function (input: 2 matrices, output: kernel matrix)'; [CKFstr, n_pars(end+1)] = nk_ConcatParamstr(CustomKernel);
%                                      PX = nk_AddParam(CustomKernel, ['ML-' CKFnamestr], 2, PX);
%                                      menustr = sprintf('%s|Define %s [ %s ]', menustr, CKFnamestr, CKFstr);                        menuact = [ menuact 101];
%                              
%                              end
                            
                      
                     end             
            end
            
            %% Other param setup
            switch SVM.prog
                case {'LIBSVM','LIBLIN'}
                    if SVM.(SVM.prog).Weighting
                        Weightparstr = 'Weighting exponents'; [Weightstr, n_pars(end+1)] = nk_ConcatParamstr(Weightdefs);
                        PX = nk_AddParam(Weightdefs, ['ML-' Weightparstr], 2, PX);
                        menustr = sprintf('%s|Define %s [ %s ]', menustr, Weightparstr, Weightstr);                     menuact = [ menuact 15 ];
                    end
                case 'MEXELM'
                    Neuronparstr = 'Hidden neurons'; [Neuronstr, n_pars(end+1)] = nk_ConcatParamstr(Neurondefs);
                    PX = nk_AddParam(Neurondefs, ['ML-' Neuronparstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Neuronparstr, Neuronstr);                         menuact = [ menuact 12 ];
                case 'kNNMEX'
                    Kparstr = 'Nearest neighbors'; [Kstr, n_pars(end+1)] = nk_ConcatParamstr(Kdefs);
                    PX = nk_AddParam(Kdefs, ['ML-' Kparstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Kparstr, Kstr);                                   menuact = [ menuact 13 ];
                case 'BLOREG'
                    Tparstr = 'Tolerance parameter(s)'; [Tstr, n_pars(end+1)] = nk_ConcatParamstr(Toldefs);
                    PX = nk_AddParam(Toldefs, ['ML-' Tparstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Tparstr, Tstr);                                   menuact = [ menuact 14 ];
                case 'DECTRE'
                    Lfparstr = 'Leafness parameter(s)'; [Lfstr, n_pars(end+1)] = nk_ConcatParamstr( Leafdefs );
                    PX = nk_AddParam(Leafdefs, ['ML-' Lfstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Lfparstr, Lfstr);                                 menuact = [ menuact 16 ];
                case 'RNDFOR'
                    Dtparstr = 'Number of decision trees'; [Dtstr, n_pars(end+1)] = nk_ConcatParamstr( Treedefs );
                    PX = nk_AddParam(Treedefs, ['ML-' Dtstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Dtparstr, Dtstr);                                 menuact = [ menuact 17 ];
                    
                    Dnparstr = 'Maximum number of features'; [Dnstr, n_pars(end+1)] = nk_ConcatParamstr( NumDdefs );
                    PX = nk_AddParam(NumDdefs, ['ML-' Dnstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Dnparstr, Dnstr);                                 menuact = [ menuact 18 ];

          
                    Critparstr = 'Function to measure the quality of a split'; ; [Critstr, n_pars(end+1)] = nk_ConcatParamstr( Critdefs );
                    PX = nk_AddParam(Critdefs, ['ML-' Critstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Critparstr, Critstr);                                 menuact = [ menuact 31 ];
                    
                    MaxDparstr = 'Maximum depth of the tree'; ; [MaxDstr, n_pars(end+1)] = nk_ConcatParamstr( MaxDdefs );
                    PX = nk_AddParam(MaxDdefs, ['ML-' MaxDstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MaxDparstr, MaxDstr);                                 menuact = [ menuact 32 ];
                    
                    MinSSplitparstr = 'Minimum number of samples to split'; ; [MinSSplitstr, n_pars(end+1)] = nk_ConcatParamstr( MinSSplitdefs );
                    PX = nk_AddParam(MinSSplitdefs, ['ML-' MinSSplitstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MinSSplitparstr, MinSSplitstr);                                 menuact = [ menuact 33 ];
                    
                    MinSLeafparstr = 'Minimum number of samples to be at a leaf'; ; [MinSLeafstr, n_pars(end+1)] = nk_ConcatParamstr( MinSLeafdefs );
                    PX = nk_AddParam(MinSLeafdefs, ['ML-' MinSLeafstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MinSLeafparstr, MinSLeafstr);                                 menuact = [ menuact 34 ];
                    
                    MinWeightFractLeafparstr = 'Minimum weighted fraction of the sum total of weights at a leaf'; ; [MinWeightFractLeafstr, n_pars(end+1)] = nk_ConcatParamstr( MinWeightFractLeafdefs );
                    PX = nk_AddParam(MinWeightFractLeafdefs, ['ML-' MinWeightFractLeafstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MinWeightFractLeafparstr, MinWeightFractLeafstr);                                 menuact = [ menuact 35 ];
                    
                    MaxLeafNodesparstr = 'Maximum number of leaf nodes'; ; [MaxLeafNodesstr, n_pars(end+1)] = nk_ConcatParamstr( MaxLeafNodesdefs );
                    PX = nk_AddParam(MaxLeafNodesdefs, ['ML-' MaxLeafNodesstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MaxLeafNodesparstr, MaxLeafNodesstr);                                 menuact = [ menuact 36 ];
                    
                    MinImpDecrparstr = 'Minimum decrease of impurity'; ; [MinImpDecrstr, n_pars(end+1)] = nk_ConcatParamstr( MinImpDecrdefs );
                    PX = nk_AddParam(MinImpDecrdefs, ['ML-' MinImpDecrstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MinImpDecrparstr, MinImpDecrstr);                                 menuact = [ menuact 37 ];
                    
                    Bootstrapparstr = 'Bootstrap samples yes/no'; ; [Bootstrapstr, n_pars(end+1)] = nk_ConcatParamstr( Bootstrapdefs );
                    PX = nk_AddParam(Bootstrapdefs, ['ML-' Bootstrapstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, Bootstrapparstr, Bootstrapstr);                                 menuact = [ menuact 38 ];
                    
                    OobScoreparstr = 'Out-of-bag samples yes/no'; ; [OobScorestr, n_pars(end+1)] = nk_ConcatParamstr( OobScoredefs );
                    PX = nk_AddParam(OobScoredefs, ['ML-' OobScorestr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, OobScoreparstr, OobScorestr);                                 menuact = [ menuact 39 ];
                    
%                     NJobsparstr = 'Number of jobs to run in parallel'; ; [NJobsstr, n_pars(end+1)] = nk_ConcatParamstr( NJobsdefs );
%                     PX = nk_AddParam(NJobsdefs, ['ML-' NJobsstr], 2, PX);
%                     menustr = sprintf('%s|Define %s [ %s ]', menustr, NJobsparstr, NJobsstr);                                 menuact = [ menuact 40 ];
%                     
%                     RandomStateparstr = 'Randomness of bootstrapping'; ; [RandomStatestr, n_pars(end+1)] = nk_ConcatParamstr( RandomStatedefs );
%                     PX = nk_AddParam(RandomStatedefs, ['ML-' RandomStatestr], 2, PX);
%                     menustr = sprintf('%s|Define %s [ %s ]', menustr, RandomStateparstr, RandomStatestr);                                 menuact = [ menuact 41 ];
%                     
%                     WarmStartparstr = ''; ; [WarmStartstr, n_pars(end+1)] = nk_ConcatParamstr( WarmStartdefs );
%                     PX = nk_AddParam(WarmStartdefs, ['ML-' WarmStartstr], 2, PX);
%                     menustr = sprintf('%s|Define %s [ %s ]', menustr, WarmStartparstr, WarmStartstr);                                 menuact = [ menuact 17 ];
%                     
                    ClassWeightparstr = 'Class weights'; ; [ClassWeightstr, n_pars(end+1)] = nk_ConcatParamstr( ClassWeightdefs );
                    PX = nk_AddParam(ClassWeightdefs, ['ML-' ClassWeightstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, ClassWeightparstr, ClassWeightstr);                                 menuact = [ menuact 42 ];
                    
                    CcpAlphaparstr = 'Complexity parameter for Minimal Cost-Complexity Pruning'; ; [CcpAlphastr, n_pars(end+1)] = nk_ConcatParamstr( CcpAlphadefs );
                    PX = nk_AddParam(CcpAlphadefs, ['ML-' CcpAlphastr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, CcpAlphaparstr, CcpAlphastr);                                 menuact = [ menuact 43 ];
                    
                    MaxSampparstr = 'Number of samples to draw from X to train base estimators (if bootstrap)'; ; [MaxSampstr, n_pars(end+1)] = nk_ConcatParamstr( MaxSampdefs );
                    PX = nk_AddParam(MaxSampdefs, ['ML-' MaxSampstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, MaxSampparstr, MaxSampstr);                                 menuact = [ menuact 44 ];
                    
        
        
        
        
        
   
        ClassWeightdefs                 = "None";
        CcpAlphadefs                    = 0.0;
        MaxSampdefs                     = "None";

                case 'SEQOPT'
                    CutOffparstr = 'No. of threshold for ambiguous case propagation'; 
                    [CutOffstr, n_pars(end+1)] = nk_ConcatParamstr( SEQOPTstepsdefs );
                    PX = nk_AddParam(SEQOPTstepsdefs, ['ML-' CutOffstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, CutOffparstr, CutOffstr);                         menuact = [ menuact 22 ];
                    LimsLparstr = 'Lower population percentage(s) (- from anchor) for ambiguous case propagation'; 
                    [LimsLstr, n_pars(end+1)] = nk_ConcatParamstr( SEQOPTlimsLdefs );
                    PX = nk_AddParam(SEQOPTlimsLdefs, ['ML-' LimsLstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, LimsLparstr, LimsLstr);                           menuact = [ menuact 23 ];
                    LimsUparstr = 'Upper population percentage(s) (+ from anchor) for ambiguous case propagation'; 
                    [LimsUstr, n_pars(end+1)] = nk_ConcatParamstr( SEQOPTlimsUdefs );
                    PX = nk_AddParam(SEQOPTlimsUdefs, ['ML-' LimsUstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, LimsUparstr, LimsUstr);                           menuact = [ menuact 24 ];
                case 'WBLCOX'
                    CoxCutOffparstr = 'Percentile thresholds for class assignment';
                    [CutOffstr, n_pars(end+1)] = nk_ConcatParamstr( CoxCutoffsdefs );
                    PX = nk_AddParam(CoxCutoffsdefs, ['ML-' CutOffstr], 2, PX);
                    menustr = sprintf('%s|Define %s [ %s ]', menustr, CoxCutOffparstr, CutOffstr);                      menuact = [ menuact 25 ];  
            end
            
            GRD.GridMaxType = GridMaxType;
            if prod(n_pars) * nP > 1 
                if OptRegul.flag, regstr = 'Yes'; else, regstr = 'No'; end
                menustr = sprintf('%s|Enable regularization of model selection [ %s ]', menustr, regstr);               menuact = [ menuact 7 ];
                if OptRegul.flag
                    if isfield(TrainParam,'RFE') && ...
                            (isfield(TrainParam.RFE,'Filter') && ...
                            isfield(TrainParam.RFE.Filter,'EnsembleStrategy') && ...
                            isfield(TrainParam.RFE.Filter.EnsembleStrategy,'type') && ...
                            TrainParam.RFE.Filter.EnsembleStrategy.type ~= 9 || ...
                        isfield(TrainParam.RFE,'Wrapper') && ...
                            isfield(TrainParam.RFE.Wrapper,'EnsembleStrategy') && ...
                            isfield(TrainParam.RFE.Wrapper.EnsembleStrategy,'type') && ...
                            TrainParam.RFE.Wrapper.EnsembleStrategy.type ~= 9)
                        switch OptRegul.flag
                            case 1
                                crossstr = 'Model complexity';
                            case 2
                                crossstr = 'Ensemble diversity';
                            case 3
                                crossstr = 'Mixed Criterion (Model Complexity & Ensemble diversity)';
                        end
                        menustr = sprintf('%s|Criterion for cross-parameter model selection [ %s ]', menustr, crossstr);  menuact = [ menuact 8 ];
                    else
                        GRD.OptRegul.type = OptRegul.type;
                    end
                    menustr = sprintf('%s|Define weight (lambda) of SV ratio [ %g ]', menustr, OptRegul.lambda);           menuact = [ menuact 9 ];
                    menustr = sprintf('%s|Define non-linearity (big gamma) of SV ratio [ %g ]', menustr, OptRegul.big_gamma);  menuact = [ menuact 10 ];

                end
                switch NodeSelect.mode
                    case 1
                        modelselstr = 'Single optimum model (no ensemble)';
                    case 2
                        modelselstr = sprintf('Aggregated cross-parameter ensemble %s %g-percentile', CompStr, NodeSelect.perc);
                    case 3
                        modelselstr = sprintf('Optimized cross-parameter ensemble %s %g-percentile', CompStr, NodeSelect.perc);
                    case 4
                        modelselstr = sprintf('Variable-threshold cross-parameter ensemble');
                end
                menustr = sprintf('%s|Specify cross-parameter model selection process [ %s ]', menustr, modelselstr);     menuact = [ menuact 11 ]; 
            end
            % =============================================================
            nk_PrintLogo
            if prod(n_pars)>3, fprintf('\n');fprintf('-> No. of ML parameter combinations: %g ', prod(n_pars)); end
            fprintf('\n\n'); mestr = 'Model optimization parameters'; navistr = [parentstr ' >>> ' mestr]; fprintf('You are here: %s >>> ', parentstr); 
            act = nk_input(mestr, 0,'mq', menustr, menuact );
            
            switch act
                case 1
                    Cdefs         = nk_input([Cparstr ' range'],0,'e',Cdefs);                           PX = nk_AddParam(Cdefs, ['ML-' Cparstr], 2, PX);
                case 2
                    Epsdefs       = nk_input([Rparstr ' range'],0,'e',Epsdefs);                         PX = nk_AddParam(Epsdefs, ['ML-' Rparstr], 2, PX);
                case 3
                    Nudefs        = nk_input([Rparstr ' range [0 ... 1]'],0,'e',Nudefs);                PX = nk_AddParam(Nudefs, ['ML-' Rparstr], 2, PX);
                case 4
                    Gdefs         = nk_input([Gparstr ' range'],0,'e',Gdefs);                           PX = nk_AddParam(Gdefs, ['ML-' Gparstr], 2, PX);
                case 5
                    PolyCoefdefs  = nk_input([Pcparstr ' for polynomial kernels'],0,'e',PolyCoefdefs);  PX = nk_AddParam(PolyCoefdefs, ['ML-' Pcparstr], 2, PX);
                case 6
                    PolyDegrdefs  = nk_input([Pdparstr ' for polynomial kernels'],0,'e',PolyDegrdefs);  PX = nk_AddParam(PolyDegrdefs, ['ML-' Pdparstr], 2, PX);
                case 7
                    if ~OptRegul.flag, OptRegul.flag = 2; end
                    OptRegul.flag   = nk_input('Enable model selection across CV1 parameters ?', 0, ...
                                                    'yes|no',[1,0], OptRegul.flag);
                case 8
                    OptRegul.type   = nk_input('Regularize with',0,'mq', ...
                                              ['Model complexity|' ...
                                               'Ensemble diversity|' ...
                                               'Mixed Criterion (Model Complexity & Ensemble diversity)'],1:3,OptRegul.type);
                    switch OptRegul.type 
                        case 1
                            OptRegul.RegulTypeComplexity = nk_RegulFunc_config;
                        case 2
                            OptRegul.RegulTypeDiversity = nk_RegulFunc_config;
                        case 3
                            OptRegul.RegulTypeComplexity = nk_RegulFunc_config;
                            OptRegul.RegulTypeDiversity = nk_RegulFunc_config;
                    end
                case 9
                    OptRegul.lambda = nk_input('Lambda (weight of SV ratio in score => 1 = equal weight): ', ...
                                               0, 'e', OptRegul.lambda, 1);
                case 10                
                    OptRegul.big_gamma = nk_input('Big gamma (non-linearity for SV ratio (.5 = sqrt, 1 = lin, 2 = square): ', ...
                                               0, 'e', OptRegul.big_gamma, 1);
                case 11
                    Nodemode = nk_input('Specify optimum selection mode', 0, 'mq', ...
                                               ['Select a single optimum parameter node|' ...
                                                'Generate cross-node ensemble by aggregating base learners above predefined percentile|' ... %'Generate cross-node ensemble by applying ensemble strategy on base learners above predefined threshold|' ...
                                                'Automatically determine optimal percentile for optimum cross-node ensemble performance'], [1,2,4], NodeSelect.mode);
                    if Nodemode, NodeSelect.mode = Nodemode; end
                    switch NodeSelect.mode
                        case 2
                            NodeSelect.perc = nk_input(sprintf('Models will be selected %s X%% threshold',CompStr), 0, 'e', NodeSelect.perc, 1);
                        case 3 % This option is not implemented yet in in nk_ModelNodeSelector!!!
%                             NodeSelect.perc = nk_input('Define percentage cutoff (>=%)', 0, 'e', NodeSelect.perc, 1);
%                             if ~isfield(NodeSelect,'EnsembleStrategy'), 
%                                 NodeSelect = nk_CostFun_config(NodeSelect, NM, 1);                      % Default target population
%                                 NodeSelect.SubSpaceFlag = 0;
%                                 NodeSelect = nk_SubSpaceStrategy_config(NodeSelect, NM, 1);
%                                 NodeSelect = nk_EnsembleStrategy2_config(NodeSelect, NM, 1); 
%                             end
%                             NodeSelect = nk_EnsembleStrategy2_config(NodeSelect, NM, [], navistr);
                        case 4
                            NodeSelect.percvec = nk_input('Define [lower : stepping : upper ] percentile search vector ', 0, 'e', NodeSelect.percvec);
                    end
                case 12
                    Neurondefs    = nk_input([Neuronparstr ' range'],0,'e',Neurondefs);                 PX = nk_AddParam(Neurondefs, ['ML-' Neuronparstr], 2, PX);
                case 13
                    Kdefs         = nk_input([Kparstr ' range'],0,'e',Kdefs);                           PX = nk_AddParam(Kdefs, ['ML-' Kparstr], 2, PX);
                case 14
                    Toldefs       =  nk_input([Tparstr ' range'],0,'e',Toldefs);                        PX = nk_AddParam(Toldefs, ['ML-' Tparstr], 2, PX);
                case 15
                    Weightdefs    =  nk_input([Weightparstr ' range'],0,'e',Weightdefs);                PX = nk_AddParam(Weightdefs, ['ML-' Weightparstr], 2, PX);
                case 16
                    Leafdefs    =  nk_input([Lfparstr ' range'],0,'e',Leafdefs);                        PX = nk_AddParam(Leafdefs, ['ML-' Lfparstr], 2, PX);
                case 17
                    Treedefs    =  nk_input([Dtparstr ' range'],0,'e',Treedefs);                        PX = nk_AddParam(Treedefs, ['ML-' Dtparstr], 2, PX);
                case 18
                    NumDdefs    = nk_input('Specify maximum number of features to consider', 0, 'mq', ...
                                               ['Square root (default)|' ...
                                                'Logarithm (log2)|' ... %'Generate cross-node ensemble by applying ensemble strategy on base learners above predefined threshold|' ...
                                                'Fraction|' ...
                                                'Absolute N|' ...
                                                'N features'], [-1,-2, -3, -3, 0], NumDdefs);
                   
                    switch NumDdefs
                        case -3
                            NumDdefs    =  nk_input([Dnparstr ' range'],0,'e',NumDdefs); 
                        case 0
                            
                    end
                            PX = nk_AddParam(NumDdefs, ['ML-' Dnparstr], 2, PX);
                case 31
                    Critdefs    = nk_input('Specify function to measure quality of split', 0, 'mq', ...
                                               ['Gini impurity (default)|' ...
                                                'Log loss|' ...
                                                'Entropy'], [1,2,3], Critdefs);
                    
                    %PX = nk_AddParam(Critdefs, ['ML-' Critparstr], 2, PX);
                case 32
                    MaxDdefs    = nk_input([MaxDparstr ' range'], 'e', MaxDdefs);
                    PX = nk_AddParam(MaxDdefs, ['ML-' MaxDparstr], 2, PX);
                case 33
                    MinSSplitdefs    = nk_input([MinSSplitparstr ' range'], 'e', MinSSplitdefs);
                    PX = nk_AddParam(MinSSplitdefs, ['ML-' MinSSplitparstr], 2, PX);
                case 34
                    MinSLeafdefs    = nk_input([MinSLeafparstr ' range'], 'e', MinSLeafdefs);
                    PX = nk_AddParam(MinSLeafdefs, ['ML-' MinSLeafparstr], 2, PX);
               
                case 35
                    MinWeightFractLeafdefs    = nk_input([MinWeightFractLeafparstr ' range'], 'e', MinWeightFractLeafdefs);
                    PX = nk_AddParam(MinWeightFractLeafdefs, ['ML-' MinWeightFractLeafparstr], 2, PX);
               
                case 36
                    MaxLeafNodesdefs    = nk_input([MaxLeafNodesparstr ' range'], 'e', MaxLeafNodesdefs);
                    PX = nk_AddParam(MaxLeafNodesdefs, ['ML-' MaxLeafNodesparstr], 2, PX);
               
                case 37
                    MinImpDecrdefs    = nk_input([MinImpDecrparstr ' range'], 'e', MinImpDecrdefs);
                    PX = nk_AddParam(MinImpDecrdefs, ['ML-' MinImpDecrparstr], 2, PX);
               
                case 38
                    Bootstrapdefs    = nk_input('Bootstrap yes or no', 0, 'mq', ...
                                               ['Yes|' ...
                                                'No'], [1,0], Bootstrapdefs);
                    PX = nk_AddParam(Bootstrapdefs, ['ML-' Bootstrapparstr], 2, PX);
               
                    
                case 39
                    OobScoredefs    = nk_input('Out-of-bag yes or no', 0, 'mq', ...
                                               ['Yes|' ...
                                                'No'], [1,0], OobScoredefs);
                    PX = nk_AddParam(OobScoredefs, ['ML-' OobScoreparstr], 2, PX);
               
%                 case 40
%                     NJobsdefs    = nk_input([MinSSplitparstr ' range'], 'e', NJobsdefs);
%                     PX = nk_AddParam(NJobsdefs, ['ML-' NJobsparstr], 2, PX);
               
%                 case 41
%                     RandomStatedefs    = nk_input([RandomStateparstr ' range'], 'e', RandomStatedefs);
%                     PX = nk_AddParam(NJobsdefs, ['ML-' NJobsparstr], 2, PX);
                case 42
                    ClassWeightdefs    = nk_input('Class weights', 0, 'mq', ...
                                               ['All equal (default)|' ...
                                                'Balanced|' ...
                                                'Balanced subsample|'...
                                                'Define yourself (will be added in the future)'], [0,1,2,3], ClassWeightdefs);
                    PX = nk_AddParam(ClassWeightdefs, ['ML-' ClassWeightparstr], 2, PX);
                case 43
                    CcpAlphadefs    = nk_input([CcpAlphaparstr ' range'], 'e', CcpAlphadefs);
                    PX = nk_AddParam(CcpAlphadefs, ['ML-' CcpAlphaparstr], 2, PX);
                case 44
                    MaxSampdefs    = nk_input([MMaxSampparstr ' range'], 'e', MaxSampdefs);
                    PX = nk_AddParam(MaxSampdefs, ['ML-' MaxSampparstr], 2, PX);

                case 20
                    GRD.matLearn = nk_matLearn_config(GRD.matLearn,SVM.matLRN.learner.framework,2);
                    if isfield(GRD.matLearn,'Params')
                        for j=1:numel(GRD.matLearn.Params)
                            PX = nk_AddParam(GRD.matLearn.Params(j).range, ['ML-' GRD.matLearn.Params(j).name], 2, PX);
                        end
                    end
                case 21
                     if isfield(GRD,(SVM.prog)), PXX = GRD.(SVM.prog); else, PXX=[]; end
                    switch SVM.prog
                        case {'GLMNET','GRDBST','ELASVM'}
                            GRD.(SVM.prog) = nk_GLMNET_config(SVM.prog, PXX, 0);
                        case 'ROBSVM'
                            t_act = 1; while t_act > 0, [ PXX, t_act ] = nk_ROBSVM_config(SVM.prog,PXX, NM.modeflag,0); end
                            GRD.(SVM.prog) = PXX; 
                    end
                case 22
                    SEQOPTstepsdefs =  nk_input([CutOffparstr ' range'],0,'e',SEQOPTstepsdefs);         PX = nk_AddParam(SEQOPTstepsdefs, ['ML-' CutOffparstr], 2, PX);
                case 23
                    SEQOPTlimsLdefs =  nk_input([LimsLparstr ' range'],0,'e',SEQOPTlimsLdefs);          PX = nk_AddParam(SEQOPTlimsLdefs, ['ML-' LimsLparstr], 2, PX);
                case 24
                    SEQOPTlimsUdefs =  nk_input([LimsLparstr ' range'],0,'e',SEQOPTlimsUdefs);          PX = nk_AddParam(SEQOPTlimsUdefs, ['ML-' LimsUparstr], 2, PX);
                case 25
                    CoxCutoffsdefs =  nk_input([CoxCutOffparstr ' range'],0,'e',CoxCutoffsdefs);        PX = nk_AddParam(CoxCutoffsdefs, ['ML-' CoxCutOffparstr], 2, PX);
            
                case 100
                    WLiters         = nk_input([WLparstr ' range'],0,'e',WLiters);                      PX = nk_AddParam(WLiters, ['ML-' WLparstr], 2, PX);
                case 101
                    for n = 1:NM.TrainParam.SVM.kernel.customfunc_nargin
                        argParstr = sprintf('CFuncparstr_%d', n); 
                        arg = sprintf('customkernel_arg%d', n);
                        eval(sprintf("%s = nk_input([%s ' range'], 0, 'e', %s)", arg, argParstr, arg)); 
                        eval(sprintf("PX = nk_AddParam(%s, ['ML-' %s], 2, PX)", arg, argParstr));
                    end
                          
            end
            if ~isempty(PX) && ~isempty(PX.opt), n_pars = size(PX.opt,1); else, n_pars = 0; end
        else
            n_pars = numel(Cdefs);
        end

        GRD.Cparams             = Cdefs;
        GRD.Gparams             = Gdefs;
        GRD.Epsparams           = Epsdefs;
        GRD.Nuparams            = Nudefs;
        GRD.PolyCoefparams      = PolyCoefdefs;
        GRD.PolyDegrparams      = PolyDegrdefs;
        GRD.Neuronparams        = Neurondefs;
        GRD.Kparams             = Kdefs;
        GRD.Tolparams           = Toldefs;
        GRD.Weightparams        = Weightdefs;
        GRD.Leafparams          = Leafdefs;
        % random forest
        GRD.Treeparams          = Treedefs;
        GRD.NumDparams          = NumDdefs;
        GRD.RFMaxDparams     = MaxDdefs;
        GRD.RFCritdefsparams    = Critdefs;
        GRD.RFMinSSplitparams   = MinSSplitdefs;
        GRD.RFMinSLeafparams    = MinSLeafdefs;
        GRD.RFMinWFLeafparams   = MinWeightFractLeafdefs;
        GRD.RFMaxLeafNparams    = MaxLeafNodesdefs;
        GRD.RFMinImpDecrparams  = MinImpDecrdefs;
        GRD.RFBootstrapparams   = Bootstrapdefs;
        GRD.RFOobScoreparams    = OobScoredefs;
        %GRD.RFNJobsparams       = NJobsdefs;
        GRD.RFClassWeightparams = ClassWeightdefs;
        GRD.RFCcpAlphaparams    = CcpAlphadefs;
        GRD.RFMaxSampparams     = MaxSampdefs;
     
        


        GRD.CutOffparams        = SEQOPTstepsdefs;
        GRD.LimsLparams         = SEQOPTlimsLdefs;
        GRD.LimsUparams         = SEQOPTlimsUdefs;
        GRD.CoxCutoffparams     = CoxCutoffsdefs;
        GRD.OptRegul            = OptRegul;
        GRD.NodeSelect          = NodeSelect;
        GRD.n_params            = n_pars;
        GRD.WLiters             = WLiters;
        if isfield(TrainParam.SVM.kernel, 'customfunc_nargin')
            if NM.TrainParam.SVM.kernel.customfunc_nargin > 0
                for n = 1:NM.TrainParam.SVM.kernel.customfunc_nargin
                    argName = sprintf('customkernel_arg%d', n); 
                    eval(sprintf("GRD.%s = %s;", argName, argName)); 
                end
            end
        end
    case 2
        
        % ***************** Setup for simulated annealing *****************
        
        % Model complexity influence:
        % ===========================
        GRD.lambda  = ...
            nk_input(['Lambda: influence of model complexity' ...
                    '(lambda=1: CV accuracy and complexity have same weigth)'],0,'e',1);
                
        GRD.big_gamma = nk_input('Penalizing gamma',0,'e',0.5);
        
        % SA params:
        % ==========
        GRD.init_T  = nk_input('Initial temperature',0,'e',100);
        GRD.end_T   = nk_input('End temperature',0,'e',1);
        GRD.dt      = nk_input('Default drop in temperature (0.1=fast; 0.01=slow)',0,'e',0.01);
        
        % Grid range params:
        % ==================        
        if isfield(GRD,'Cparams')
            defs = GRD.Cparams';
        else
            defs = [-5 2 15];
        end

        GRD.Cparams =  nk_input('C parameter exponent => St:Step:End',0,'e',defs,[3,1]);

        if isfield(GRD,'Gparams')
            defs = GRD.Gparams';
        else
            defs = [-15 2 3];
        end

        switch NM.SVM.kernel.typ
            case ' -t 1'
            case {' -t 2', ' --t 2'}
                GRD.Gparams =  nk_input('Gamma parameter exponent => St:Step:End',0,'e',defs,[3,1]);
            case ' -t 3'
        end
        
end

GRD.opt_probability = false;
GRD.PX = PX; 
TrainParam.GRD = GRD;

