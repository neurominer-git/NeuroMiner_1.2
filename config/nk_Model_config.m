function [SVM, act] = nk_Model_config(SVM, TrainParam, parentstr, varind, act)
global NM EXPERT

d = nk_GetParamDescription2(NM,TrainParam,'GridParam');
if (isfield(NM,'TrainParam') && isfield(NM.TrainParam,'RAND') && NM.TrainParam.RAND.InnerFold == -1 )
    SVM.GridParam = 1; optimstr = ''; menusel = [];
else
    optimstr = d.GridParam; optimstr = ['Define ML model performance criterion [ ' optimstr ' ]|']; menusel = 3;
end

if isfield(NM.TrainParam, 'LABEL') && NM.TrainParam.LABEL.flag
    modeflag = NM.TrainParam.LABEL.newmode;
else
    modeflag = NM.modeflag;
end

if isfield(SVM,'prog')
    d = nk_GetParamDescription2(NM,TrainParam,'SVMprog',d);
    d = nk_GetParamDescription2(NM,TrainParam,'classifier',d);
    d = nk_GetParamDescription2(NM,TrainParam,'kernel',d);
    if isempty(strfind(d.prog,'RVM'))
        progstr = [d.prog '=>' d.classifier];
    else
        progstr = d.prog;
    end
    
    progstr = [ 'Select learning algorithm [ ' progstr ' ]|'];
    if isfield(SVM,'prog') && ~any(strcmp(SVM.prog,{'kNNMEX','GLMFIT','DECTRE','RNDFOR','MEXELM','ELASVM','ROBSVM'}))
        if any(strcmp(SVM.prog,{'GLMNET','GRDBST'}))&& isfield(NM.TrainParam.GRD, SVM.prog)
            learnparstr = [];
            mnuconf = [];
        else
            if isfield(SVM, SVM.prog) || ( strcmp(SVM.prog,'MikRVM') && isfield(SVM,'RVM') )
                defstr = 'modify parameters'; 
            else
                defstr = 'parameters undefined'; 
            end
            learnparstr = [ 'Configure ' d.prog '[ ' defstr ' ]|'];
            mnuconf = 2;
        end
    else
        learnparstr = [];
        mnuconf = [];
    end
    kernelstr = ['Define kernel type [ ' d.kernel ' ]|']; 
    
    % Some learners do not operate in kernel space - skip kernel menu item
    % in these cases
    if any(strcmp(SVM.prog,{'IMRELF','kNNMEX','LIBLIN','ROBSVM','GLMFIT','MEXELM', 'DECTRE','RNDFOR','GLMNET','GRDBST','SEQOPT','WBLCOX', 'ELASVM'}))
        mnuact = [ optimstr progstr learnparstr ];
        mnusel = [ menusel 1 mnuconf];
    elseif strcmp(SVM.prog,'matLRN')
        fprintf('\nDefine possible kernel configuration required for a matLearn algorithm in the Learning algorithm parameters setup')
        mnuact = [ optimstr progstr learnparstr ];
        mnusel = [ menusel 1 mnuconf];
    else
        mnuact = [ optimstr progstr learnparstr kernelstr] ;
        mnusel = [ menusel 1 mnuconf 4];
    end
    
else
    
    if isempty(SVM)
        fprintf('No model configuration found for modality %g', varind)
        copyflag = nk_input('Do you want to use another modality''s configuration as template for current model configuration',0,'yes|no',[1,0],1);
        if copyflag
            varind_copy = nk_input('Specify modality index',0,'e');
            if varind_copy <= numel(NM.TrainParam.SVM) && ~isempty(NM.TrainParam.SVM{varind_copy})
                SVM = NM.TrainParam.SVM{varind_copy};
            end
        end
    end
    progstr ='NA';
    optimstr = 'NA';
    mnuact = [optimstr progstr];
    mnusel = [3 1];
end

if EXPERT && ~any(strcmp(SVM.prog, {'SEQOPT','WBLCOX','DECTRE'}))
    if isfield(SVM,'Post') && isfield(SVM.Post,'Detrend') && SVM.Post.Detrend
        detrendstr = 'enabled';
    else
        detrendstr = 'Disabled / NA';
    end
    switch modeflag
        case 'regression'
            mnuact = [ mnuact 'Detrend predicted targets [ ' detrendstr ']|' ]; 
            mnusel = [ mnusel 5 ];

        case 'classification'
            mnuact = [ mnuact 'Optimize decision threshold based on ROC [ ' detrendstr ' ]|' ]; 
            mnusel = [ mnusel 6 ];
    end
end

if isfield(SVM,'kernel'), kerntype = SVM.kernel.kernstr; else, kerntype = []; end

switch modeflag

    case 'classification'
        if EXPERT
            sftmenu = ['SVM --------------> LIBSVM|' ...
                       'SVM / LR ---------> LIBLINEAR|' ...
                       'SVEN -------------> Support Vector Elastic Network|' ...
                       'NEURAL -----------> Extreme Learning Machine (you have to scale data to -1/1 and prune low-var or NaN vectors) |' ...
                       'RVM / BAYES ------> Psorakis & Damoulas''s implementation with Multiple Kernel Learning|' ...
                       'RVM / BAYES ------> Relevance vector classification (MC_RVM)|' ...
                       'RVM / BAYES ------> Bayesian Sparse Logistic Regression (for classification in high-D spaces)|' ...
                       'LOCAL LEARNING ---> k-Nearest Neighbors (kNN)|' ...
                       'LOCAL LEARNING ---> IMRelief|' ...
                       'UNIVARIATE -------> GLM Logistic Regression|' ...
                       'DECISION TREE ----> MATLAB Statistics Toolbox implementation of the Decision Tree algorithm|' ...
                       'RANDOM FOREST ----> Python''s sklearn.ensemble.RandomForestClassifier algorithm|' ...
                       'matLearn ---------> Mark Schmidt''s matLearn library|'...
                       'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ...
                       'GRADIENT BOOST----> Python''s sklearn.ensemble.GradientBoostingRegressor algorithm'];
           sftval   = [ 'LIBSVM'; ...
                        'LIBLIN'; ...
                        'ELASVM'; ...
                        'MEXELM'; ...
                        'MKLRVM'; ...
                        'MVTRVR'; ...
                        'BLOREG'; ...
                        'kNNMEX'; ...
                        'IMRELF'; ...
                        'GLMFIT'; ...
                        'DECTRE'; ...
                        'RNDFOR'; ...
                        'matLRN'; ...
                        'GLMNET'; ...
                        'GRDBST'];
            if isfield(TrainParam,'STACKING') && TrainParam.STACKING.flag==1
                sftmenu = [sftmenu ...
                        '|SEQOPT -----------> Predictive sequence optimization algorithm (Stacking only)'];
                sftval = [sftval; ...
                        'SEQOPT'];
            end
            if isfield(NM,'time') && license('test', 'optimization_toolbox')
                sftmenu = [sftmenu ...
                        '|SURVIVAL MODELS --> Weibull-Cox Proportional Harzards Regression'];
                sftval = [sftval; ...
                        'WBLCOX'];
            end
        else
            sftmenu = ['SVM --------------> LIBSVM|' ...
                       'SVM/LR -----------> LIBLINEAR|' ...
                       'SVEN -------------> Support Vector Elastic Net|' ...
                       'LOCAL LEARNING ---> k-Nearest Neighbors (kNN)|' ...
                       'RVM / BAYES ------> Relevance vector classification (MC_RVM)|' ...
                       'UNIVARIATE -------> GLM Logistic Regression|' ...
                       'DECISION TREE ----> MATLAB Statistics Toolbox implementation of the Decision Tree algorithm|' ...
                       'RANDOM FOREST ----> Python''s sklearn.ensemble.RandomForestClassifier algorithm|'...
                       'matLearn ---------> Mark Schmidt''s matLearn library|'...
                       'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ...
                       'GRADIENT BOOST----> Carlos Becker''s gradient boosting algorithm'];
           
            sftval   = ['LIBSVM'; ...
                        'LIBLIN'; ...
                        'ELASVM'; ...
                        'kNNMEX'; ...
                        'MVTRVR'; ...
                        'GLMFIT'; ...
                        'DECTRE'; ...
                        'RNDFOR'; ...
                        'matLRN'; ...
                        'GLMNET'; ...
                        'GRDBST'];
            if isfield(NM,'time')
                sftmenu = [sftmenu ...
                        '|SURVIVAL MODELS --> Weibull-Cox Proportional Harzards Regression'];
                sftval = [sftval; ...
                        'WBLCOX'];
            end
        end
        
              
    case 'regression'
        if EXPERT
            sftmenu = [ 'SVM --------------> LIBSVM|' ...
                        'SVM --------------> LIBLINEAR|' ...
                        'RVM / BAYES ------> Mike Tipping''s algorithm|' ...
                        'RVM / BAYES ------> Multinomial relevance vector regression (MVRVR) (for multiple target label prediction)|' ...
                        'UNIVARIATE -------> GLM Linear Regression|' ...
                        'RANDOM FOREST ----> Python''s sklearn.ensemble.RandomForestRegressor algorithm|' ...
                        'matLearn ---------> Select a machine learning strategy from Mark Schmidt''s matLearn toolbox|' ...
                        'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ... 
                        'GRADIENT BOOST----> Python''s sklearn.ensemble.GradientBoostingRegressor algorithm'];
            sftval = ['LIBSVM'; ...
                      'LIBLIN'; ...
                      'MikRVM'; ...
                      'MVTRVR'; ...
                      'GLMFIT'; ...
                      'RNDFOR'; ...
                      'matLRN'; ...
                      'GLMNET'; ...
                      'GRDBST'];
        else
            sftmenu = [ 'SVM --------------> LIBSVM|' ...
                        'SVM --------------> LIBLINEAR|' ...
                        'RVM / BAYES ------> Mike Tipping''s algorithm|' ...
                        'UNIVARIATE -------> GLM Linear Regression|' ...
                        'RANDOM FOREST ----> Python''s sklearn.ensemble.RandomForestRegressor algorithm|' ...
                        'matLearn ---------> Select a machine learning strategy from Mark Schmidt''s matLearn toolbox|' ...
                        'GLMNET -----------> Hastie''s library for LASSO/Elastic-net regularized GLMs|' ... 
                        'GRADIENT BOOST----> Python''s sklearn.ensemble.GradientBoostingRegressor algorithm'];
            sftval = ['LIBSVM'; ...
                      'LIBLIN'; ...
                      'MikRVM'; ...
                      'GLMFIT'; ...
                      'RNDFOR'; ...
                      'matLRN'; ...
                      'GLMNET'; ...
                      'GRDBST'];
        end
end

nk_PrintLogo
mestr = 'ML algorithm parameters'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
if ~exist('act','var') || isempty(act)
    act = nk_input(mestr, 0,'mq', mnuact, mnusel,1);
end

switch act
 
    case 1

        nk_PrintLogo
        mestr = ['Define model algorithm for ' modeflag]; fprintf('\nYou are here: %s >>>',navistr);
        selProg = nk_input(mestr,0,'mq', sftmenu, 1:size(sftval,1), 1);
        if selProg, SVM.prog = sftval(selProg,:); end
        
        SVM = nk_Model_config(SVM, TrainParam, parentstr, varind, 2);
         
    case 2
        
        switch SVM.prog

            case {'LIBSVM','LIBLIN'}

                SVM = nk_SVM_config(NM, SVM, SVM.prog, kerntype, navistr);
                if strcmp(SVM.prog,'LIBLIN') && SVM.LIBLIN.b
                    SVM.RVMflag = true;
                end
                SVM = nk_Kernel_config(SVM,1); 

            case 'MikRVM'
                SVM.RVMflag = true;
                SVM.RVM.UserOpt = SB2_UserOptions;
                SVM.RVM.ParamSet = SB2_ParameterSettings;
                switch modeflag
                    case 'regression'
                       sftval = char(nk_input('Select likelihood model',0,'m', ...
                                ['Gaussian (for real valued functions)|' ...
                                 'Poisson (for positive integer value function)'],[1 2],1));
                end
                if sftval
                    sftlst = {'GAUSS','POISS'};
                    SVM.RVM.LikelihoodModel = sftlst{sftval};
                end
                SVM = nk_Kernel_config(SVM);

            case 'MKLRVM'
                SVM = nk_MKLRVM_config(SVM);
                SVM = nk_Kernel_config(SVM);

            case 'IMRELF'
                SVM.IMRELF = nk_IMRelief_config(SVM);
                SVM.kernel.kerndesc = 'other';
                SVM.kernel.kerndef = 1;
                SVM.kernel.kernstr = ' -t 2';
                % Resetting hyperparameters when IMRELIEF is chosen.
                NM.TrainParam.GRD.Cparams = [0.5 1 2 4 8 16 32]; % Sigma
                NM.TrainParam.GRD.Gparams = [0.5 .1 0.05 0.01 0.005 0.001]; % Lambda

            case 'GLMFIT'
                SVM.RVMflag = true; % Probability flag
                SVM = nk_Kernel_config(SVM);

            case 'kNNMEX'
                SVM.RVMflag = true;
                SVM = nk_Kernel_config(SVM);

            case 'KPCSVM'
                SVM = nk_KPCA_LINSVM(SVM);

            case 'BLOREG'
                SVM = nk_Kernel_config(SVM);

            case 'LSTSVM'
                SVM = nk_LSTSVM_config(SVM);

            case 'MVTRVR'
                SVM.MVTRVR.iter = nk_input('# of iterations (Less # give sparser, but less precise models)',0,'e',20);
                SVM = nk_Kernel_config(SVM);

            case 'FAMRVR'
                SVM.FAMRVR.iter = nk_input('# of iterations of the EM algorithm',0,'e',1000);
                SVM.FAMRVR.tolerance = nk_input('Tolerance',0,'e',.1);
                SVM = nk_Kernel_config(SVM);

            case {'MEXELM','DECTRE','RNDFOR'}
                switch SVM.prog
                    case 'DECTRE'
                        if ~isfield(SVM,SVM.prog), SVM.(SVM.prog) = []; end
                        SVM.Post.Detrend = 0;
                    otherwise
                        SVM = nk_Kernel_config(SVM);
                end
                
            case 'matLRN'
                if ~isfield(SVM,'matLRN'), SVM.matLRN = []; end
                switch modeflag
                    case 'classification'
                        matLRN = nk_matLearn_config(SVM.matLRN,'binaryclass',1);

                    otherwise
                        matLRN = nk_matLearn_config(SVM.matLRN,'regression',1);
                end
                if isempty(SVM.matLRN) || isempty(SVM.matLRN.algo{1}) || ~strcmp(matLRN.algo{1},SVM.matLRN.algo{1})
                    NM.TrainParam.GRD.matLearn = nk_matLearn_config(matLRN,matLRN.learner.framework,3);
                end
                SVM.matLRN = matLRN;
                SVM = nk_Kernel_config(SVM);

            case {'GLMNET','GRDBST'}
                if ~isfield(SVM,SVM.prog), SVM.(SVM.prog) = []; end
                if ~isfield(NM.TrainParam.GRD,SVM.prog), NM.TrainParam.GRD.(SVM.prog) = nk_GLMNET_config(SVM.prog, [],1, modeflag); end
                switch SVM.prog
                    case 'GLMNET'
                        switch modeflag
                            case 'classification'
                                SVM.GLMNET.family = 'binomial';
                            otherwise
                                SVM.GLMNET.family = 'gaussian';
                        end
                end

            case 'ROBSVM'
                if ~isfield(SVM,SVM.prog), SVM.(SVM.prog) = []; end
                if ~isfield(NM.TrainParam.GRD,SVM.prog), NM.TrainParam.GRD.(SVM.prog) = nk_ROBSVM_config(SVM.prog, [], modeflag, 1); end

            case 'SEQOPT'
                act = 1; if ~isfield(SVM,'SEQOPT'), [~, SVM.SEQOPT ] = nk_SEQOPT_config([],1); end
                while act 
                  [ act, SVM.SEQOPT ] = nk_SEQOPT_config( SVM.SEQOPT);             
                end

            case 'WBLCOX'
                act = 1; if ~isfield(SVM,'WBLCOX'), [~, SVM.WBLCOX] = nk_WBLCOX_config(NM, [],1); end
                while act 
                  [ act, SVM.WBLCOX ] = nk_WBLCOX_config( NM, SVM.WBLCOX );             
                end
        end
        
        
    case 3
        
        if isfield(NM.TrainParam,'LABEL') && NM.TrainParam.LABEL.flag == 1
            modeflag = NM.TrainParam.LABEL.newmode;
        else
            modeflag = [];
        end
        SVM.GridParam = nk_EvalFunc_config(NM, SVM, navistr, modeflag);
            
    case 4
        
        SVM = nk_Kernel_config(SVM);
        
    case 5
        
        SVM.Post.Detrend = nk_input('Enable post-hoc linear prediction detrending (using C1 test data)',0,'yes|no',[1,0],1);
        
    case 6
        
        SVM.Post.Detrend = nk_input('Enable post-hoc optimization of decision threshold (using ROC analysis of C1 test data)',0,'yes|no',[1,0],1);
    
end

