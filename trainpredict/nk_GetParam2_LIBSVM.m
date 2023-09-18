% =========================================================================
% FORMAT function [param, model] = nk_GetParam_LIBSVM(Y, label, ModelOnly, 
%                                                                 ...cmdstr)
% =========================================================================
% Train LIBSVM models and evaluate their performance using Y & label, 
% SlackParam, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 11/2020

function [param, model] = nk_GetParam2_LIBSVM(Y, label, ModelOnly, cmd)
         
global SVM EVALFUNC LIBSVMTRAIN MODEFL CMDSTR GK 

param = []; flw = 0; GK = struct('gkernelBool',0); oneclass = false;
if contains(cmd,'-s 2'), oneclass = true; end
if strcmp(LIBSVMTRAIN,'svmtrain312'), flw = 1; end
if flw
    switch SVM.kernel.kerndef 
        case 5
            GK.gkernelBool = 1; GK.gkernelFunction = 'WL'; GK.iter = CMDSTR.WLiter; 
        case 6
            GK.gkernelBool = 1; GK.gkernelFunction = 'WLspdelta'; GK.iter = CMDSTR.WLiter; 
        case 7
            GK.gkernelBool = 1; GK.gkernelFunction = 'WLedge'; GK.iter = CMDSTR.WLiter; 
        case 8 
            GK.evalStr = strings; 
            if SVM.kernel.customfunc_nargin > 0
                for n = 1:SVM.kernel.customfunc_nargin
                    argName = sprintf('customkernel_arg%d', n); 
                    arg_i = CMDSTR.(argName);
                    if CMDSTR.(argName)
                        GK.evalStr = sprintf('%s, %d' , GK.evalStr, arg_i);
                    end
                end
            end
    end 
end

% Check if weighting is necessary
cmd = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmd);

% Check if sample weighting is necessary (currently regression only)
if flw
    W = ones(numel(label),1);
    if strcmp(MODEFL,'regression') && SVM.LIBSVM.Weighting
        W = nk_WeigthDataInstanceHisto(label);
    end
end

if iscell(Y) 
   
    % MKL-based learning not implemented yet
   
else % Univariate case
    if GK.gkernelBool == 1
        Y = GraphKernel_matrixInput(Y,Y, GK.gkernelFunction, GK.iter);
        numTrain = size(Y,1);
        Y = [(1:numTrain)' , Y];
    end
    if SVM.kernel.kerndef == 8
        Y = eval(sprintf('feval(SVM.kernel.customfunc, Y, Y %s)', GK.evalStr));
        Y = [(1:numTrain)' , Y];
    end
    if size(label,1) ~= size(Y,1), label = label'; end
    if SVM.RVMflag, label(label==-1)=2; end
    if ~SVM.LIBSVM.Weighting && (isfield(SVM,'AdaBoost') && SVM.AdaBoost.flag)        
        N = length(label); % X training labels
        W = 1/N * ones(N,1); %Weights initialization
        for m=1:SVM.AdaBoost.BoostIter
            %Calculate the error and alpha in adaBoost with cross validation
            model               = svmtrain312( W, label, Y, cmd ); 
            Xout                = svmpredict312( label, Y, model );
            s1                  = sum( (label ==  1) .* (Xout) .* W);
            s2                  = sum( (label == -1)  .* (Xout) .* W);
            if s1 == 0 && s2 == 0, break; end
            % Compute Alpha
            alpha              = 0.5*log((s1 + eps) / (s2+eps));  
            % update the weight
            W                   = exp( -1 * (label .* Xout .* alpha));
            W                   = W/norm(W);
        end  
    else
        if flw
            if oneclass
                idx = label == -1;
                model = feval( LIBSVMTRAIN, W(idx), label(idx), Y(idx,:), cmd);
            else
                model = feval( LIBSVMTRAIN, W, label, Y, cmd);
            end
            model = nk_BuildCalibrationModel(SVM, MODEFL, model, Y, label);
        else
            if oneclass
                idx = label == -1;
                model = feval( LIBSVMTRAIN, label(idx), Y(idx,:), cmd);
            else
                model = feval( LIBSVMTRAIN, label, Y, cmd );
            end
            model = nk_BuildCalibrationModel(SVM, MODEFL, model, Y, label);
        end
    end
    if ~ModelOnly
        [param.target, param.dec_values] = nk_GetTestPerf_LIBSVM([], Y, label, model);
        if SVM.RVMflag, param.dec_values = nk_CalibrateProbabilities(param.dec_values); end
        param.val = EVALFUNC(label, param.target);
    end
end
