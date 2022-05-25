function nkdev_MakeTransparent(Y, L, IN)
global SVM

if ~isfield(IN, 'RAND') || isempty(IN.RAND)
	% Setup permutation parameters
	IN.RAND.OuterPerm = 1;
	IN.RAND.InnerPerm = 1;
	IN.RAND.OuterFold = 5;
	IN.RAND.InnerFold = 5;
	IN.RAND.Decompose = 1;
end

if ~isfield(IN,'algorithm') || isempty(IN.algorithm)
	IN.algorithm = 'LINKERNSVM';
end

% Define algorithm to use for the simulation
LIBSVMTRAIN = @svmtrain312;
LIBSVMPREDICT = @svmpredict312;

switch IN.algorithm
	case 'LINKERNSVM'
		SVM.prog = 'LIBSVM';
		SVM.(SVM.prog).Weighting = 1;
		SVM.(SVM.prog).Optimization.b = 0;
		SCALEMODE = 'scale';
	case 'LINSVM'
		SVM.prog = 'LIBLIN';
		SVM.(SVM.prog).b = 0;
		SVM.(SVM.prog).classifier = 3;
        SVM.(SVM.prog).tolerance = 0.01;
		CMDSTR.WeightFact = 1;
		SCALEMODE = 'std';
	case 'L2LR'
		SVM.prog = 'LIBLIN';
		SVM.(SVM.prog).b = 0;
		SVM.(SVM.prog).classifier = 0;
        SVM.(SVM.prog).tolerance = 0.01;
		CMDSTR.WeightFact = 1;
		SCALEMODE = 'std';
	case 'L1SVC'
		SVM.prog = 'LIBLIN';
		SVM.(SVM.prog).b = 0;
		SVM.(SVM.prog).classifier = 5;
        SVM.(SVM.prog).tolerance = 0.01;
		SCALEMODE = 'std';
    case 'L1LR'
		SVM.prog = 'LIBLIN';
		SVM.(SVM.prog).b = 0;
		SVM.(SVM.prog).classifier = 6;
        SVM.(SVM.prog).tolerance = 0.01;
		SCALEMODE = 'std';
	case 'RF' % not implemented yet
end

% We basically always want to weight the decision boundary
SVM.(SVM.prog).Weighting = 1;
CMDSTR.WeightFact = 1;

% At the moment only binary classification is supported 
MODEFL = 'classification';

% We want to repeat this ten times for producing more stable results
reps = 10;
pCV = RAND.OuterPerm;
nCV = RAND.OuterFold;
R = zeros(reps,pCV*nCV);

xSVM = SVM;
Lx = L; Lx(L==-1) = 2;
cv = nk_MakeCrossFolds(Lx, RAND, 'classification', [], {'A','B'} );
ncm = any(isnan(Y(:)));
for j=1:pCV
	
		% We could use here a multi-core to parallelize the CV cycle
        for i=1:nCV
            
            Tr = cv.TrainInd{j,i};
            Ts = cv.TestInd{j,i};					
			Lr = L(Tr);
            Ls = L(Ts);
            W = ones(numel(Lr),1);

            % Scale / standardize matrix using NM and do it properly within the cross-validation cycle
			switch SCALEMODE
				case 'scale'
                    if verbose, fprintf('\nScale data matrix.'); end
					[Tr_Y, INi] = nk_PerfScaleObj(Y(Tr,:)); % Train scaling model
					Ts_Y = nk_PerfScaleObj(Y(Ts,:), INi); % Apply it to the test data
				case 'std'
                    if verbose, fprintf('\nStandardize data matrix.'); end
					[Tr_Y, INi] = nk_PerfStandardizeObj(Y(Tr,:)); % Train standardization model
					Ts_Y = nk_PerfStandardizeObj(Y(Ts,:), INi); % Apply it to the test data
            end
            
            if ncm
                if verbose, fprintf('\nImpute data matrix.'); end
                [Tr_Y, INi] = nk_PerfImputeObj(Tr_Y); % Train imputation model
                Ts_Y = nk_PerfImputeObj(Ts_Y, INi); % Apply it to the test data
            end

            % Here we need to create modified version of the test data 
            switch algorithm
                case 'LINKERNSVM'
                    % standard LIBSVM params:
                    cmd = '-s 0 -t 0 -c 1';
                    % set weighting!
                    cmd = nk_SetWeightStr(xSVM, MODEFL, CMDSTR, Lr, cmd);
                    % Train model
                    model = LIBSVMTRAIN(W, Lr, Tr_Y, cmd);
                   
                case {'LINSVM', 'L2LR', 'L1LR', 'L1SVC'}
                    %Define command string
                    cmd = nk_DefineCmdStr(xSVM, MODEFL); 
                    cmd = [ ' -c 1' cmd.simplemodel cmd.quiet ];
                    % set weighting!
                    cmd = nk_SetWeightStr(xSVM, MODEFL, CMDSTR, Lr, cmd);
                    % Train model
                    model = train_liblin244(Lr, sparse(Tr_Y), cmd);
            end

            for h=1:size(Ts_Y)
                h_Ls = Ls(h,:);
                [h_Ts_Y_pos, h_Ts_Y_neg] = nkdev_CreateData4ModelInterpreter(Tr_Y, Ts_Y(h,:), IN.nperms, IN.frac, model); 
                switch algorithm
                    case  'LINKERNSVM'
                         [ ~, ~, ds_pos ] = LIBSVMPREDICT( h_Ls, h_Ts_Y_pos, model, sprintf(' -b %g',xSVM.LIBSVM.Optimization.b));
                         [ ~, ~, ds_neg ] = LIBSVMPREDICT( h_Ls, h_Ts_Y_neg, model, sprintf(' -b %g',xSVM.LIBSVM.Optimization.b));
                    case {'LINSVM', 'L2LR', 'L1LR', 'L1SVC'}
                         [ ~, ~, ds_pos ] = predict_liblin244(h_Ls, sparse(h_Ts_Y_pos), model, sprintf(' -b %g -q',xSVM.LIBLIN.b)); 
                         [ ~, ~, ds_neg ] = predict_liblin244(h_Ls, sparse(h_Ts_Y_neg), model, sprintf(' -b %g -q',xSVM.LIBLIN.b)); 
                         ds_pos = ds_pos(:,1);
                         ds_neg = ds_neg(:,1);
                end
            end

        end
end