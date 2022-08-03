function [ YYY, LLL, CCC, fl ] = nk_PerfADASYN(Y, L, IN, C, noconcatfl)

adasyn_beta                     = [];   %let ADASYN choose default
adasyn_kDensity                 = [];   %let ADASYN choose default
adasyn_kSMOTE                   = [];   %let ADASYN choose default
adasyn_normalized               = false;%false lets ADASYN handle normalization
cfl                             = false;
fl                              = false;

if exist('IN','var') && ~isempty('IN')
    if isfield(IN,'beta') && ~isempty(IN.beta), adasyn_beta = IN.beta; end
    if isfield(IN,'kDensity') && ~isempty(IN.kDensity), adasyn_kDensity = IN.kDensity; end
    if isfield(IN,'kSMOTE') && ~isempty(IN.kSMOTE), adasyn_kSMOTE = IN.kSMOTE; end
    if isfield(IN,'normalized') && ~isempty(IN.normalized), adasyn_normalized = IN.normalized; end
end

if exist('C','var') && ~isempty(C)
    cfl = true; nC = size(C,2); 
else
    CC = []; 
end

if ~exist('noconcatfl','var') || isempty(noconcatfl)
    noconcatfl = false;
end

U = nk_CountUniques(L(~isnan(L)));
[~,maxLidx] = max(U.UN{1});
refcl = U.UX{1}(maxLidx);
uL = U.UX{1}; uL(uL==refcl) = []; nL = numel(uL);
    
if iscell(Y)
    
    N = numel(Y);
    LLL = cell(N,1);
    YYY = cell(N,1);
    if cfl, CCC = cell(N,1); end
    % Loop through training data shelves
    for i=1:N
        
        YY = []; LL = [];
        if cfl, CC=[]; end
    
        for cl = 1:nL
        
            idx1 = L==uL(cl);
            idx2 = L==refcl;
                
            % If covariates are available add them to training matrix
            if cfl
                nY = size(Y{i},2); tY = [ [Y{i}(idx1,:) C(idx1,:)]; [Y{i}(idx2,:) C(idx2,:)] ];
            else
                tY = [Y{i}(idx1,:); Y{i}(idx2,:)];
            end
            
            % Perform ADASYN
            [ tYsyn, Lsyn ] = ADASYN(tY, L, adasyn_beta, adasyn_kDensity, adasyn_kSMOTE, adasyn_normalized); 
            
            % Check whether covariates were integrated in training data and
            % split synthetic data into synthetic training and covariate matrices
            if ~isempty(tYsyn)
                if cfl 
                    Ysyn = tYsyn(:,1:nY); Csyn = tYsyn(:,nY+1:nY+nC);
                    CC  = [CC Csyn];
                else
                    Ysyn = tYsyn;
                end
                Lsyn(Lsyn==1) = uL(cl);
                Lsyn(Lsyn==0) = refcl;
    
                YY = [YY; Ysyn]; LL = [LL; Lsyn]; 
                fl = true;
            end
        end
        if ~noconcatfl
            YYY{i} = [Y; YY];
            LLL{i} = [L; LL];
            if cfl, CCC{i} = [C; CC]; end
        else
            YYY{i} = YY;
            LLL{i} = LL;
            if cfl, CCC{i} = CC; end
        end
    end

else
    
    YY = []; LL = []; CCC=[];
    if cfl, CC=[];  end

    for cl=1:nL
        
        idx1 = L==uL(cl);
        idx2 = L==refcl;

        % If covariates are available add them to training matrix
        if cfl
            nY = size(Y,2); tY = [ [Y(idx1,:) C(idx1,:)]; [Y(idx2,:) C(idx2,:)] ];
        else
            tY = [Y(idx1,:); Y(idx2,:)];
        end
        tL = [ones(sum(idx1),1); zeros(sum(idx2),1)];
        
        % Perform ADASYN
        [tYsyn, Lsyn] = ADASYN(tY, tL, adasyn_beta, adasyn_kDensity, adasyn_kSMOTE, adasyn_normalized);
        
        % Check whether covariates were integrated in training data and
        % split synthetic data into synthetic training and covariate matrices
        if ~isempty(tYsyn)
            if cfl 
                Ysyn = tYsyn(:,1:nY); Csyn = tYsyn(:,nY+1:nY+nC);
                CC  = [CC; Csyn];
            else
                Ysyn = tYsyn;
            end
            Lsyn(Lsyn==1) = uL(cl);
            Lsyn(Lsyn==0) = refcl;

            YY = [YY; Ysyn]; LL = [LL; Lsyn]; 
            fl = true;
        end
    end
    if ~noconcatfl
        YYY = [Y; YY];
        LLL = [L; LL];
        if cfl, CCC = [C; CC]; end
    else
        YYY = YY;
        LLL = LL;
        if cfl, CCC = CC; end
    end
end