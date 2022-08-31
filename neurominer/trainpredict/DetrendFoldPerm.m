function OUT = DetrendFoldPerm(IN, OUT)

global SVM MODEFL EVALFUNC

if isfield(SVM,'Post') && isfield(SVM.Post,'Detrend') && SVM.Post.Detrend
    OUT.detrend.flag = true;
else
    OUT.detrend.flag = false;
end

if OUT.detrend.flag

    PermVec = 1:IN.nperms;
    FoldVec = 1:IN.nfolds;
    ClassVec = 1:IN.nclass;

    PermNum = numel(PermVec);
    FoldNum = numel(FoldVec);
    ClassNum = numel(ClassVec);
    thresh = zeros(PermNum * FoldNum,ClassNum); aucs=thresh;
    beta = thresh; p=thresh;
    ll=1;
    for ii=1:PermNum % Loop through CV1 permutations

        for jj=1:FoldNum % Loop through CV1 folds

            i = PermVec(ii); j= FoldVec(jj);

            for ccurclass=1:ClassNum % Loop through dichotomizers

                curclass = ClassVec(ccurclass);
                CVInd = IN.Y.CVInd{i,j}{curclass};
                zu = size(OUT.Trtargs{i,j,curclass},2);
                
                P = []; D = []; L = [];

                for z = 1 : zu
                    
                    D = [ D; OUT.CVdecs{i,j,curclass}(CVInd,z) ] ;
                    P = [ P; OUT.CVtargs{i,j,curclass}(CVInd,z) ] ;
                    L = [ L; IN.Y.CVL{i,j}{curclass} ] ; 
                    
                end

                switch MODEFL
                    case 'regression'
                        [~, beta(ll), p(ll)] = nk_DetrendPredictions2([], [], P, L);
                    case 'classification'
                        L(L==-1)=0; X = [D L]; I=sum(isnan(X),2);
                        ROCout = roc2(X(~I,:),[],[],0);
                        aucs(ll,curclass)= ROCout.AUC; thresh(ll,curclass)=ROCout.co;
                end
                ll=ll+1;

            end
        end
    end
    
    switch MODEFL
        case 'regression'
            OUT.detrend.beta = mean(beta);
            OUT.detrend.p = mean(p);
            
        case 'classification'
            OUT.detrend.thresh = mean(thresh);

             for ii=1:PermNum % Loop through CV1 permutations
                for jj=1:FoldNum % Loop through CV1 folds
                    i = PermVec(ii); j= FoldVec(jj);
                    for ccurclass=1:ClassNum % Loop through dichotomizers
                        curclass = ClassVec(ccurclass);
                        zu = size(OUT.Trtargs{i,j,curclass},2);
                        TrInd = IN.Y.TrInd{i,j}{curclass};
                        CVInd = IN.Y.CVInd{i,j}{curclass};
                        for z = 1 : zu
                            OUT.Trdecs{i,j,curclass}(:,z) =  OUT.Trdecs{i,j,curclass}(:,z) - OUT.detrend.thresh(curclass) ;
                            OUT.Trtargs{i,j,curclass}(:,z) = sign( OUT.Trdecs{i,j,curclass}(:,z) );
                            OUT.CVdecs{i,j,curclass}(:,z) =  OUT.CVdecs{i,j,curclass}(:,z) - OUT.detrend.thresh(curclass) ;
                            OUT.CVtargs{i,j,curclass}(:,z) = sign( OUT.CVdecs{i,j,curclass}(:,z) );
                            OUT.tr{i,j,curclass}(:,z) = EVALFUNC( IN.Y.TrL{i,j}{curclass}, OUT.Trdecs{i,j,curclass}(TrInd,z) );
                            OUT.ts{i,j,curclass}(:,z) = EVALFUNC( IN.Y.CVL{i,j}{curclass}, OUT.CVdecs{i,j,curclass}(CVInd,z) );
                        end
                    end
                end
            end
    end

end

end