function [ Lperm, Yperm ] = nk_VisXPermY(Y, L, IND, permmode, indpermrows, indpermcols, permind, analysis, inp, paramfl, BINMOD, n, h, k, l, pnt, FullPartFlag, F, u)
global MODEFL

switch permmode

    case 1

        Lperm = L(indpermrows(:,permind));
        Yperm = Y;

    case {2,3}

        if permmode == 3
            Lperm = L(indpermrows(:,permind));
        else
            Lperm = L(IND);
        end
        L = L(IND); Y=Y(IND,:);
        indN = isnan(L);
        uL = unique(L(~indN)); nuL = numel(uL);
        Yperm = zeros(size(Y));
        for i=1:nuL
            indi = L == uL(i);
            Yperm(indi,:) = Y(indi,indpermcols(i,:,permind));
        end

        if strcmp(MODEFL,'classification')
            tL = Lperm;
            for i=1:nuL
                if i==1
                    tL( Lperm == uL(1) ) = 1;
                else
                    tL( Lperm == uL(i) ) = -1;
                end
            end
            Lperm = tL;
        end
        if sum(indN)
            Lperm(indN)=[]; Yperm(indN,:)=[];
        end

    case 4
        % MD: hard coded covar indices
        chosen_covInd = inp.VIS.PERM.covars_idx;
        % MD: permute the site labels
        covars_p = md_CovarsPerm_covInd(analysis.datadescriptor.input_settings.covars,chosen_covInd);
        analysis_p = analysis;
        analysis_p.datadescriptor.input_settings.covars =covars_p;
        inp_p = inp;
        inp_p.covars = covars_p;

        % MD: apply preprocessing with the permuted site labels
        [ ~, ~, ~, mapY_p, ~, ~, Param_p, ~, ~ ] = nk_ApplyTrainedPreproc(analysis_p, inp_p, paramfl);

        % MD: extract the feature data
        % get CV1 training and test data
        CVInd   = mapY_p.CVInd{k,l}{h};
        TrInd   = mapY_p.TrInd{k,l}{h};

        if BINMOD, hix = h; else, hix =1; end
        [ modelTr , modelCV, ~ ] = nk_ReturnAtOptPos(mapY_p.Tr{k,l}{hix},  mapY_p.CV{k,l}{hix}, mapY_p.Ts{k,l}{hix}, [], Param_p{1}(k,l,hix), pnt);
        switch inp.FUSION.flag
            case 2
                for n = 1:nM
                    [ modelTr, modelCV, ~ ] = nk_ReturnAtOptPos(mapY_p.Tr{k,l}{hix},  mapY_p.CV{k,l}{hix}, mapY_p.Ts{k,l}{hix}, [], Param_p{n}(k,l,hix), pnt);
                end
            otherwise
                [ modelTr, modelCV, ~ ] = nk_ReturnAtOptPos(mapY_p.Tr{k,l}{hix},  mapY_p.CV{k,l}{hix}, mapY_p.Ts{k,l}{hix}, [], Param_p{1}(k,l,hix), pnt);
        end

        % Concatenate Training and CV data if needed
        if FullPartFlag
            modelTr = [modelTr; modelCV ];
            TrInd = [TrInd; CVInd];
        end
        modelTr = modelTr(TrInd,:);

        % MD: Extract feature data
        try
            Ymodel = nk_ExtractFeatures(modelTr, F, [], u);
        catch
            error('Dimensionality mismatch between training data matrix and feature selection mask. Check your settings')
        end
        Yperm = Ymodel;
        Lperm = L(IND);
end

function covars = md_CovarsPerm_covInd(covars,chosen_covInd)
% create permutation indices
cn = size(covars,1);
p = randperm(cn);
% Permute each covar column according to the
% permutation indices
for ck = 1:size(chosen_covInd,2)
    covars(:,chosen_covInd(ck)) = covars(p,chosen_covInd(ck));
end

