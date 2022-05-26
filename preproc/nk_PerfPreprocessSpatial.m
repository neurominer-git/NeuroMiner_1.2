function [sY, sYocv, sCocv, inp, optfl, ukbin, uBINMOD, BINMOD] = ...
    nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, paramfl, BINMOD, kbin, ukbin)

global PREPROC
sYocv = []; sCocv = []; 

%Eventually, apply spatial operations on image data
if isfield(paramfl,'PREPROC') && isfield(paramfl.PREPROC,'SPATIAL') && paramfl.PREPROC.SPATIAL.cubetype>1 
    % Filter the data (imaging only)
    smoothfl = true; if isfield(inp,'issmoothed') && inp.issmoothed, smoothfl = false; end
    sY = cell(kbin,1); 
    if ~isempty(Yocv), sYocv = cell(kbin,1); end
    if ~isempty(Cocv), sCocv = cell(kbin,1); end
    optfl = true; ukbin = kbin; 

    if smoothfl
        if inp.multiflag
            uPREPROC = nk_SetParamChain(paramfl, 1, PREPROC);
            tP = nk_ReturnParamChain(uPREPROC, 1); 
            fprintf('\nAll Predictors: Smoothing Training & CV data');
            tsY = nk_PerfSpatFilt2( Y, uPREPROC, paramfl.PV ); 
            if isfield(inp.X,'Yw'), inp.Yw = cell(kbin,1); tsYw = nk_PerfSpatFilt2( inp.X.Yw, uPREPROC, paramfl.PV ); end
            % Processing of out-of-crossvalidation data (can be multiple containers)
            if ~isempty(Yocv)
                if iscell(Yocv)
                    tsYocv = cell(numel(Yocv),1);
                    for n=1:numel(Yocv)
                        fprintf('\nSmoothing independent test data (%s)', inp.desc_oocv{n});
                        tsYocv{n} = nk_PerfSpatFilt2( Yocv{n}, uPREPROC, paramfl.PV ); 
                    end
                else
                    fprintf('\nSmoothing independent test data (%s)', inp.desc_oocv);
                    tsYocv = nk_PerfSpatFilt2( Yocv, uPREPROC, paramfl.PV ); 
                end
            end
            % Processing of calibration data (can be multiple containers)
            if ~isempty(Cocv)
                if iscell(Cocv)
                    tsCocv = cell(numel(Yocv),1);
                    for n=1:numel(Yocv)
                        fprintf('\nSmoothing calibration data (%g)',u, n)
                        tsCocv{n} = nk_PerfSpatFilt2( Cocv{n}, uPREPROC, paramfl.PV ); 
                    end
                else
                    fprintf('\nSmoothing calibration data (%g)',u)
                    tsCocv = nk_PerfSpatFilt2( Cocv, uPREPROC, paramfl.PV ); 
                end
            end

            for u=1:kbin
                paramfl.P{u} = tP;
                sY{u} = tsY;
                if isfield(inp.X,'Yw'), inp.Yw{u} = tsYw; end
                if ~isempty(Yocv), sYocv{u} = tsYocv; end
                if ~isempty(Cocv), sCocv{u} = tsCocv; end
            end
        else
            for u=1:kbin
                uPREPROC = nk_SetParamChain(paramfl, u, PREPROC);
                paramfl.P{u} = nk_ReturnParamChain(uPREPROC, 1); 
                fprintf('\nPredictor #%g: Smoothing Training & CV data',u)
                sY{u} = nk_PerfSpatFilt2( Y, uPREPROC, paramfl.PV ); 
                if isfield(inp,'Yw') && ~isempty(inp.Yw) 
                    fprintf('\nSmoothing weighting map')
                    inp.Yw{u} = nk_PerfSpatFilt2( inp.iYw, uPREPROC, paramfl.PV ); 
                else
                    I = arrayfun( @(j) isfield(uPREPROC.ACTPARAM{j},'RANK'), 1:numel( uPREPROC.ACTPARAM ));
                    if any(I), 
                        Ix = find(I);
                        for qx = 1:numel(Ix)
                            if isfield(uPREPROC.ACTPARAM{Ix(qx)}.RANK,'EXTERN')
                                inp.Yw{u} = uPREPROC.ACTPARAM{Ix(qx)}.RANK.EXTERN;
                                inp.Yw{u} = nk_PerfSpatFilt2( inp.Yw{u}, uPREPROC, paramfl.PV ); 
                                %here, we assume that there is only one
                                %weighting map to be smoothed alongside the
                                %data. This will obviously not work for
                                %multiple weighting maps...
                                break
                            end
                        end
                    end
                end
                
                if ~isempty(Yocv), 
                    if iscell(Yocv)
                        for n=1:numel(Yocv)
                           fprintf('\nPredictor #%g: Smoothing independent test data (%s)', u, inp.desc_oocv{n});
                           sYocv{u,n} = nk_PerfSpatFilt2( Yocv{n}, uPREPROC, paramfl.PV ); 
                        end
                    else
                        fprintf('\nPredictor #%g: Smoothing independent test data (%s)', u, inp.desc_oocv);
                        sYocv{u} = nk_PerfSpatFilt2( Yocv, uPREPROC, paramfl.PV ); 
                    end
                end
                if ~isempty(Cocv), 
                    if iscell(Yocv)
                        for n=1:numel(cocv)
                            fprintf('\nPredictor #%g: Smoothing calibration data (%g)',u, n)
                            sCocv{u,n} = nk_PerfSpatFilt2( Cocv{n}, uPREPROC, paramfl.PV ); 
                        end
                    else
                        fprintf('\nPredictor #%g: Smoothing calibration data',u)
                        sCocv{u} = nk_PerfSpatFilt2( Cocv, uPREPROC, paramfl.PV ); 
                    end
                end
            end
        end
    else
       
        sY = Y;
        if ~isempty(Yocv), sYocv = Yocv; end
        if ~isempty(Cocv), sCocv = Cocv; end
        if isfield(inp,'X') && isfield(inp.X,'Yw') && ~isfield(inp,'Yw'), inp.Yw = inp.X.Yw; end
    end
else
    %uBINMOD = BINMOD;
    optfl = false;
    sY = Y;
    if ~isempty(Yocv), sYocv = Yocv; end
    if ~isempty(Cocv), sCocv = Cocv; end
    if isfield(inp,'X') && isfield(inp.X,'Yw') && ~isfield(inp,'Yw'), inp.Yw = inp.X.Yw; end
end
uBINMOD = BINMOD;