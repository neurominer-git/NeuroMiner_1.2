function [sY, sYocv, sCocv, sYw, inp, optfl] = nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, paramfl, kbin)

global PREPROC
sYocv = []; sCocv = []; sYw = [];

%Eventually, apply spatial operations on image data
if isfield(paramfl,'PREPROC') && isfield(paramfl.PREPROC,'SPATIAL') && paramfl.PREPROC.SPATIAL.cubetype>1 
    % Filter the data (imaging only)
    smoothfl = true; if isfield(inp,'issmoothed') && inp.issmoothed, smoothfl = false; end
    sY = cell(kbin,1); 
    if ~isempty(Yocv), sYocv = cell(kbin,1); end
    if ~isempty(Cocv), sCocv = cell(kbin,1); end
    if isfield(inp,'iYw'), sYw = cell(kbin,1); end
    optfl = true; 

    if smoothfl
        if inp.multiflag
            uPREPROC = nk_SetParamChain(paramfl, 1, PREPROC);
            tP = nk_ReturnParamChain(uPREPROC, 1); 
            fprintf('\nAll Predictors: Smoothing Training & CV data');
            % Smooth data
            tsY = nk_PerfSpatFilt( Y, uPREPROC, paramfl.PV ); 
            % Is there a weighting mask that must be smoothed, too?
            if isfield(inp,'iYw')
                % Here a weighting mask for a given modality "i" has been
                % extracted from inp.X(n).Yw to inp.iYw by a parent function
                tsYw = nk_PerfSpatFilt( inp.iYw, uPREPROC, paramfl.PV ); 
            else
                % Here a weighting mask needs to be found in the
                % preprocessing parameters (e.g. as part of the ranking
                % function)
                tsYw = nk_SmoothMaskInActParam(uPREPROC, paramfl.PV);
            end
            % Processing of out-of-crossvalidation data (can be multiple containers)
            if ~isempty(Yocv)
                if iscell(Yocv)
                    tsYocv = cell(numel(Yocv),1);
                    for n=1:numel(Yocv)
                        fprintf('\nSmoothing independent test data (%s)', inp.desc_oocv{n});
                        tsYocv{n} = nk_PerfSpatFilt( Yocv{n}, uPREPROC, paramfl.PV ); 
                    end
                else
                    fprintf('\nSmoothing independent test data (%s)', inp.desc_oocv);
                    tsYocv = nk_PerfSpatFilt( Yocv, uPREPROC, paramfl.PV ); 
                end
            end
            % Processing of calibration data (can be multiple containers)
            if ~isempty(Cocv)
                if iscell(Cocv)
                    tsCocv = cell(numel(Yocv),1);
                    for n=1:numel(Yocv)
                        fprintf('\nSmoothing calibration data (%g)',n)
                        tsCocv{n} = nk_PerfSpatFilt( Cocv{n}, uPREPROC, paramfl.PV ); 
                    end
                else
                    fprintf('\nSmoothing calibration data')
                    tsCocv = nk_PerfSpatFilt( Cocv, uPREPROC, paramfl.PV ); 
                end
            end

            for u=1:kbin
                paramfl.P{u} = tP;
                sY{u} = tsY;
                if isfield(inp,'iYw'), sYw{u} = tsYw; end
                if ~isempty(Yocv), sYocv{u} = tsYocv; end
                if ~isempty(Cocv), sCocv{u} = tsCocv; end
            end
        else
            for u=1:kbin
                uPREPROC = nk_SetParamChain(paramfl, u, PREPROC);
                paramfl.P{u} = nk_ReturnParamChain(uPREPROC, 1); 
                fprintf('\nPredictor #%g: Smoothing Training & CV data',u)
                inp.smoothing_kernels = uPREPROC.SPATIAL.PX.opt;
                sY{u} = nk_PerfSpatFilt( Y, uPREPROC, paramfl.PV ); 

                if isfield(inp,'iYw') && ~isempty(inp.iYw) 
                    fprintf('\nSmoothing weighting map')
                    sYw{u} = nk_PerfSpatFilt( inp.iYw, uPREPROC, paramfl.PV ); 
                else
                    sYw{u} = nk_SmoothMaskInActParam(uPREPROC, paramfl.PV);
                end
                
                if ~isempty(Yocv)
                    if iscell(Yocv)
                        for n=1:numel(Yocv)
                           fprintf('\nPredictor #%g: Smoothing independent test data (%s)', u, inp.desc_oocv{n});
                           sYocv{u,n} = nk_PerfSpatFilt( Yocv{n}, uPREPROC, paramfl.PV ); 
                        end
                    else
                        fprintf('\nPredictor #%g: Smoothing independent test data (%s)', u, inp.desc_oocv);
                        sYocv{u} = nk_PerfSpatFilt( Yocv, uPREPROC, paramfl.PV ); 
                    end
                end
                if ~isempty(Cocv)
                    if iscell(Yocv)
                        for n=1:numel(cocv)
                            fprintf('\nPredictor #%g: Smoothing calibration data (%g)',u, n)
                            sCocv{u,n} = nk_PerfSpatFilt( Cocv{n}, uPREPROC, paramfl.PV ); 
                        end
                    else
                        fprintf('\nPredictor #%g: Smoothing calibration data',u)
                        sCocv{u} = nk_PerfSpatFilt( Cocv, uPREPROC, paramfl.PV ); 
                    end
                end
            end
        end
    else
        % Check whether only some of the pre-smoothed data shelves are needed. 
        % If so extract, extract only the smoothed data shelves needed for the 
        % current CV2 partition
        for u=1:kbin
            uPREPROC = nk_SetParamChain(paramfl, u, PREPROC);
            if numel(Y{u}) ~= numel(uPREPROC.SPATIAL.PX.opt)
                idx=ismember(inp.smoothing_kernels,uPREPROC.SPATIAL.PX.opt);
            else
                idx= true(numel(Y{u}),1);
            end
            sY{u} = Y{u}(idx);
            if ~isempty(Yocv), sYocv{u} = Yocv{u}(idx); end
            if ~isempty(Cocv), sCocv{u} = Cocv{u}(idx); end
            if isfield(inp,'iYw'), sYw{u} = inp.iYw{u}(idx); end
        end
    end
else
    optfl = false;
    sY = Y;
    if ~isempty(Yocv), sYocv = Yocv; end
    if ~isempty(Cocv), sCocv = Cocv; end
    if isfield(inp,'iYw'), sYw = inp.iYw; end
end