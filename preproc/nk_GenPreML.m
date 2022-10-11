function [PX, PreML] = nk_GenPreML(PREPROC)

PX = nk_ReturnParamChain(PREPROC, 1); PreML = [];

if ~isempty(PX)
    nP = numel(PX);
    if nP>1
        Params = []; Params_desc=[]; ModalityVec = [];
        for n=1:nP
            Params = [Params PX(n).Params];
            for m=1:numel(PX(n).Params_desc)
                PX(n).Params_desc{m} = sprintf('%s_{M%g}',PX(n).Params_desc{m},n);
            end
            Params_desc = [Params_desc PX(n).Params_desc] ;
            ModalityVec = [ModalityVec n*ones(1,numel(PX(n).Params))];
        end
        PreML.Params_desc   = Params_desc;
        PreML.Params        = Params;
        PreML.ModalityVec   = ModalityVec;
    else
        PreML.Params_desc   = PX.Params_desc;
        PreML.Params        = PX.Params;
        PreML.ModalityVec   = ones(1,numel(PX.Params));
    end
end