% =========================================================================
% FORMAT NM = nk_OOCVDataIO_config(NM, O, parentstr)
% =========================================================================
% This function manages the input of independent test data into NM
%
% Inputs:
% -------
% NM        : The NM workspace
% O         : The OOCV workspace
% parentstr : Name of the calling function
%
% Outputs:
% --------
% NM (see above)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 07/2023

function NM = nk_OOCVDataIO_config(NM, O, parentstr)
    
if (~exist('O','var') || isempty(O)) || isnumeric(O)
    [NM, Y, inp.oocvind, inp.fldnam, inp.dattype] = nk_SelectOOCVdata(NM, O, 1);
else
    inp.fldnam      = O.fldnam;
    inp.oocvind     = O.ind;
    Y               = NM.(inp.fldnam){inp.oocvind};
    Y.desc          = O.desc;
    Y.date          = O.date;
    Y.label_known   = false;
    switch O.fldnam
        case 'OOCV'
            inp.dattype     = 'independent test data';
        case 'C'
            inp.dattype     = 'calibration data';
    end
end
if isempty(inp.oocvind), return; end
NM.(inp.fldnam){inp.oocvind}.desc = Y.desc;
NM.(inp.fldnam){inp.oocvind}.date = Y.date;
if ~isfield(NM.(inp.fldnam){inp.oocvind},'n_subjects_all'), NM.(inp.fldnam){inp.oocvind}.n_subjects_all = Inf; end
if ~isfield(NM.(inp.fldnam){inp.oocvind},'labels_known') && isfield(Y,'labels_known') || (NM.(inp.fldnam){inp.oocvind}.labels_known ~=  Y.labels_known)
    NM.(inp.fldnam){inp.oocvind}.labels_known = Y.labels_known; 
end
inp.currmodal = 1;
inp.nummodal = numel(NM.Y);
inp.na_str = '?';
inp.covflag = false;
inp.desc = Y.desc;
inp.Ydims = cellfun('size', NM.Y , 2);
inp.datadescriptor = NM.datadescriptor;
inp.brainmask = NM.brainmask;
inp.badcoords = NM.badcoords;
inp.analysis = NM.analysis;
inp.oocvmode = true;
inp.id = NM.id;
if isfield(NM,'covars'), inp.covflag = true; inp.covars = NM.covars; end

act = 1; while act, [act, Y, inp] = nk_ManageDataContainers(act, inp, Y, parentstr); end

NM.(inp.fldnam){inp.oocvind} = Y;

