function [handles, oocvind] = get_oocvind(handles)

% Get index to OOCV data
oocvind = handles.selCVoocv.Value - 1;

if ~isfield(handles,'oocvind'), handles.oocvind = oocvind; end

% Try to find OOCV data pointer in newly selected analysis.
if handles.curranal ~= handles.prevanal && handles.OOCVinfo.Analyses{handles.prevanal}.OOCVdone
    descprev = handles.OOCVinfo.Analyses{handles.prevanal}.descriptor{handles.oocvind};
    descnew = handles.OOCVinfo.Analyses{handles.curranal}.descriptor;
    I = find(strcmp(descnew,descprev));
    if ~isempty(I) && I ~= oocvind 
        handles.selCVoocv.Value = I+1;
        oocvind = I;
    end
end
handles.prevanal = handles.curranal;