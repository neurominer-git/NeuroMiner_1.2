function [handles, visdata, oocvdata, mlidata] = switch_analysis(handles)

visdata = []; oocvdata = []; mlidata = [];
analind = handles.curranal;
handles.curmodal = 1; if strcmp(handles.selModal.Enable,'on'), handles.curmodal = handles.selModal.Value; end

if handles.multi_modal && strcmp(handles.selModal.String{handles.selModal.Value},'Bagged predictor')
    GDdims = handles.NM.analysis{analind}.META; handles.METAstr = 'bagged';
else
    GDdims = handles.NM.analysis{analind}.GDdims{handles.curmodal}; handles.METAstr = 'none';
end

% Remove subgroup indices
if isfield(handles,'SubIndex') 
    handles = rmfield(handles,'SubIndex');
end

% Set current label
if size(handles.NM.label,2)>1, handles.multilabel = true; else, handles.multilabel=false; end
handles.curlabel = get(handles.selLabel,'Value');

% Check whether subindex is activated

% set alternative label
if isfield(handles.NM.analysis{1,analind}.params,'label')
    handles.label = handles.NM.analysis{1,analind}.params.label.label;
    handles.modeflag = handles.NM.analysis{1,analind}.params.label.modeflag;
else
    handles.label = handles.NM.label;
    handles.modeflag = handles.NM.modeflag;
end

% Check whether selected analysis has visualisation data
if isfield(handles.NM.analysis{analind},'visdata')
    visdata = handles.NM.analysis{analind}.visdata; 
    
elseif isfield(handles,'visdata')
    handles = rmfield(handles,'visdata');
end

% Check whether selected analysis has OOCV data
if isfield(handles.NM.analysis{analind},'OOCV') && handles.NM.defs.analyses_locked
    idx = ~cellfun(@isempty,handles.NM.analysis{analind}.OOCV);
    oocvdata = handles.NM.analysis{analind}.OOCV(idx); 
    handles.oocvview = true;
    [handles, oocvind] = get_oocvind(handles);
elseif isfield(handles,'OOCV')
    handles = rmfield(handles,'OOCV');
    oocvdata = [];
    handles.oocvview = false;
    handles.selCVoocv.Value=1;
end

load_selCVoocv(handles)

% Check whether selected analysis has MLI data
mlifl='off';
if isfield(handles,'oocvview') && handles.oocvview && oocvind
    if isfield(handles.NM.analysis{analind}.OOCV{oocvind},'MLI')
        mlidata = handles.NM.analysis{analind}.OOCV{oocvind}.MLI;
        mlifl = 'on';
    elseif isfield(handles,'MLI')
        handles = rmfield(handles,'MLI');
    end
else
    if isfield(handles.NM.analysis{analind},'MLI')
        mlidata = handles.NM.analysis{analind}.MLI;
        mlifl = 'on';
    elseif isfield(handles,'MLI')
        handles = rmfield(handles,'MLI');
    end
end
handles.thisMLIresult.Visible = mlifl;

handles = load_analysis(handles, ...
                    'Subjects', handles.NM.cases, ...
                    'Params', handles.NM.analysis{analind}.params, ...
                    'Analysis', handles.NM.analysis{analind}.params.cv, handles.label, GDdims, ...
                    'Visdata', visdata, ...
                    'OOCVdata', oocvdata, ...
                    'MLIdata', mlidata);
