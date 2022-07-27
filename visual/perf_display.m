function handles = perf_display(handles)
% Set current analysis
analind = get(handles.selAnalysis,'Value'); 
handles.prevanal = handles.curranal;
handles.curranal = analind;

if ~handles.NM.analysis{analind}.status || ~isfield(handles.NM.analysis{analind},'GDdims')
    set_visibility(handles)
    set_panel_visibility(handles,'off');
    return; 
end
% Set current modality
% Set multi-modal flag
if numel(handles.NM.analysis{analind}.GDdims)>1
    handles.multi_modal = 1;
    set(handles.selModal,'Enable','on')
    load_modal(handles, handles.NM.analysis{analind}.GDdims);
    if isfield(handles.NM.analysis{analind},'META')
       popuplist = handles.selModal.String;
       popuplist{end+1}='Bagged predictor';
       handles.selModal.String = popuplist;
    end
else
    handles.multi_modal = 0;
    set(handles.selModal,'Enable','off')  
end

if handles.multilabel
    if isfield(handles.NM.analysis{analind}.params.TrainParam,'MULTILABEL')
        if handles.curlabel > numel(handles.NM.analysis{handles.curranal}.params.TrainParam.MULTILABEL.sel)
            handles.selLabel.Value = 1;
        end
        handles.selLabel.String = handles.NM.labelnames(handles.NM.analysis{handles.curranal}.params.TrainParam.MULTILABEL.sel);
        handles.curlabel = handles.NM.analysis{handles.curranal}.params.TrainParam.MULTILABEL.sel(get(handles.selLabel,'Value'));
    else
        handles.selLabel.String = handles.NM.labelnames;
        handles.curlabel= get(handles.selLabel,'Value');
    end
end
if isfield(handles,'MLIapp') && ~isfield(handles.NM.analysis{analind},'MLI')
    handles.MLIapp.delete;
    handles = rmfield(handles,'MLIapp');
elseif isfield(handles,'MLIapp') && ~isnumeric(handles.MLIapp) && isfield(handles.NM.analysis{analind},'MLI')
    appmli.updateFcn(appmli,handles);
end

[handles, visdata, oocvdata, mlidata] = switch_analysis(handles);

handles.lbStartup.String = 'Customize menus ...';
set_panel_visibility(handles,'on')
set_visibility(handles)
load_selYAxis(handles)
load_popupmenu1(handles)
load_selCVoocv(handles)
load_selModelMeasures(handles)
load_selSubParams(handles)

if ~isempty(visdata), 
    load_selModality(handles); 
    load_selPager(handles); 
else
    handles.selModelMeasures.Value = 1;
end



if isfield(handles,'MultiClass'), load_selOneVsAll_Info(handles); end
handles = display_main(handles);



