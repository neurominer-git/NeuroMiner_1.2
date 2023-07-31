function handles = binarize_regr(handles)

m = str2double(get(handles.txtBinarize,'String')); 
if isempty(m)
    errordlg('Enter numeric threshold')    
elseif m < min(handles.labels) || m > max(handles.labels)
    errordlg(sprintf('Threshold out of target range [%g %g]',min(handles.curRegr.labels),max(handles.curRegr.labels)));
end

%% Display ROC
[handles.hroc, handles.hroc_random] = display_roc(handles);

%% Display pie charts
[handles.h1pie, handles.h2pie] = display_piecharts(handles);