function handles = binarize_regr(handles)

m = str2double(get(handles.txtBinarize,'String')); 
if isempty(m)
    errordlg('Enter numeric threshold')    
elseif m < min(handles.labels) || m > max(handles.labels)
    errordlg(sprintf('Threshold out of target range [%g %g]',min(handles.labels),max(handles.labels)));
else
    handles.curRegr.b_label = handles.labels; 
    handles.curRegr.b_label(handles.labels>=m) = 1; 
    handles.curRegr.b_label(handles.labels<m) = -1;
    handles.curRegr.b_pred = handles.Regr.mean_predictions; 
    handles.curRegr.b_pred = handles.curRegr.b_pred - m;
    [handles.curRegr.X, ...
     handles.curRegr.Y, ...
     handles.curRegr.T, ...
     handles.curRegr.AUC] = perfcurve2(handles.curRegr.b_label, handles.curRegr.b_pred, 1);
    handles.curRegr.contigmat = ALLPARAM(handles.curRegr.b_label, handles.curRegr.b_pred);
    handles.curRegr.contigmat.BINARIZATION_THRESHOLD = m;
    handles = display_contigmat(handles);
end

%% Display ROC
[handles.hroc, handles.hroc_random] = display_roc(handles);

%% Display pie charts
[handles.h1pie, handles.h2pie] = display_piecharts(handles);