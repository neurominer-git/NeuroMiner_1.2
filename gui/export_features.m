function [ERR, TBL] = export_features(handles, batchmode)

warning off

curclass = get(handles.popupmenu1,'Value');
filename = sprintf('%s_A%g_Features', handles.params.TrainParam.SAV.matname, handles.curranal);
[nlabels, nmodal] = size(handles.visdata_table);

if strcmp(handles.modeflag,'regression')
    ModeStr = 'Rgr';
else
    ModeStr = sprintf('Cl%g',curclass);
end

for i=1:nlabels
    for j=1:nmodal
        if nmodal > 1 && nlabels >1
             SpreadSheetName = sprintf('%s-M%g-L%g', ModeStr, j, i);
        elseif nmodal > 1
             SpreadSheetName = sprintf('%s-M%g', ModeStr, j);
        elseif nlabels > 1
             SpreadSheetName = sprintf('%s-L%g', ModeStr, i);
        else
             SpreadSheetName = sprintf('%s', ModeStr);
        end
        TBL = handles.visdata_table(i, j).tbl(curclass);
        ERR = tbl2file(TBL, filename, SpreadSheetName);
    end
end

if isempty(ERR)
    if ~batchmode
        msgbox(['Feature data successfully exported to file ' filename]);
    end
else
    warndlg(['Feature data NOT successfully exported to file ' filename]);
end

warning on