function load_selCase(handles,cases)

popupstr = cellstr(join(string([cellstr(num2str((1:numel(cases))')) cases]),' | '));
handles.selCase.String = popupstr;
handles.selCase.Value = 1;
