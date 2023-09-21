function load_selPager(handles, perpage)

if ~exist('perpage','var') || isempty(perpage), perpage = 50; end
if strcmp(handles.popupmenu1.String{handles.popupmenu1.Value},'Multi-group classifier')
    multiflag = true;
else
    multiflag = false;
end

if handles.curmodal<= size(handles.visdata,1) && handles.visdata{handles.curmodal, handles.curlabel}.params.visflag ~= 1 && ~multiflag 
    nfeats = handles.visdata{handles.curmodal,handles.curlabel}.params.nfeats;
    if nfeats > perpage
        vec = 1:perpage:nfeats;
        if vec(end) < nfeats, vec(end+1) = nfeats; end
        for i=1:numel(vec)-1
            popuplist{i} = sprintf('%g:%g',vec(i),vec(i+1));
        end
        popuplist{end} = sprintf('%g:%g',vec(1),vec(end));
        set(handles.selPager, 'String', popuplist); 
        set(handles.selPager,'Enable','on');
    else
        popuplist{1} = sprintf('1:%g',nfeats);
        set(handles.selPager, 'String', popuplist); 
        set(handles.selPager,'Enable','off');
    end
    set(handles.tglSortFeat,'Enable','on');
    set(handles.cmdExportFeats,'Enable','on');
else
    set(handles.selPager,'Enable','off');
    set(handles.tglSortFeat,'Enable','off');
    set(handles.cmdExportFeats,'Enable','off');
end