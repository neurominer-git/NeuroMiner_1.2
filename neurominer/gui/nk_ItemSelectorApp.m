function varargout = nk_ItemSelectorApp(varargin)
app = nk_ItemSelector_App(varargin{:});
%waitfor(app.cmdFinish,'UserData')
if isvalid(app)
    varargout = app.output;
    app.delete
end
