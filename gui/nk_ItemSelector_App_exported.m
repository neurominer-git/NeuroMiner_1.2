classdef nk_ItemSelector_App_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        figure1            matlab.ui.Figure
        uipanel1           matlab.ui.container.Panel
        tglOpCases         matlab.ui.control.StateButton
        tglOpFeats         matlab.ui.control.StateButton
        lbCaseFeat         matlab.ui.control.Label
        cmdResetFilt       matlab.ui.control.Button
        cmdFiltCases       matlab.ui.control.Button
        cmdFiltFeats       matlab.ui.control.Button
        edit2              matlab.ui.control.EditField
        edit1              matlab.ui.control.EditField
        tglScale           matlab.ui.control.StateButton
        selMeas            matlab.ui.control.DropDown
        axMat              matlab.ui.control.UIAxes
        axFeats            matlab.ui.control.UIAxes
        axCases            matlab.ui.control.UIAxes
        pnManFeatSel       matlab.ui.container.Panel
        text7              matlab.ui.control.Label
        lstCasesSelected   matlab.ui.control.ListBox
        text6              matlab.ui.control.Label
        cmdAddCasesAll     matlab.ui.control.Button
        cmdRemoveCasesAll  matlab.ui.control.Button
        cmdRemoveCases     matlab.ui.control.Button
        cmdAddCases        matlab.ui.control.Button
        lstCasesPool       matlab.ui.control.ListBox
        cmdAddAll          matlab.ui.control.Button
        cmdRemoveAll       matlab.ui.control.Button
        cmdRemoveItems     matlab.ui.control.Button
        cmdAdd             matlab.ui.control.Button
        text4              matlab.ui.control.Label
        text3              matlab.ui.control.Label
        lstSrc             matlab.ui.control.ListBox
        lstDst             matlab.ui.control.ListBox
        cmdOK              matlab.ui.control.Button
    end

    
    methods (Access = private)
        function figdata = DefineMatrixTextBox(app, handles, cases, list)
            
            nC = numel(cases);
            nL = numel(list);
            
            % Define textbox info data
            pss = cell(nC,nL);
            
            for j = 1:nC
                for i = 1:nL
                    pss{j,i} = sprintf('Case [%g]: %s  |  Feature [%g]: %s', j, cases{j}, i, list{i});
                end
            end
            
            figdata.x = repmat(1:numel(list),nC,1);
            figdata.y = repmat((1:numel(cases))',1,nL);
            figdata.patterntext = pss;
            figdata.lbCaseFeat = handles.lbCaseFeat;
        end
        
        function ModItems(app, cmd, handles, act)
            
            switch cmd
                case {'rem','remall'}
                    switch act
                        case 'feats'
                            AddInd  = 'SrcIndex';
                            RemInd  = 'DstIndex';
                            LstSrc  = 'lstDst';
                            LstDst  = 'lstSrc';
                            T='Items';
                            TT='selFeats';
                            D = 2;
                        case 'cases'
                            AddInd  = 'AllCasesIndex';
                            RemInd  = 'SelectedCasesIndex';
                            LstSrc  = 'lstCasesSelected';
                            LstDst  = 'lstCasesPool';
                            T='Cases';
                            TT='selCases';
                            D = 1;
                    end
            
                case {'add','addall'}
                    switch act
                        case 'feats'
                            AddInd  = 'DstIndex';
                            RemInd  = 'SrcIndex';
                            LstSrc  = 'lstSrc';
                            LstDst  = 'lstDst';
                            T='Items';
                            TT='selFeats';
                            D = 2;
                        case 'cases'
                            RemInd  = 'AllCasesIndex';
                            AddInd  = 'SelectedCasesIndex';
                            LstDst  = 'lstCasesSelected';
                            LstSrc  = 'lstCasesPool';
                            T='Cases';
                            TT='selCases';
                            D = 1;
                    end
            end
            
            if strcmp(handles.(LstSrc).String,''); return; end
            
            if strfind(char(cmd),'all')
                handles.(LstSrc).Value=1:numel(handles.(LstSrc).String);
            end
            
            selItems = get(handles.(LstSrc),'Value');
            
            if ~isempty(selItems)
                switch act
                    case 'feats'
                        handles.(AddInd) = sort([ handles.(AddInd); handles.(RemInd)(selItems) ]);
                    case 'cases'
                        handles.(AddInd) = sort([ handles.(AddInd); handles.(RemInd)(selItems) ]);
                end
            else
                handles.(AddInd) = sort( handles.(AddInd) );
            end
            handles.(RemInd)(selItems)=[];
            if ~numel(handles.(AddInd))
                set(handles.(LstDst),'String',{''});
                set(handles.(LstDst),'Max',1)
                set(handles.(LstDst),'Value',1);
            else
                handles.(LstDst).String = handles.(T)(handles.(AddInd));
                handles.(LstDst).Max = numel(handles.(T)(handles.(AddInd)));
            end
            if ~numel(handles.(RemInd))
                set(handles.(LstSrc),'String',{''});
                set(handles.(LstSrc),'Max',1);
                set(handles.(LstSrc),'Value',1);
            else
                handles.(LstSrc).String = handles.(T)(handles.(RemInd));
                handles.(LstSrc).Max = numel(handles.(T)(handles.(RemInd)));
            end
            nS = numel(handles.(LstSrc).String);
            nD = numel(handles.(LstDst).String);
            if nD > 1,
                handles.(LstDst).Value = numel(handles.(AddInd))-1 ;
            else
                handles.(LstDst).Value = 1;
            end
            if nS > 1,
                handles.(LstSrc).Value = numel(handles.(RemInd))-1 ;
            else
                handles.(LstSrc).Value = 1;
            end
            
            switch act
                case 'feats'
                    handles.(TT) = false(1,size(handles.M, D));
                    handles.(TT)(handles.DstIndex) = true;
                case 'cases'
                    handles.(TT) = false(size(handles.M, D),1);
                    handles.(TT)(handles.SelectedCasesIndex) = true;
            end
            if ~any(handles.(TT))
                cla(handles.axMat);
                cla(handles.axFeats);
                cla(handles.axCases);
            else
                handles = display_matrix(app, handles);
                handles = display_eval(app, handles);
            end
            
            % Update handles structure
            guidata(handles.figure1, handles);
        end
        
        function handles = display_eval(app, handles)
            
            handles.M_plot = filter_matrix(app, handles.M, handles);
            str = handles.selMeas.String{handles.selMeas.Value};
            hold(handles.axFeats,'off');
            hold(handles.axCases,'off');
            
            switch str
                case 'mean'
                    mF = nm_nanmean(handles.M_plot,1);
                    mC = nm_nanmean(handles.M_plot,2);
                case 'median'
                    mF = nm_nanmedian(handles.M_plot,1);
                    mC = nm_nanmedian(handles.M_plot,2);
                case 'variance'
                    mF = nm_nanvar(handles.M_plot,1);
                    mC = nm_nanvar(handles.M_plot,2);
                case 'interquartile range'
                    mF = iqr(handles.M_plot,1);
                    mC = iqr(handles.M_plot,2);
                case 'identical'
                    mC = zeros(size(handles.M_plot,1),1);
                    mF = zeros(size(handles.M_plot,2),1);
                    for i=1:size(handles.M_plot,2)
                        mF(i) = numel(unique(handles.M_plot(:,i)))==1;
                    end
                    for i=1:size(handles.M_plot,1)
                        mC(i) = numel(unique(handles.M_plot(i,:)))==1;
                    end
                case '%(nan)'
                    mF = sum(isnan(handles.M_plot),1)*100/size(handles.M_plot,1);
                    mC = sum(isnan(handles.M_plot),2)*100/size(handles.M_plot,2);
                case 'skewness'
                    mF = skewness(handles.M_plot,1,1);
                    mC = skewness(handles.M_plot,1,2);
                case 'kurtosis'
                    mF = kurtosis(handles.M_plot,1,1);
                    mC = kurtosis(handles.M_plot,1,2);
                case 'cross-correlation'
                    IN.X = handles.M_plot;
                    X = nk_PerfImputeObj(handles.M_plot,IN);
                    nF = size(handles.M_plot,2); mF = zeros(1,nF);
                    for i=1:nF
                        ind = true(1,nF); ind(i)=false;
                        Fi = X(:,i); Fni = X(:,ind);
                        mF(i) = mean(abs(nk_CorrMat(Fni,Fi)));
                    end
                    nC = size(handles.M_plot,1); mC = zeros(1,nC);
                    for i=1:nC
                        ind = true(1,nC); ind(i)=false;
                        Ci = X(i,:); Cni = X(ind,:);
                        mC(i) = mean(abs(nk_CorrMat(Cni',Ci')));
                    end
            end
            
            handles.feats_plot = plot(handles.axFeats,mF,'ko-','LineWidth',1,'MarkerSize', 3);
            if size(handles.M_plot,2)>1,xlim(handles.axFeats,[1 size(handles.M_plot,2)]);end
            handles.cases_plot = plot(handles.axCases,mC,'ko-','LineWidth',1,'MarkerSize', 3);
            if size(handles.M_plot,1)>1,xlim(handles.axCases,[1 size(handles.M_plot,1)]);end
            handles.axCases.XDir = 'reverse';
            hold(handles.axFeats,'on');
            hold(handles.axCases,'on');
            
            if isfield(handles,'thr_feats')
                if isfield(handles,'thrCases'), delete(handles.thrCases);end
                handles.thrFeats = plotpatch(app, handles.M_plot, handles.thr_feats, handles.axFeats, handles.tglOpFeats);
            end
            if isfield(handles,'thr_cases')
                if isfield(handles,'thrCases'), delete(handles.thrCases);end
                handles.thrCases = plotpatch(app, handles.M_plot, handles.thr_cases, handles.axCases, handles.tglOpCases);
            end
            handles.mF = mF;
            handles.mC = mC;
            view(handles.axCases,[-90,90])
            handles.axCases.YDir = 'reverse';
            handles.axCases.XTick = [];
            handles.axFeats.XTick = [];
            
            % Update handles structure
            guidata(handles.figure1, handles)
        end
        
        function handles = display_matrix(app, handles)
            
            if ~isempty(handles.M_plot)
               handles.M_plot = filter_matrix(app, handles.M, handles);
               axes(handles.axMat); cla
               handles.mat_plot = imagesc(handles.M_plot, 'Parent', handles.axMat);
               set(handles.mat_plot,'ButtonDownFcn', {@axMat_ButtonDownFcn,handles});
               %set(handles.mat_plot,'HitTest','off');
               if isfield(handles,'selCases') && ~isempty(handles.selCases), Cases = handles.Cases(handles.selCases); else Cases = handles.Cases; end
               if isfield(handles,'selFeats') &&~isempty(handles.selFeats), Feats = handles.Items(handles.selFeats); else Feats = handles.Items; end
               figdata = DefineMatrixTextBox(app, handles,Cases,Feats);
               set(handles.axMat,'UserData',figdata);
            end
        end
        
        function [sel_ind, thr] = eval_filter(app, thr, ax, lims, vals, sel_prev, op)
            
            sel_ind = true(lims(2)-lims(2)+1,1);
            
            switch op
                case '>='
                    op = @le;
                case '<='
                    op = @ge;
            end
            
            hold(ax,'on');
            try
               sel_ind = op(vals,thr);
               if ~isempty(sel_prev),
                   sel_prev_f = find(sel_prev);
                   sel_new_f = sel_prev_f(sel_ind);
                   sel_ind = false(lims(2)-lims(2)+1,1);
                   sel_ind(sel_new_f)=true;
               end
               if ~any(sel_ind)
                   errordlg('Your filtering was too aggressive. Returning to previous setting')
                   sel_ind = true(lims(2)-lims(2)+1,1);
                   thr=max(vals);
               end
            catch
               errordlg('Invalid input into field!','Error')
            end
        end
        
        function M = filter_matrix(app, M, handles)
            
            if isfield(handles,'selFeats') && ~isempty(handles.selFeats)
                M = M(:,handles.selFeats);
            end
            if isfield(handles,'selCases') && ~isempty(handles.selCases)
                M = M(handles.selCases,:);
            end
        end
        
        function hoverCallback(app, src, evt)
            
                %handles = guidata(gcf);
                %handles=[];
                % Grab the x & y axes coordinate where the mouse is
                axesHdl = get(gcf,'CurrentAxes');
                mousePoint = get(axesHdl, 'CurrentPoint');
                figdata = get(axesHdl,'UserData');
                if ~isempty(figdata)
            
                    mouseX = mousePoint(1,1);
                    mouseY = mousePoint(1,2);
            
                    % Compare where data and current mouse point to find the data point
                    % which is closest to the mouse point
                    distancesToMouse = hypot(figdata.x - mouseX, figdata.y - mouseY);
                    [~, indpat] = min(abs(distancesToMouse(:)));
            
                    % If the distance is less than some threshold, set the text
                    % object's string to show the data at that point.
                    xrange = nk_Range(get(axesHdl, 'Xlim'),2);
                    yrange = nk_Range(get(axesHdl, 'Ylim'),2);
                    if abs(mouseX - figdata.x(indpat)) < 0.02*xrange && abs(mouseY - figdata.y(indpat)) < 0.02*yrange
                        figdata.lbCaseFeat.String =  figdata.patterntext(indpat);
                        figdata.lbCaseFeat.Visible='on';
                    else
                        figdata.lbCaseFeat.Visible='off';
                    end
                    set(axesHdl,'UserData',figdata);
            
                end
        end
        
        function handles = load_list(app, handles, list, act)
            
            switch act
                case 'feats'
                    handles.Items = list;
                    if isfield(handles,'selFeats') && ~isempty(handles.selFeats)
                        selFeats = handles.selFeats;
                    else
                        selFeats = true(numel(list),1);
                    end
                    set(handles.lstDst, 'String', list(selFeats));
                    set(handles.lstDst, 'Max', numel(list(selFeats)),'Min',0);
                    if ~any(selFeats==0)
                        set(handles.lstSrc, 'String', {''});
                        selNFeats = [];
                    else
                        set(handles.lstSrc, 'String', list(~selFeats));
                        selNFeats = ~selFeats;
                    end
                    handles.DstIndex = find(selFeats);
                    handles.SrcIndex = find(selNFeats);
            
                case 'cases'
                    handles.Cases = list;
                    if isfield(handles,'selCases') && ~isempty(handles.selCases)
                        selCases = handles.selCases;
                    else
                        selCases = true(numel(list),1);
                    end
                    set(handles.lstCasesSelected, 'String', list(selCases));
                    set(handles.lstCasesSelected, 'Max', numel(list(selCases)),'Min',0);
                    if ~any(selCases==0)
                        set(handles.lstCasesPool, 'String', {''});
                        selNCases = [];
                    else
                        set(handles.lstCasesPool, 'String', list(~selCases));
                        selNCases = ~selCases;
                    end
                    handles.SelectedCasesIndex = find(selCases);
                    handles.AllCasesIndex = find(selNCases);
            end
        end
        
        function C = orange(app)
            C = [1 .5 0];
        end
        
        function handle = plotpatch(app, M, thr, ax, tgl)
            
            if strcmp(tgl.String,'>=')
                st = thr; hg = ax.YLim(end);
            else
                st = ax.YLim(1); hg = thr;
            end
            handle = patch(ax, [1 size(M,1) size(M,1) 1], [ st st hg hg ], 'red', 'FaceAlpha', 0.5, 'LineStyle', 'none');
        end
        
    end
    

    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function nk_ItemSelector_OpeningFcn(app, varargin)
            % Ensure that the app appears on screen when run
            movegui(app.figure1, 'onscreen');
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app); %#ok<ASGLU>
            
            % This function has no output args, see OutputFcn.
            % hObject    handle to figure
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % varargin   command line arguments to nk_ItemSelector (see VARARGIN)
            
            % Choose default command line output for nk_ItemSelector
            handles.output = hObject;
            
            % Update handles structure
            guidata(hObject, handles);
            handles.M = []; handles.Cases = [];
            if(nargin > 3)
                for index = 1:2:(nargin-3),
                    if nargin-3==index, break, end
                    switch lower(varargin{index})
                        case 'list'
                            list = varargin{index+1};
                        case 'matrix'
                            handles.M = varargin{index+1};
                            handles.YLimCases = [min(min(handles.M,[],2)) max(max(handles.M,[],2))];
                            handles.YLimFeats = [min(min(handles.M)) max(max(handles.M))];
                            handles.M_plot = varargin{index+1};
                            set(handles.axFeats,'XTick',[]);
                            set(handles.axFeats,'YTick',[]);
                        case 'cases'
                            handles.Cases = varargin{index+1};
                        case 'selfeats'
                            handles.selFeats = varargin{index+1};
                        case 'selcases'
                            handles.selCases = varargin{index+1};
                        case 'mode'
                            handles.mode = varargin{index+1};
                    end
                end
            end
            
            switch handles.mode
                case 'feats'
                    handles.edit2.Enable = 'off';
                    handles.cmdFiltCases.Enable = 'off';
                    handles.tglOpCases.Enable = 'off';
                    handles.lstCasesPool.Enable = 'off';
                    handles.lstCasesSelected.Enable = 'off';
                    handles.cmdAddCases.Enable = 'off';
                    handles.cmdRemoveCases.Enable = 'off';
                    handles.cmdAddCasesAll.Enable = 'off';
                    handles.cmdRemoveCasesAll.Enable = 'off';
                case 'cases'
                    handles.edit1.Enable = 'off';
                    handles.cmdFiltFeats.Enable = 'off';
                    handles.tglOpFeats.Enable = 'off';
                    handles.lstSrc.Enable = 'off';
                    handles.lstDst.Enable = 'off';
                    handles.cmdAdd.Enable = 'off';
                    handles.cmdRemoveItems.Enable = 'off';
                    handles.cmdAddAll.Enable = 'off';
                    handles.cmdRemoveAll.Enable = 'off';
            end
            
            handles = load_list(app, handles, list, 'feats');
            handles = load_list(app, handles, handles.Cases, 'cases');
            handles = display_matrix(app, handles);
            handles = display_eval(app, handles);
            
            figdata = DefineMatrixTextBox(app, handles,handles.Cases,handles.Items);
            set(handles.axMat,'UserData',figdata);
            set(gcf, 'WindowButtonMotionFcn', @app.hoverCallback);
            
            % Update handles structure
            guidata(handles.figure1, handles);
            
            % Make the GUI modal
            %set(handles.figure1,'WindowStyle','modal')
            
            % UIWAIT makes nk_OOCVDataIO wait for user response (see UIRESUME)
            uiwait(handles.figure1);
        end

        % Button down function: axMat
        function axMat_ButtonDownFcn(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to axMat (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            figdata = get(handles.axMat,'UserData');
            mousePoint = get(handles.axMat, 'CurrentPoint');
            mouseX = mousePoint(1,1); mouseY = mousePoint(1,2);
            distancesToMouse = hypot(figdata.x - mouseX, figdata.y - mouseY);
            [~, indpat] = min(abs(distancesToMouse(:)));
            xrange = nk_Range(get(handles.axMat, 'Xlim'),2);
            yrange = nk_Range(get(handles.axMat, 'Ylim'),2);
            if abs(mouseX - figdata.x(indpat)) < 0.02*xrange && abs(mouseY - figdata.y(indpat)) < 0.02*yrange
              handles.lstDst.Value = figdata.x(indpat);
            end
        end

        % Button pushed function: cmdAddAll
        function cmdAddAll_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdAddAll (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'addall',handles,'feats')
        end

        % Button pushed function: cmdAddCasesAll
        function cmdAddCasesAll_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdAddCasesAll (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'addall',handles,'cases')
        end

        % Button pushed function: cmdAddCases
        function cmdAddCases_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdAddCases (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'add',handles,'cases')
        end

        % Button pushed function: cmdAdd
        function cmdAdd_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdAdd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'add',handles, 'feats')
        end

        % Button pushed function: cmdFiltCases
        function cmdFiltCases_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdFiltCases (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            if isfield(handles,'thrCases'), delete(handles.thrCases); end
            if isfield(handles,'thr_cases'),
                thr_cases=handles.thr_cases;
            else
                thr_cases = str2double(char(handles.edit2.String));
            end
            if isfield(handles,'selCases'), selCases = handles.selCases; else, selCases = []; end
            [handles.selCases, handles.thr_cases] = eval_filter(app, thr_cases, handles.axCases, [1 size(handles.M,1)], handles.mC, selCases, handles.tglOpCases.String);
            set(handles.edit2,'BackgroundColor','red');
            
            handles = display_matrix(app, handles);
            handles = display_eval(app, handles);
            handles = load_list(app, handles, handles.Items,'feats');
            handles = load_list(app, handles, handles.Cases,'cases');
            
            % Update handles structure
            guidata(handles.figure1, handles)
        end

        % Button pushed function: cmdFiltFeats
        function cmdFiltFeats_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdFiltFeats (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            if isfield(handles,'selFeats'), selFeats = handles.selFeats; else, selFeats = []; end
            if isfield(handles,'thr_feats'),
                thr_feats=handles.thr_feats;
            else
                thr_feats = str2double(char(handles.edit1.String));
            end
            [handles.selFeats, handles.thr_feats] = eval_filter(app, thr_feats, handles.axFeats, [1 size(handles.M,2)], handles.mF, selFeats, handles.tglOpFeats.String);
            set(handles.edit1,'BackgroundColor','red');
            
            handles = display_matrix(app, handles);
            handles = display_eval(app, handles);
            handles = load_list(app, handles, handles.Items,'feats');
            handles = load_list(app, handles, handles.Cases,'cases');
            % Update handles structure
            guidata(handles.figure1, handles)
        end

        % Button pushed function: cmdOK
        function cmdOK_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdOK (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Use UIRESUME instead of delete because the OutputFcn needs
            % to get the updated handles structure.
            
            % Update handles structure
            guidata(hObject, handles);
            
            % Use UIRESUME instead of delete because the OutputFcn needs
            % to get the updated handles structure.
            uiresume(handles.figure1);
        end

        % Button pushed function: cmdRemoveAll
        function cmdRemoveAll_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdRemoveAll (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'remall',handles,'feats')
        end

        % Button pushed function: cmdRemoveCasesAll
        function cmdRemoveCasesAll_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdRemoveCasesAll (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'remall',handles,'cases')
        end

        % Button pushed function: cmdRemoveCases
        function cmdRemoveCases_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdRemoveCases (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'rem',handles,'cases')
        end

        % Button pushed function: cmdRemoveItems
        function cmdRemoveItems_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdRemoveItems (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            ModItems(app, 'rem',handles, 'feats')
        end

        % Button pushed function: cmdResetFilt
        function cmdResetFilt_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cmdResetFilt (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            if isfield(handles,'selFeats'),
                handles = rmfield(handles,'selFeats');
            end
            if isfield(handles,'selCases'),
                handles = rmfield(handles,'selCases');
            end
            if isfield(handles,'thrFeats')
                delete(handles.thrFeats); handles = rmfield(handles,'thrFeats'); handles = rmfield(handles,'thr_feats');
            end
            if isfield(handles,'thrCases')
                delete(handles.thrCases); handles = rmfield(handles,'thrCases'); handles = rmfield(handles,'thr_cases');
            end
            
            set(handles.edit1,'BackgroundColor','white');
            set(handles.edit2,'BackgroundColor','white');
            handles.edit1.String = {'0'};
            handles.edit2.String = {'0'};
            
            handles = display_matrix(app, handles);
            handles = display_eval(app, handles);
            handles = load_list(app, handles, handles.Items,'feats');
            handles = load_list(app, handles, handles.Cases,'cases');
            
            % Update handles structure
            guidata(handles.figure1, handles)
        end

        % Value changed function: edit1
        function edit1_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to edit1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'String') returns contents of edit1 as text
            %        str2double(get(hObject,'String')) returns contents of edit1 as a double
            handles.thr_feats = str2double(char(hObject.String));
            if isfield(handles,'thrFeats'), delete(handles.thrFeats); end
            set(hObject,'BackgroundColor',orange(app));
            handles = display_eval(app, handles);
            
            % Update handles structure
            guidata(handles.figure1, handles)
        end

        % Value changed function: edit2
        function edit2_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to edit2 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'String') returns contents of edit2 as text
            %        str2double(get(hObject,'String')) returns contents of edit2 as a double
            
            handles.thr_cases = str2double(char(hObject.String));
            if isfield(handles,'thrCases'), delete(handles.thrCases); end
            set(hObject,'BackgroundColor',orange(app));
            handles = display_eval(app, handles);
            % Update handles structure
            guidata(handles.figure1, handles)
        end

        % Value changed function: selMeas
        function selMeas_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to selMeas (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: contents = cellstr(get(hObject,'String')) returns selMeas contents as cell array
            %        contents{get(hObject,'Value')} returns selected item from selMeas
            
            handles = display_matrix(app, handles);
            handles = display_eval(app, handles);
            
            % Update handles structure
            guidata(handles.figure1, handles)
        end

        % Value changed function: tglOpCases
        function tglOpCases_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to tglOpFeats (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hint: get(hObject,'Value') returns toggle state of tglOpFeats
            
            if get(hObject,'Value')
               hObject.String = '<=';
            else
               hObject.String = '>=';
            end
            if isfield(handles,'thr_cases')
                if isfield(handles,'thrCases'), delete(handles.thrCases);end
                handles.thrCases = plotpatch(app, handles.M_plot, handles.thr_cases, handles.axCases, handles.tglOpCases);
                % Update handles structure
                guidata(handles.figure1, handles)
            end
        end

        % Value changed function: tglOpFeats
        function tglOpFeats_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to tglOpFeats (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hint: get(hObject,'Value') returns toggle state of tglOpFeats
            
            if get(hObject,'Value')
               hObject.String = '<=';
            else
               hObject.String = '>=';
            end
            if isfield(handles,'thr_feats')
                if isfield(handles,'thrFeats'), delete(handles.thrFeats);end
                handles.thrFeats = plotpatch(app, handles.M_plot, handles.thr_feats, handles.axFeats, handles.tglOpFeats);
                % Update handles structure
                guidata(handles.figure1, handles)
            end
        end

        % Value changed function: tglScale
        function tglScale_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to tglScale (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hint: get(hObject,'Value') returns toggle state of tglScale
            
            if get(hObject,'Value')
               handles.M_unscaled = handles.M;
               handles.M = nk_PerfScaleObj(handles.M,[]);
            else
               handles.M = handles.M_unscaled;
            end
            
            handles = display_matrix(app, handles);
            handles = display_eval(app, handles);
            
            % Update handles structure
            guidata(handles.figure1, handles)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create figure1 and hide until all components are created
            app.figure1 = uifigure('Visible', 'off');
            app.figure1.Colormap = [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0];
            app.figure1.Position = [624 -18 490 843];
            app.figure1.Name = 'NM Matrix Inspector';
            app.figure1.Resize = 'off';
            app.figure1.HandleVisibility = 'callback';
            app.figure1.Tag = 'figure1';

            % Create pnManFeatSel
            app.pnManFeatSel = uipanel(app.figure1);
            app.pnManFeatSel.Tag = 'pnManFeatSel';
            app.pnManFeatSel.FontSize = 8;
            app.pnManFeatSel.Position = [8 15 472 320];

            % Create cmdOK
            app.cmdOK = uibutton(app.pnManFeatSel, 'push');
            app.cmdOK.ButtonPushedFcn = createCallbackFcn(app, @cmdOK_Callback, true);
            app.cmdOK.Tag = 'cmdOK';
            app.cmdOK.BackgroundColor = [0.67843137254902 0.92156862745098 1];
            app.cmdOK.FontSize = 10;
            app.cmdOK.FontWeight = 'bold';
            app.cmdOK.Position = [9.25 9.33333333333333 450.75 32.5];
            app.cmdOK.Text = 'OK';

            % Create lstDst
            app.lstDst = uilistbox(app.pnManFeatSel);
            app.lstDst.Items = {'Destination Items'};
            app.lstDst.Tag = 'lstDst';
            app.lstDst.FontSize = 8;
            app.lstDst.Position = [264.25 184.333333333333 195.75 115.833333333333];
            app.lstDst.Value = 'Destination Items';

            % Create lstSrc
            app.lstSrc = uilistbox(app.pnManFeatSel);
            app.lstSrc.Items = {'Source Items'};
            app.lstSrc.Tag = 'lstSrc';
            app.lstSrc.FontSize = 8;
            app.lstSrc.Position = [9.25 184.333333333333 195.75 115.833333333333];
            app.lstSrc.Value = 'Source Items';

            % Create text3
            app.text3 = uilabel(app.pnManFeatSel);
            app.text3.Tag = 'text3';
            app.text3.HorizontalAlignment = 'center';
            app.text3.VerticalAlignment = 'top';
            app.text3.WordWrap = 'on';
            app.text3.FontSize = 8;
            app.text3.FontWeight = 'bold';
            app.text3.Position = [10 301.833333333333 193.5 14.1666666666667];
            app.text3.Text = 'Feature pool';

            % Create text4
            app.text4 = uilabel(app.pnManFeatSel);
            app.text4.Tag = 'text4';
            app.text4.HorizontalAlignment = 'center';
            app.text4.VerticalAlignment = 'top';
            app.text4.WordWrap = 'on';
            app.text4.FontSize = 8;
            app.text4.FontWeight = 'bold';
            app.text4.Position = [265 301.833333333333 194.25 14.1666666666667];
            app.text4.Text = 'Features to import';

            % Create cmdAdd
            app.cmdAdd = uibutton(app.pnManFeatSel, 'push');
            app.cmdAdd.ButtonPushedFcn = createCallbackFcn(app, @cmdAdd_Callback, true);
            app.cmdAdd.Tag = 'cmdAdd';
            app.cmdAdd.FontSize = 14;
            app.cmdAdd.FontWeight = 'bold';
            app.cmdAdd.Tooltip = 'Add selected features to import';
            app.cmdAdd.Position = [210.8665 274.333333333333 47.25 26.6666666666667];
            app.cmdAdd.Text = '>';

            % Create cmdRemoveItems
            app.cmdRemoveItems = uibutton(app.pnManFeatSel, 'push');
            app.cmdRemoveItems.ButtonPushedFcn = createCallbackFcn(app, @cmdRemoveItems_Callback, true);
            app.cmdRemoveItems.Tag = 'cmdRemoveItems';
            app.cmdRemoveItems.FontSize = 14;
            app.cmdRemoveItems.FontWeight = 'bold';
            app.cmdRemoveItems.Tooltip = 'Remove selected features from import';
            app.cmdRemoveItems.Position = [210.8665 244.333333333333 47.25 26.6666666666667];
            app.cmdRemoveItems.Text = '<';

            % Create cmdRemoveAll
            app.cmdRemoveAll = uibutton(app.pnManFeatSel, 'push');
            app.cmdRemoveAll.ButtonPushedFcn = createCallbackFcn(app, @cmdRemoveAll_Callback, true);
            app.cmdRemoveAll.Tag = 'cmdRemoveAll';
            app.cmdRemoveAll.FontSize = 14;
            app.cmdRemoveAll.FontWeight = 'bold';
            app.cmdRemoveAll.Tooltip = 'Remove all features from import';
            app.cmdRemoveAll.Position = [211 184.333333333333 47.25 26.6666666666667];
            app.cmdRemoveAll.Text = '<<';

            % Create cmdAddAll
            app.cmdAddAll = uibutton(app.pnManFeatSel, 'push');
            app.cmdAddAll.ButtonPushedFcn = createCallbackFcn(app, @cmdAddAll_Callback, true);
            app.cmdAddAll.Tag = 'cmdAddAll';
            app.cmdAddAll.FontSize = 14;
            app.cmdAddAll.FontWeight = 'bold';
            app.cmdAddAll.Tooltip = 'Add all features to import';
            app.cmdAddAll.Position = [210.8665 214.333333333333 47.25 26.6666666666667];
            app.cmdAddAll.Text = '>>';

            % Create lstCasesPool
            app.lstCasesPool = uilistbox(app.pnManFeatSel);
            app.lstCasesPool.Items = {'Source Items'};
            app.lstCasesPool.Tag = 'lstCasesPool';
            app.lstCasesPool.FontSize = 8;
            app.lstCasesPool.Position = [9.25 46.8333333333333 195.75 115.833333333333];
            app.lstCasesPool.Value = 'Source Items';

            % Create cmdAddCases
            app.cmdAddCases = uibutton(app.pnManFeatSel, 'push');
            app.cmdAddCases.ButtonPushedFcn = createCallbackFcn(app, @cmdAddCases_Callback, true);
            app.cmdAddCases.Tag = 'cmdAddCases';
            app.cmdAddCases.FontSize = 14;
            app.cmdAddCases.FontWeight = 'bold';
            app.cmdAddCases.Tooltip = 'Add selected cases to import';
            app.cmdAddCases.Position = [211 136.833333333333 47.25 26.6666666666667];
            app.cmdAddCases.Text = '>';

            % Create cmdRemoveCases
            app.cmdRemoveCases = uibutton(app.pnManFeatSel, 'push');
            app.cmdRemoveCases.ButtonPushedFcn = createCallbackFcn(app, @cmdRemoveCases_Callback, true);
            app.cmdRemoveCases.Tag = 'cmdRemoveCases';
            app.cmdRemoveCases.FontSize = 14;
            app.cmdRemoveCases.FontWeight = 'bold';
            app.cmdRemoveCases.Tooltip = 'Remove selected cases from import';
            app.cmdRemoveCases.Position = [211 106.833333333333 47.25 26.6666666666667];
            app.cmdRemoveCases.Text = '<';

            % Create cmdRemoveCasesAll
            app.cmdRemoveCasesAll = uibutton(app.pnManFeatSel, 'push');
            app.cmdRemoveCasesAll.ButtonPushedFcn = createCallbackFcn(app, @cmdRemoveCasesAll_Callback, true);
            app.cmdRemoveCasesAll.Tag = 'cmdRemoveCasesAll';
            app.cmdRemoveCasesAll.FontSize = 14;
            app.cmdRemoveCasesAll.FontWeight = 'bold';
            app.cmdRemoveCasesAll.Tooltip = 'Remove all cases from import';
            app.cmdRemoveCasesAll.Position = [211 46.8333333333333 47.25 26.6666666666667];
            app.cmdRemoveCasesAll.Text = '<<';

            % Create cmdAddCasesAll
            app.cmdAddCasesAll = uibutton(app.pnManFeatSel, 'push');
            app.cmdAddCasesAll.ButtonPushedFcn = createCallbackFcn(app, @cmdAddCasesAll_Callback, true);
            app.cmdAddCasesAll.Tag = 'cmdAddCasesAll';
            app.cmdAddCasesAll.FontSize = 14;
            app.cmdAddCasesAll.FontWeight = 'bold';
            app.cmdAddCasesAll.Tooltip = 'Add all cases to import';
            app.cmdAddCasesAll.Position = [211 76.8333333333333 47.25 26.6666666666667];
            app.cmdAddCasesAll.Text = '>>';

            % Create text6
            app.text6 = uilabel(app.pnManFeatSel);
            app.text6.Tag = 'text6';
            app.text6.HorizontalAlignment = 'center';
            app.text6.VerticalAlignment = 'top';
            app.text6.WordWrap = 'on';
            app.text6.FontSize = 8;
            app.text6.FontWeight = 'bold';
            app.text6.Position = [7 166 193.5 14.1666666666667];
            app.text6.Text = 'Cases pool';

            % Create lstCasesSelected
            app.lstCasesSelected = uilistbox(app.pnManFeatSel);
            app.lstCasesSelected.Items = {'Destination Items'};
            app.lstCasesSelected.Tag = 'lstCasesSelected';
            app.lstCasesSelected.FontSize = 8;
            app.lstCasesSelected.Position = [264.25 46.8333333333333 195.75 115.833333333333];
            app.lstCasesSelected.Value = 'Destination Items';

            % Create text7
            app.text7 = uilabel(app.pnManFeatSel);
            app.text7.Tag = 'text7';
            app.text7.HorizontalAlignment = 'center';
            app.text7.VerticalAlignment = 'top';
            app.text7.WordWrap = 'on';
            app.text7.FontSize = 8;
            app.text7.FontWeight = 'bold';
            app.text7.Position = [265 166 193.5 14.1666666666667];
            app.text7.Text = 'Cases to be imported';

            % Create uipanel1
            app.uipanel1 = uipanel(app.figure1);
            app.uipanel1.TitlePosition = 'centertop';
            app.uipanel1.Title = 'Matrix Inspector';
            app.uipanel1.Tag = 'uipanel1';
            app.uipanel1.FontWeight = 'bold';
            app.uipanel1.FontSize = 8;
            app.uipanel1.Position = [8 342 474 497];

            % Create axCases
            app.axCases = uiaxes(app.uipanel1);
            app.axCases.Colormap = [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0];
            app.axCases.FontSize = 13.968253968254;
            app.axCases.NextPlot = 'replace';
            app.axCases.Tag = 'axCases';
            app.axCases.Position = [330 23 128 294];

            % Create axFeats
            app.axFeats = uiaxes(app.uipanel1);
            app.axFeats.Colormap = [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0];
            app.axFeats.FontSize = 13.940329218107;
            app.axFeats.NextPlot = 'replace';
            app.axFeats.Tag = 'axFeats';
            app.axFeats.Position = [4 306 346 165];

            % Create axMat
            app.axMat = uiaxes(app.uipanel1);
            app.axMat.Colormap = [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0];
            app.axMat.FontSize = 13.9596602972399;
            app.axMat.NextPlot = 'replace';
            app.axMat.ButtonDownFcn = createCallbackFcn(app, @axMat_ButtonDownFcn, true);
            app.axMat.Tag = 'axMat';
            app.axMat.Position = [4 24 346 293];

            % Create selMeas
            app.selMeas = uidropdown(app.uipanel1);
            app.selMeas.Items = {'mean', 'median', 'variance', 'interquartile range', 'identical', '%(nan)', 'kurtosis', 'skewness', 'cross-correlation'};
            app.selMeas.ValueChangedFcn = createCallbackFcn(app, @selMeas_Callback, true);
            app.selMeas.Tag = 'selMeas';
            app.selMeas.FontSize = 8;
            app.selMeas.BackgroundColor = [1 1 1];
            app.selMeas.Position = [355.913112164297 434.147512864494 95.0394944707741 29.2667238421955];
            app.selMeas.Value = 'mean';

            % Create tglScale
            app.tglScale = uibutton(app.uipanel1, 'state');
            app.tglScale.ValueChangedFcn = createCallbackFcn(app, @tglScale_Callback, true);
            app.tglScale.Tag = 'tglScale';
            app.tglScale.Text = 'Scale';
            app.tglScale.FontSize = 8.91938250428815;
            app.tglScale.Position = [355.913112164297 415.7512864494 95.781990521327 23.4133790737564];

            % Create edit1
            app.edit1 = uieditfield(app.uipanel1, 'text');
            app.edit1.ValueChangedFcn = createCallbackFcn(app, @edit1_Callback, true);
            app.edit1.Tag = 'edit1';
            app.edit1.HorizontalAlignment = 'center';
            app.edit1.FontSize = 8.83722060459171;
            app.edit1.Position = [386.355450236967 385.648370497427 34.8973143759873 23.4133790737564];
            app.edit1.Value = '0';

            % Create edit2
            app.edit2 = uieditfield(app.uipanel1, 'text');
            app.edit2.ValueChangedFcn = createCallbackFcn(app, @edit2_Callback, true);
            app.edit2.Tag = 'edit2';
            app.edit2.HorizontalAlignment = 'center';
            app.edit2.FontSize = 8.83722060459171;
            app.edit2.Position = [386.355450236967 357.217838765009 34.8973143759873 23.4133790737564];
            app.edit2.Value = '0';

            % Create cmdFiltFeats
            app.cmdFiltFeats = uibutton(app.uipanel1, 'push');
            app.cmdFiltFeats.ButtonPushedFcn = createCallbackFcn(app, @cmdFiltFeats_Callback, true);
            app.cmdFiltFeats.Tag = 'cmdFiltFeats';
            app.cmdFiltFeats.FontSize = 11.7066895368782;
            app.cmdFiltFeats.Tooltip = 'Filter Features';
            app.cmdFiltFeats.Position = [425.707740916272 385.648370497427 25.9873617693523 23.4133790737564];
            app.cmdFiltFeats.Text = 'Ff';

            % Create cmdFiltCases
            app.cmdFiltCases = uibutton(app.uipanel1, 'push');
            app.cmdFiltCases.ButtonPushedFcn = createCallbackFcn(app, @cmdFiltCases_Callback, true);
            app.cmdFiltCases.Tag = 'cmdFiltCases';
            app.cmdFiltCases.FontSize = 11.7066895368782;
            app.cmdFiltCases.Tooltip = 'Filter Cases';
            app.cmdFiltCases.Position = [425.707740916272 357.217838765009 26.7298578199052 23.4133790737564];
            app.cmdFiltCases.Text = 'Fc';

            % Create cmdResetFilt
            app.cmdResetFilt = uibutton(app.uipanel1, 'push');
            app.cmdResetFilt.ButtonPushedFcn = createCallbackFcn(app, @cmdResetFilt_Callback, true);
            app.cmdResetFilt.Tag = 'cmdResetFilt';
            app.cmdResetFilt.FontSize = 8.83722060459171;
            app.cmdResetFilt.Position = [355.913112164297 327.951114922813 95.781990521327 23.4133790737564];
            app.cmdResetFilt.Text = 'Reset Filter';

            % Create lbCaseFeat
            app.lbCaseFeat = uilabel(app.uipanel1);
            app.lbCaseFeat.Tag = 'lbCaseFeat';
            app.lbCaseFeat.HorizontalAlignment = 'center';
            app.lbCaseFeat.VerticalAlignment = 'top';
            app.lbCaseFeat.WordWrap = 'on';
            app.lbCaseFeat.FontSize = 8;
            app.lbCaseFeat.Position = [29.4848484848485 1.82177014390404 424.274322169059 25.1288659793814];
            app.lbCaseFeat.Text = '';

            % Create tglOpFeats
            app.tglOpFeats = uibutton(app.uipanel1, 'state');
            app.tglOpFeats.ValueChangedFcn = createCallbackFcn(app, @tglOpFeats_Callback, true);
            app.tglOpFeats.Tag = 'tglOpFeats';
            app.tglOpFeats.Text = '>=';
            app.tglOpFeats.FontSize = 11.7066895368782;
            app.tglOpFeats.FontWeight = 'bold';
            app.tglOpFeats.Position = [355.913112164297 385.648370497427 25.2448657187994 23.4133790737564];

            % Create tglOpCases
            app.tglOpCases = uibutton(app.uipanel1, 'state');
            app.tglOpCases.ValueChangedFcn = createCallbackFcn(app, @tglOpCases_Callback, true);
            app.tglOpCases.Tag = 'tglOpCases';
            app.tglOpCases.Text = '>=';
            app.tglOpCases.FontSize = 11.7066895368782;
            app.tglOpCases.FontWeight = 'bold';
            app.tglOpCases.Position = [355.913112164297 357.217838765008 25.2448657187994 23.4133790737564];

            % Show the figure after all components are created
            app.figure1.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = nk_ItemSelector_App_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.figure1)

                % Execute the startup function
                runStartupFcn(app, @(app)nk_ItemSelector_OpeningFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.figure1)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.figure1)
        end
    end
end