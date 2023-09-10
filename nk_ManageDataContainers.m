function [act, datacontainer, inp ] = nk_ManageDataContainers(act, inp, datacontainer, parentstr)
nk_PrintLogo
if inp.oocvmode 
    containerstr = 'OOCV ';
else
    containerstr = 'Training/CV ';
end

% check whether an alternative label was used in one of the locked analyses
altlabels = [];
if isfield(inp,'analysis')
    for i=1:length(inp.analysis)
        if isfield(inp.analysis{1,i}.params, 'label') && inp.analysis{1,i}.params.label.altlabelflag
            altlabels = [altlabels, inp.analysis{1,i}.params.label.labelname];
        end
    end
end

navistr = sprintf('%s\t>>> %s',parentstr, [containerstr 'data input']); fprintf('\nYou are here: %s >>>',navistr); 
if isempty(act), act = 1; end
if act < 10
    fprintf('\n\n\t============================')
    fprintf('\n\t    AVAILABLE MODALITIES  ') 
    fprintf('\n\t============================')
    fprintf('\n')
    for i=1:inp.nummodal
        availstr = 'available';

        if ~isfield(datacontainer,'Y')  || i > numel(datacontainer.Y) || (iscell(datacontainer.Y) && isempty(datacontainer.Y{i}))
            availstr = 'not loaded';
        elseif isfield(datacontainer,'Y') 
            if iscell(datacontainer.Y) 
                if ischar(datacontainer.Y{i})
                    containerpth = datacontainer.Y{i};
                    availstr = 'linked';
                elseif isnumeric(datacontainer.Y{i})
                    availstr = 'loaded';
                else
                    availstr = 'not loaded';
                end
            else 
                if ischar(datacontainer.Y) 
                    containerpth = datacontainer.Y;
                    availstr = 'linked';
                elseif isnumeric(datacontainer.Y)
                    availstr = 'loaded';
                else
                    availstr = 'not loaded';
                end
            end
        end
        str = sprintf('Modality %g [ %s ]: ', i, inp.datadescriptor{i}.desc);
        fprintf('\n'); 
        if i==inp.currmodal
            fprintf('==> %s ', str); 
        else
            fprintf('\t%s', str); 
        end
        if strcmp(availstr,'linked')
            fprintf('%s [ %s ] ', availstr, containerpth);
        else
            fprintf( '%s ', availstr);
        end   
    end

    fprintf('\n')
    if inp.nummodal > 1
        mnuact = [ sprintf('Select modality [ currently: modality %g ]',inp.currmodal) ...
               sprintf('|Add/Modify %sdata in modality %g', containerstr, inp.currmodal) ... 
               sprintf('|Delete %sdata in modality %g', containerstr, inp.currmodal) ]; 
    
        mnusel = 1:3;
    else
        mnuact = ['Add/Modify ' containerstr 'data|Delete ' containerstr 'data' ]; 
        mnusel = 2:3;
    end
    
    switch availstr
        case 'loaded'
            mnuact = [mnuact ['|Export ' containerstr 'modalities to files and create links in NM structure|Replace ' containerstr 'modalities with links to files']];
            mnusel = [mnusel  4 5];
        case 'linked'
            mnuact = [mnuact ['|Update ' containerstr 'data link ']];
            mnusel = [mnusel  5];
            if exist(datacontainer.Y,'file')
                mnuact = [mnuact ['|Re-import ' containerstr 'modalities from file into NM structure']];
                mnusel = [mnusel  6];
            end
        case 'not loaded'
            mnuact = [mnuact ['|Fill ' containerstr 'modalities with links to files']];
            mnusel = [mnusel 5];
    end
    
    if inp.oocvmode
    
        if inp.covflag && isfield(datacontainer,'Y') && datacontainer.n_subjects_all>0
            if isfield(datacontainer,'covars')
                mnuact = [ mnuact '|Modify covariates in OOCV data' ];    
            else
                mnuact = [ mnuact '|Add covariate data in OOCV data' ];   
            end
            mnusel = [mnusel 7];
        end

        if ~isempty(altlabels)
            mnuact = [mnuact '|Add alternative label(s) to OOCV data'];
            mnusel = [mnusel 8]; 
        end

        if isfield(datacontainer,'groups')
            mnuact = [mnuact '|Modify subgroup information to OOCV data'];
            mnuact = [mnuact '|Remove subgroup information from OOCV data'];
            mnusel = [mnusel 9 10];
        else
            mnuact = [mnuact '|Add subgroup information to OOCV data'];
            mnusel = [mnusel 9];
        end

    end
    
    act = nk_input(sprintf('Select action for %scontainer %s', containerstr, inp.desc),0,'mq',mnuact,mnusel);

end

switch act
    
    case {1, 11} % Define active modality
        inp.currmodal = nk_input('Select active modality',0,'i',inp.currmodal);
        if inp.currmodal > inp.nummodal, inp.currmodal = inp.nummodal; end
    case {2, 12} % Read-in data
        datacontainer = InputDataModality(inp, datacontainer, navistr);
    case {3, 13} % Clear data from modality
        datacontainer = ClearDataModality(datacontainer);
        if ~isfield(datacontainer,'Y'), act = 0; end
    case {4 , 14}
        datacontainer = LinkData2Disk(inp, datacontainer, 'export&link');
    case {5 , 15}
        datacontainer = LinkData2Disk(inp, datacontainer, 'replacelink');
    case {6, 16}
        datacontainer = ReimportDatafromDisk(datacontainer);
    case {7, 17}
        % Don't forget the covariates if they are present in the discovery data
        datacontainer.covars = nk_DefineCovars_config(datacontainer.n_subjects_all, inp.covars); 
    case {8, 18}
        % if alternative labels were used in any of the locked analyses,
        % new labels have to be input for the validation data too 
        datacontainer.label = nk_DataLabel_config(datacontainer.n_subjects_all, altlabels);
    case {9, 19}
        if isfield(datacontainer,'groups')
            groups = datacontainer.groups;
        else
            groups = [];
        end
        datacontainer.groups = nk_input('Specify logical subgroup matrix (Each column indicates one subgroup)', 0, 'e', [], [datacontainer.n_subjects_all, inf]);
        if isfield(datacontainer,'grpnames') && numel(datacontainer.grpnames) == size(datacontainer.groups,2)
            grpnames = datacontainer.grpnames;
        else
            grpnames = [];
        end
        datacontainer.grpnames = nk_input('Provide subgroup names', 0, 'e', [], size(datacontainer.groups,2));
        if strcmp(inp.modeflag,'classification')
            if isfield(datacontainer,'refgroup')
                if ~datacontainer.refgroup, refgroup = 2; else, refgroup = 1; end
            else
                refgroup = 2;
            end
            refgroupflag = nk_input('Do you want to specify a reference group among the subgroups',0,'yes|no',[1 0], refgroup);
            if refgroupflag
                datacontainer.refgroup = nk_input('Select reference subgroup', 0, 'm', strjoin(datacontainer.grpnames,'|'), 1:numel(datacontainer.grpnames));  
            end
        end
    case {10, 20}
        datacontainer = rmfield(datacontainer,'groups');
        datacontainer = rmfield(datacontainer,'grpnames');
        if isfield(datacontainer,'refgroup')
            datacontainer = rmfield(datacontainer,'refgroup');
        end
end

% _________________________________________________________________________
function datacontainer = ReimportDatafromDisk(datacontainer)
fprintf('\nRe-importing linked data container into NM: %s',datacontainer.Y);
if iscell(datacontainer.Y) 
    for i=1:numel(datacontainer.Y)
        if ischar(datacontainer.Y{i}) && exist(datacontainer.Y{i},"file")
            datafile = datacontainer.Y{i};
            datacontainer.Y{i} = load(datafile,'Y');
        end
    end
else
    load(datacontainer.Y);
    if exist("OOCV","var")
       datacontainer.Y = OOCV;
    else
       datacontainer.Y = Y;
    end
end

% _________________________________________________________________________
function datacontainer = LinkData2Disk(inp, datacontainer, act)
global OCTAVE

if inp.oocvmode
    prefix = ['OOCV_NM' inp.id '_'];
else
    prefix = ['TrCV_NM' inp.id '_'];
end

switch act
    case 'export&link'
        [filename, pathname] = uiputfile({'.mat'},'Save data container to disk');
        if strcmp(filename,'.mat'), filename = '_.mat'; end
        if isempty(filename), return, end
        for i=1:numel(datacontainer.Y)
            pth = regexprep(filename,'.mat', sprintf('M%g.mat',i));
            pth = fullfile(pathname, [prefix pth] );
            Y = datacontainer.Y{i};
            Modality = i;
            fprintf('\nExporting Modality %g to %s', Modality, pth)
            try
                save(pth, 'Y', 'Modality');
            catch
                if OCTAVE
                    save(pth, 'Y', 'Modality');
                else
                    save(pth, 'Y', 'Modality', '-v7.3');
                end
            end
            datacontainer.Y{i} = pth;
        end
    case 'replacelink'
        [filename, pathname] = uigetfile({'.mat'},'Link data container to data file on disk','MultiSelect','on');
        if iscell(filename)
            for i=1:numel(filename)
                suff = sprintf('_M%g',i);
                idx = strcmp(filename,suff);
                if ~any(idx)
                    errordlg(sprintf('Your file selection does not contain a valid NM Data Container for Modality #%g',i),'Error')
                    break
                end
                pth = fullfile(pathname, filename{idx});
                datacontainer.Y{i} = pth;
            end
        else
            pth = fullfile(pathname, filename);
            load(pth)
            if exist("OOCV","var")
                datacontainer = OOCV;
                datacontainer.Y = pth;
            elseif exist("Modality","var")
                datacontainer.Y{Modality} = pth;
            else
                errordlg(sprintf('Cannot assign file %s to Modality', pth),'Error')
            end
        end
end

% _________________________________________________________________________
function datacontainer = InputDataModality(inp, datacontainer, parentstr)

nk_PrintLogo
fprintf('\n\n'); mestr = sprintf('Input %s for Modality %g', inp.dattype, inp.currmodal);  
navistr = sprintf('%s >>> %s',parentstr, mestr); fprintf('\nYou are here: %s >>>',navistr); 

% Retrieve input settings from the discovery data
IO = inp.datadescriptor{inp.currmodal}.input_settings;
if isfield(IO,'selCases'), IO = rmfield(IO,'selCases'); end

% Remove setting for non-labeled subjects
IO.nangroup=false;
IO.nan_subjects=0;
if isfield(IO,'Pnan')
    IO = rmfield(IO,'Pnan');
    IO = rmfield(IO,'Vnan');
end

% Activate independent test data input
IO.oocvflag = true;
IO.labels_known = datacontainer.labels_known;
IO.badcoords = inp.badcoords{inp.currmodal}; 
IO.brainmask = inp.brainmask{inp.currmodal}; 
IO.Ydims = inp.Ydims(inp.currmodal);

if IO.labels_known
    IO.n_subjects = IO.n_subjects/0;
    IO.n_subjects_all = Inf;
else
    IO.n_subjects = Inf;
    IO.n_subjects_all = Inf;
    IO.n_samples = 1;
end
    
if inp.currmodal>1 && isfield(datacontainer,'cases')
    IO.ID = datacontainer.cases;
else
    IO = rmfield(IO,'ID');
    if isfield(IO,'survanal_time'), IO = rmfield(IO,'survanal_time'); end
end
    
if strcmp(IO.datasource,'matrix')
    IO.matrix_edit = inp.na_str;
    IO.sheets = inp.na_str;
    IO.sheet = inp.na_str;
    IO.sheets = inp.na_str;
    IO.M_edit = inp.na_str;
    IO.featnames_cv = inp.featnames{inp.currmodal};
else
    if strcmp(IO.datasource,'spm')
        IO.datasource = 'nifti'; 
        IO.groupmode = 1;
        IO = rmfield(IO,'design');
        IO = SetFileFilter(IO,IO.groupmode,IO.datasource);
    end
    IO.globvar_edit = inp.na_str;
    if isfield(IO,'g') && ~isempty(IO.g), IO = rmfield(IO,'g'); end
    if IO.labels_known
        IO.P = repmat({[]},1,IO.n_samples);
        IO.V = IO.P;
    else
        IO.P = []; 
        IO.V = [];
    end
    IO.PP = [];
    IO = rmfield(IO,'Vinfo');
    IO = rmfield(IO,'Vvox');
    IO = rmfield(IO,'F');
    IO = rmfield(IO,'files');
    if isfield(IO,'L') && ~isempty(IO.L)
        IO.label_edit = inp.na_str;
        IO = rmfield(IO,'L'); 
    end
end
t_act = Inf; t_mess = [];while ~strcmp(t_act,'BACK'), [ IO, t_act, t_mess ] = DataIO( datacontainer , mestr, IO, t_mess, inp.currmodal);  end
if IO.completed
    datacontainer = TransferModality2NM( datacontainer, IO, inp.currmodal ); 
    datacontainer.n_subjects_all = size(datacontainer.cases,1);
end
% _________________________________________________________________________
function datacontainer = ClearDataModality(inp, datacontainer)

if iscell(datacontainer.Y)
    datacontainer.Y{inp.currmodal} = [];
else
    datacontainer.Y = [];
end
datacontainer.brainmask{inp.currmodal} = [];
datacontainer.badcoords{inp.currmodal} = [];
datacontainer.datadescriptor{inp.currmodal} = [];
datacontainer.files{inp.currmodal} = [];
datacontainer.featnames{inp.currmodal} = [];

if iscell(datacontainer.Y) &&  ~sum(~cellfun(@isempty,datacontainer.Y ))
    datacontainer = rmfield(datacontainer,'Y');
    datacontainer = rmfield(datacontainer,'brainmask');
    datacontainer = rmfield(datacontainer,'badcoords');
    datacontainer = rmfield(datacontainer,'datadescriptor');
    datacontainer = rmfield(datacontainer,'files');
    datacontainer = rmfield(datacontainer,'featnames') ;
    datacontainer = rmfield(datacontainer,'cases');
    if isfield(datacontainer,'covars')
        datacontainer = rmfield(datacontainer,'covars');
    end
    if isfield(datacontainer,'groupnames')
        datacontainer = rmfield(datacontainer,'groupnames');
    end
    datacontainer = rmfield(datacontainer,'label');  
    datacontainer.n_subjects = 0;
    datacontainer.n_subjects_all = 0;
end

