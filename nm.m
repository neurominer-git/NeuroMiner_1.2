function nm(varargin)
% 
% NeuroMiner Startup and Main Interface Function
%
% *******************************************
% ****\                                 /****
% *****\     ~~~~~~~~~~~~~~~~~~~~~     /*****
% ******\     N e u r o M i n e r     /******
% *******\   ~~~~~~~~~~~~~~~~~~~~~   /*******
% *******/                           \*******
% ******/     pattern recognition     \******
% *****/      for neurodiagnostic      \*****
% ****/           applications          \****
% ***/                                   \***
% *******************************************
%
% type 'nm <varargin>' at the MATLAB command line to start the program
% <varargin> could be e.g.
% nm nosplash 
% nm nosplash expert
% nm expert
% =========================================================================
% nm.m implements three different analysis stages;
%
% (1) Data import: Allows the user to add, modify and delete modalities from
%     the NM workspace. Data can be input into NM from a variety of
%     sources. Currently, the user can choose among 3D neuroimaging data 
%     (.img/.nii, .nii.gz/.img.gz and .mgh/.mgz formats are supported.),
%     MATLAB 2D-matrix data, and table-based /text-based formats 
%     (.xlsx/.xls,.txt/.csv). For latter data formats a matrix inspector
%     allows the user to visualize the matrix, to select features and cases
%     based on descriptive statistics criteria.
%
% (2) Model generation and cross-validation: Enables the user to define a
%     hyperparameter space for the optimization of complex machine learning
%     pipelines. The model generation and validation process is split into
%     three steps - preprocessing, model training and visualization. These
%     steps are completely wrapped into a repeated, nested cross-validation 
%     structure that is defined by the user. In addition, the user can run 
%     NM's intepretable ML module to compute the input features' saliencies 
%     for the individual predictions obtained at the CV2 level. NM currently 
%     provides classification and regression models. Classification can be 
%     conducted for binary and multi-group problems, with the latter being 
%     always decomposed into one-vs-one or one-vs-all binary classifiers.
%     Currently, a real multi-group learning process is not implemented, 
%     but it is planned for a future release. Results can be inspected with
%     the NM Results Viewer.
%     
% (3) Application to independent data: Allows the user to apply models
%     trained in (2) to new datasets, which have to comply with the input 
%     settings defined in (1). The entire set of models obtained
%     from the cross-validation setup defined in (1) and trained in (2)
%     will be applied to these new data. In addition, the user can use NM's
%     interpretation module to gain insight into which features determined 
%     the individual prediction most.
%
% Each stage has to be completed in order to move to the next stage
% User can choose to use NM in expert mode by invoking nm with the 'expert'
% option.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2022
global EXPERT DEV OCTAVE NM 

nosplash  = false; EXPERT = false; OCTAVE = false; DEV = false; 
% Initialize NM
% Show splash
if exist('varargin','var')
    for i=1:numel(varargin)
        switch varargin{i}
            case 'nosplash'
                nosplash = true; 
            case 'expert'
                EXPERT = true;
            case 'octave'
                OCTAVE = true;
            case 'dev'
                DEV = true;
        end
    end
end
if ~nosplash, splash('splash','png',1500); end
obj = onCleanup(@()QuitNeuroMiner());
init = struct('all', 1, 'addrootpath', 0, 'addDRpath', 0, 'addLIBSVMpath', 0, 'addMIpath', 0, 'addLIBLINpath', 0, 'addMikeRVMpath', 0);
 nk_Initialize(init);

% Loop into the main interface as long as the user does not decide to quit
% or an error occurs.
act = []; ERR=[];
while ~strcmp(act,'quit') && isempty(ERR)
    [ERR, act] = nm_interface; 
end

% Perform error handling if needed
if ~isempty(ERR)
    nk_SetupGlobVars2(NM,'clear')
    fprintf('\n========================================================== \n');
    fprintf('>>> OOOops... NM crashed :-(( <<<');
    fprintf('\nError: %s\nStack:', ERR.message);
    for i=1:numel(ERR.stack)
        fprintf('\nIn %s, name: %s, line: %g.', ...
            ERR.stack(i).file,ERR.stack(i).name,ERR.stack(i).line)
        switch ERR.stack(i).name
            case {'nk_PreprocessPrep','nk_MLOptimizerPrep','nk_VisModelsPrep','nk_OOCVPrep'}
                if isfield(NM,'runtime')
                    nk_NMLogFileManager('add_entry', NM, NM.analysis{NM.runtime.curanal}, 'nm', ERR);
                end
        end
    end
    fprintf('\n>>> NM structure secured to MATLAB workspace.')
    fprintf('\n\nPlease copy the above message and email it with a description to: ')
    fprintf('\nnikolaos.koutsouleris@med.uni-muenchen.de\n')
    if isfield(NM,'runtime'), NM=rmfield(NM,'runtime'); end
end

end

function [ERR, action] = nm_interface
global NM BATCH NMinfo
ERR = [];

if ~isempty(NM) && ~isstruct(NM) 
    error('The NM workspace variable does not have the correct format and thus cannot be not recognized. Clear the variable and start again!')
else
    NMfields = fieldnames(NM); action = 'quit';
end

try
    
    menutitle = 'MAIN INTERFACE';
    
    % Check NM status and define current NM mode based on status
    [s, paramstr] = nk_GetNMStatus(NM);
    
    if ~s.import_finished
        NMinfo.clback = rgb('HoneyDew');
        NMinfo.cllogo = rgb('ForestGreen');
        NMinfo.clmenu = rgb('DarkGreen');
        mn_str = 'Load data for model discovery and cross-validation|'; mn_act = 1;
        menutitle = [menutitle ' [ DATA INPUT MODE ]'];
    else
        switch NM.modeflag
            case 'classification'
                mdltypestr = 'classifiers';
            case 'regression'
                mdltypestr = 'regressors';
        end
        
        switch s.setup_ok
            
            case 1 % Parameter definitions ok
                
                if ~s.analyses_locked
                    
                    NMinfo.clback = rgb('LightCyan');
                    NMinfo.cllogo = rgb('SteelBlue');
                    NMinfo.clmenu = rgb('DarkBlue');
                    menutitle = [menutitle ' [ MODEL DISCOVERY MODE ]'];
                    mn_str = 'Inspect data used for model discovery and cross-validation|'; mn_act = 8;
                    
                    mn_str = [ mn_str 'Set up NM parameter workspace|' ]; mn_act = [ mn_act 2 ];
                    
                    if s.analyses_exist 
                        mn_str = [ mn_str ...
                            'Initialize & manage analyses|' ...
                            'Preprocess features|' ...
                            'Train supervised ' mdltypestr '|'];
                        mn_act = [mn_act 3:5];
                    else
                        mn_str = [ mn_str 'Initialize analyses|'];
                        mn_act = [mn_act 3];
                    end
                    
                    if s.analyses_ready
                        if s.analyses_nondeterministic
                             mn_str = [mn_str ...
                                'Interpret predictions of ' mdltypestr ' in discovery data|' ...
                                'Open NM Results Viewer (cross-validation results)|'];
                            mn_act = [mn_act 18 7];
                        else
                            mn_str = [mn_str ...
                                'Visualize ' mdltypestr '|' ...
                                'Interpret predictions of ' mdltypestr ' in discovery data|' ...
                                'Open NM Results Viewer (cross-validation results)|'];
                            mn_act = [mn_act 6 18 7];
                        end
                    end
                    
                    if s.analyses_completed 
                       mn_str = [ mn_str  'Lock analyses and start NM application mode|' ]; mn_act = [ mn_act 9 ];
                    end

                else
                    NMinfo.clback = rgb('Linen');
                    NMinfo.cllogo = rgb('DarkSalmon');
                    NMinfo.clmenu = rgb('IndianRed');
                    mn_str = 'Update analyses'' root paths'; mn_act = 17;
                    mn_str = [ mn_str '|Load data for model application' ]; mn_act = [mn_act 1 ];
                    if s.oocv_data_ready
                        mn_str = [ mn_str '|Set up parameters for model application' ]; mn_act = [ mn_act 2 ];
                        if s.oocv_anal_ready
                            mn_str = [ mn_str ...
                                '|Apply ' mdltypestr ' to independent data' ...
                                '|Interpret predictions of ' mdltypestr ' in independent data']; 
                            mn_act = [ mn_act 10 18]; 
                        end
                    end
                    if s.oocv_anal_ready
                        mn_str = [mn_str '|Open NM Results Viewer (cross-validation & independent test results)']; 
                    else
                        mn_str = [mn_str '|Open NM Results Viewer (cross-validation results)']; 
                    end
                    mn_act = [mn_act 7];
                    if ~isfield(NM.defs,'data_scrambled') || ~NM.defs.data_scrambled
                        menutitle = [menutitle ' [ MODEL APPLICATION MODE ]'];
                        mn_str = [mn_str '|Shred input data in NM structure (external validation)']; mn_act = [ mn_act 16 ];
                    else
                        menutitle = [menutitle ' [ MODEL EXPORT MODE ]'];
                    end
                end

            case 0
                mn_str = ['Inspect data used for model discovery and cross-validation' ...
                            '|Set up NM parameter workspace' ];
                mn_act = [8 2];
                NMinfo.clback = rgb('LightSalmon');
                NMinfo.cllogo = rgb('DarkRed');
                NMinfo.clmenu = rgb('Crimson');
        end

    end
    
    NM = nk_InitNMwindowColors(NM,NMinfo.clback);
    
    nk_PrintLogo(true)
    
    fprintf('\nCurrent working directory: %s',pwd)
    if ~s.setup_ok && s.import_finished
        fprintf('\n\n')
        cprintf('*red','Parameter setup not complete! \n')
        if iscell(paramstr), paramstr = char(paramstr); end
        for i=1:size(paramstr,1) 
            cprintf('red','%s \n',paramstr(i,:)); 
        end
    end

% 
%     if ~s.analyses_completed && s.import_finished
%         fprintf('\n\n')
%         cprintf('*red','Important: To enter external application mode, all initialized analyses need to be complete & trained.\nIncomplete analyses can otherwise be deleted in the ''initialize &  manage analysis'' menu')
%     end

    if numel(NMfields) == 1 && strcmp(NMfields{1},'defs')
        mn_str = [ mn_str ...
            '|Load NeuroMiner structure' ...
            '|Change working directory' ...
            '|Investigate sample size' ...
            '|Utilities'];
        mn_act = [ mn_act 11, 13, 99, 15 ];
    else
        mn_str = [ mn_str ...
            '|Load NeuroMiner structure' ...
            '|Save NeuroMiner structure' ...
            '|Change working directory' ...
            '|Investigate sample size' ...
            '|Utilities'];
        mn_act = [ mn_act 11 12 13 99 15 ];
    end
    
    mn_sel = nk_input(menutitle,0,'mq',mn_str,mn_act);
    
    mn_opt = {'loaddata', ...
              'config', ...
              'initanal', ...
              'preproc', ...
              'mlopt', ...
              'visual', ...
              'disptrain', ...
              'inspect', ...
              'lock', ...
              'oocv', ...
              'loadmat', ...
              'savemat', ...
              'changepwd', ...
              'help', ...
              'utilities', ...
              'export', ...
              'update', ...
              'interpret'};
    
    switch mn_sel
        case 0
            action = 'quit';
            return
        case 99 %sample size simulation
            action = 'simulateML';

        otherwise
            action = mn_opt{mn_sel};
            
    end

    switch action

        case 'loaddata'
            
            if s.analyses_locked
                NM = nk_DefineOOCVData_config(NM, 1,'MAIN');
            else
                act = Inf; 
                while act>0, 
                    [NM, act] = nk_DataIO3_config(NM, 'MAIN'); 
                end
            end
            
        case 'inspect'
            nk_SelectVariateIndex(NM,1,[],0);
            nk_input('Press enter to return to the main menu',0,'sq');
            
        case 'config'
            if s.oocv_data_ready && s.analyses_locked
                act = 1; while act>0, [NM.TrainParam, act] = nk_OOCV_config(NM.TrainParam); end
            else
                act = 1; varind = []; while act>0, [act, varind] = nk_TrainClass_config([], varind,'MAIN'); end
            end
            
        case 'initanal'
             t_act = 'loop'; A = []; while ~strcmp(t_act,'BACK'), [t_act, NM, A] = nk_InitAnalysisPrep(NM, A, 'MAIN INTERFACE >> ANALYSIS MANAGER'); end
            
        case 'preproc'
            BATCH = false; p = []; act = 1; analind = []; GridAct = []; while act>0, [act, analind, p, GridAct] = nk_PreprocessPrep(act, analind, GridAct, p, 'MAIN' ); end; 
            
        case 'mlopt'    
            if isfield(NM,'analysis'), act = 1;  inp = []; while act>0, [act, inp] = nk_MLOptimizerPrep(act, inp, 'MAIN INTERFACE >> ML TRAINING MODULE'); end; end
            delete(findobj('Name','NM Optimization Status Viewer'));
            
        case 'disptrain'
            if isfield(NM,'analysis')
                nk_PrintResults2('AnalysisIndex',1, 'NM',NM); 
            end
            
        case 'lock'
            NM.defs.analyses_locked = nk_input('Do you want to finish the model discovery phase',0,'yes|no',[1,0],2);

        case 'visual'
            if isfield(NM,'analysis')
                inp = []; act = 1; while act>0, [act, inp] = nk_VisModelsPrep(act, inp, 'MAIN INTERFACE >> VISUALIZE MODELS'); end
            end

        case 'interpret'
            if isfield(NM,'analysis') 
                if isfield(NM.defs,'analyses_locked') && NM.defs.analyses_locked
                    titstr = 'MAIN INTERFACE >> INTERPRET MODELS'' PREDICTIONS [application mode]';
                else
                    titstr = 'MAIN INTERFACE >> INTERPRET MODELS'' PREDICTIONS [discovery mode]';
                end
                inp = []; act = 1; while act>0, [act, NM, inp] = nk_MLInterpreterPrep(NM, act, inp, titstr); end
            end
            
        case 'oocv'
            if isfield(NM,'analysis')
                inp = []; act = 1;
                while act>0, 
                    [act, NM, inp ] = nk_OOCVPrep(NM, act, inp, 'MAIN INTERFACE >> APPLY MODELS TO NEW DATA'); 
                end
            end
            
        case 'loadmat'
            NM = loadmat(NM);
            NMx = NM;
            assignin('base', 'NM', NMx);
            
        case 'savemat'
            savemat(NM);
        
        case 'changepwd'
            directoryname = uigetdir(pwd, 'Pick a new working directory');
            if directoryname, cd(directoryname); end
        
        case 'utilities'
            nk_Utilities
            
        case 'export' 
            scrambleflag = nk_input('Are you sure you want to shred all input features?',0,'yes|no',[1,0],2);
            if scrambleflag
                for i=1:numel(NM.Y)
                    [m,n] = size(NM.Y{i}); NM.Y{i} = rand(m,n);
                end
                NM.defs.data_scrambled = 1;
                savemat(NM);
                msgbox('Data successfully scrambled to generate anonymized NM structure!')
            end
            
        case 'update'
            complvec = []; for z=1:numel(NM.analysis), if NM.analysis{z}.status, complvec = [ complvec z ]; end; end
            t_act = 1; brief = 1; analind = 1; showmodalvec = []; 
            while t_act>0, 
                [t_act, analind, ~, showmodalvec , brief] = nk_SelectAnalysis(NM, 0, 'MAIN INTERFACE >> UPDATE ANALYSES ROOT DIRECTORIES ', analind, [], 1, showmodalvec, brief); 
            end
            if ~isempty(analind), 
                analind = complvec(analind);
            else
                analind = complvec;
            end
            newdir = nk_DirSelector('Update analyses'' root paths');
            NM = nk_UpdateRootPaths(NM, analind, newdir);
        case 'simulateML'
            if s.analyses_locked
                % use information from your sample TO DO
                %simulationRes = simulate_sampleSize(NM);
            else
                
                
                simulationRes = simSample_App2(NM);
                %simulationRes = simulate_sampleSize();
            end
    end

catch ERR
    return
end

end

function NM = loadmat(NM)

[filename, pathname] = uigetfile('*.mat', 'Pick a NeuroMiner structure file');

if ~isequal(filename,0) && ~isequal(pathname,0)
    NM = [];
    fprintf('\nLoading %s ... please wait',filename)
    load(fullfile(pathname,filename))
    if ~exist('NM','var')
        retryflag = ...
            nk_input([filename ' does not contain a valid NeuroMiner structure! Retry ?'], ...
            0, 'yes|no',[1,0],1);
        if retryflag, loadmat, end
    else
        cd(pathname)
    end
end

end

function savemat(NM)

[filename, pathname] = uiputfile('*.mat', 'Save NeuroMiner structure file');

if isfield(NM.defs,'JTextArea')
    NM = nk_InitNMwindowColors(NM, [1 1 1 ]);
    NM.defs = rmfield(NM.defs,'JTextArea');
end
if ~isequal(filename,0) && ~isequal(pathname,0)
    fprintf('\nSaving %s ... please wait',filename)
    save(fullfile(pathname,filename),'NM')
end

end

function QuitNeuroMiner()
global NM

fprintf('\n'); %clc; 
fprintf('Good Bye... \n');
delete(findobj('Tag','PrintCVBarsBin'));
%delete(findobj('Name','NM Results Manager'));
delete(findobj('Name','Probabilistic Feature Exraction'));
clear CALIBAVAIL OOCVAVAIL
if isfield(NM,'runtime'), NM = rmfield(NM,'runtime'); end
NMx = NM;
clearvars -global NM st
assignin('base', 'NM', NMx)
if exist('temp.nii','file'); delete('temp.nii'); end
if exist('tempMLI.nii','file'); delete('tempMLI.nii'); end
try warning('off','last'); catch; end
% this does not work, but I wanted to implement sth that closes the 
% Result Viewer as well as the MLI Viewer when NM is closed 

%if isfield(handles,'figure1'); handles.figure1.delete; end 
%if isfield(handles, 'MLIapp'); handles.MLIapp.delete; end
end
