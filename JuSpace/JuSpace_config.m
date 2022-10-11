function [JUSPACE, act ] = JUSPACE_config(JUSPACE, brainmask, parentstr, defaultsfl)

% from NM structure
% brainmask 
% badcoords

% user input
Atlas = [];
Option1 = 4;
CorType = 1; % Spearman
Option3 = 0;
AutoCorCorrect = 1; % correct for spatial correlations
Option5 = 0;
PetList = [];

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if isfield(JUSPACE,'atlas') && ~isempty(JUSPACE.atlas) 
        Atlas = JUSPACE.atlas; 
    end
    if isfield(JUSPACE,'option1') && ~isempty(JUSPACE.option1) 
        Option1 = JUSPACE.option1; 
    end
    % correlation type
    if isfield(JUSPACE,'cortype') && ~isempty(JUSPACE.cortype) 
        CorType = JUSPACE.cortype; 
    end
    if isfield(JUSPACE,'option3') && ~isempty(JUSPACE.option3) 
        Option3 = JUSPACE.option3; 
    end
    % adjust for spatial correlations
    if isfield(JUSPACE,'autocorcorrect') && ~isempty(JUSPACE.autocorcorrect) 
        AutoCorCorrect = JUSPACE.autocorcorrect; 
    end
    if isfield(JUSPACE,'option5') && ~isempty(JUSPACE.option5) 
        Option5 = JUSPACE.option5; 
    end
    if isfield(JUSPACE,'petList') && ~isempty(JUSPACE.petList) 
        PetList = JUSPACE.petList; 
    end
    if ~isempty(Atlas)
        AtlasDef = 1;
        ATLASSTR = Atlas;
    else
        AtlasDef = 2;
        ATLASSTR = 'not defined';
    end
 
    if ~isempty(CorType) % correlation type
        switch CorType
            case 1
                CORTYPESTR = 'Spearman correlation';
            case 2
                CORTYPESTR = 'Pearson correlation';
            case 3
                CORTYPESTR = 'Multiple linear regression';
        end
    else
        CORTYPESTR = 'not defined';
    end

    if ~isempty(AutoCorCorrect) % adjust for autocorrelations 1 = yes
        AUTOCORCORRECTSTR = num2str(AutoCorCorrect);
        if AutoCorCorrect
            AUTOCORCORRECTSTR = 'yes';
        else
            AUTOCORCORRECTSTR = 'no';
        end
    else
        AUTOCORCORRECTSTR = 'not defined';
    end
   
    if ~isempty(PetList)
        for i= 1:size(PetList,2)
            if i == 1
                PETLISTSTR = PetList{i}.id;
            else 
                PETLISTSTR = sprintf('%s, %s', PETLISTSTR, PetList{i}.id);
        
            end
        end
    else
        PETLISTSTR = 'no neurotransmitter selected';
    end

    menustr = ['Atlas                       [' ATLASSTR ']|' ...    % 'Option 1                           [' OPTION1STR ']|' ...
        'Correlation type                   [' CORTYPESTR ']|' ...  % 'Option 3                           [' OPTION3STR ']|' ...
        'Adjust for spatial correlations    [' AUTOCORCORRECTSTR ']|' ...  % 'Option 5                           [' OPTION5STR ']|' ...
        'Pet List                           [' PETLISTSTR ']'];
    menuact = [1 2 3 4];

    nk_PrintLogo
    mestr = 'JuSpace Toolbox step setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr);
    act = nk_input(mestr,0,'mq', menustr, menuact);

    switch act
        case 1
            hdrstr = 'Select atlas';
            startDir = what('JuSpace_v1.3'); 
            Atlas = nk_FileSelector(1, 'nifti', hdrstr, '.*\.nii$', [], sprintf('%s/atlas', startDir.path));
            %JUSPACE.atlas = Atlas;
        case 2
            CorType = nk_input('Define measure to quantify similarity between variables', 0, 'mq', ...
                ['Spearman correlation |', ...
                'Pearson correlation |', ...
                'Multiple linear regression'], 1:3, 0);
            %JUSPACE.cortype = CorType;
            %PX = nk_AddParam(CorType, 'CorType', 1, []);
        case 3
            AutoCorCorrect = nk_input('Define number of parameters',0,'e');
            %JUSPACE.autocorcorrect = AutoCorCorrect;
            %PX = nk_AddParam(AutoCorCorrect, 'AutoCorCorrect', 1, []);
        case 4
            PetList = print_petmaps_quickselector();
%             hdrstr = 'Select PET maps';
%             startDir = what('JuSpace_v1.3'); 
%             PetList = nk_FileSelector(NaN, 'nifti', hdrstr, '.*\.nii$', [], sprintf('%s/PETatlas', startDir.path));
            
    end
    

else
    act = 0;
end

JUSPACE.atlas = Atlas;
JUSPACE.petList = PetList;
%JUSPACE.option1 = Option1;
JUSPACE.cortype = CorType;
%JUSPACE.option3 = Option3;
JUSPACE.autocorcorrect = AutoCorCorrect;
%JUSPACE.option5 = Option5;
JUSPACE.brainmask = brainmask;


end


function neurotransmitterSel = print_petmaps_quickselector()

neurotransmitter{1}.id = '5HT1a';
neurotransmitter{1}.listidx = 1;
neurotransmitter{2}.id = '5HT1b';
neurotransmitter{2}.listidx = 5;
neurotransmitter{3}.id = '5HT21';
neurotransmitter{3}.listidx = 7;
neurotransmitter{4}.id = '5HT4';
neurotransmitter{4}.listidx = 8;
neurotransmitter{5}.id = '5HT6';
neurotransmitter{5}.listidx = 9;
neurotransmitter{6}.id = '5HTT';
neurotransmitter{6}.listidx = 10;
neurotransmitter{7}.id = 'A4B2';
neurotransmitter{7}.listidx = 13;
neurotransmitter{8}.id = 'CB1';
neurotransmitter{8}.listidx = 15;
neurotransmitter{9}.id = 'CBF';
neurotransmitter{9}.listidx = 16;
neurotransmitter{10}.id = 'CMRGlu';
neurotransmitter{10}.listidx = 19;
neurotransmitter{11}.id = 'D1';
neurotransmitter{11}.listidx = 20;
neurotransmitter{12}.id = 'D2';
neurotransmitter{12}.listidx = 23;
neurotransmitter{13}.id = 'DAT';
neurotransmitter{13}.listidx = 27;
neurotransmitter{14}.id = 'FDOPA';
neurotransmitter{14}.listidx = 29;
neurotransmitter{15}.id = 'GABAa';
neurotransmitter{15}.listidx = 31;
neurotransmitter{16}.id = 'H3';
neurotransmitter{16}.listidx = 32;
neurotransmitter{17}.id = 'M1';
neurotransmitter{17}.listidx = 33;
neurotransmitter{18}.id = 'MU';
neurotransmitter{18}.listidx = 34;
neurotransmitter{19}.id = 'NET';
neurotransmitter{19}.listidx = 37;
neurotransmitter{20}.id = 'NMDA';
neurotransmitter{20}.listidx = 38;
neurotransmitter{21}.id = 'SV2A';
neurotransmitter{21}.listidx = 39;
neurotransmitter{22}.id = 'VAChT';
neurotransmitter{22}.listidx = 40;
neurotransmitter{23}.id = 'mGluR5';
neurotransmitter{23}.listidx = 45;
neurotransmitter{24}.id = 'rCST';
neurotransmitter{24}.listidx = 46;
neurotransmitter{25}.id = 'rPWR';
neurotransmitter{25}.listidx = 47;

nk_PrintLogo
fprintf('\n\t'); fprintf('============================================= ');
fprintf('\n\t'); fprintf('***        Neurotransmitter Selector        *** ');
fprintf('\n\t'); fprintf('============================================= ');
for i=1:numel(neurotransmitter)
    fprintf('\n\t** [ %2g ]: ID : %s', i, neurotransmitter{i}.id);
end
fprintf('\n');
petind = nk_input('Type sequence of neurotransmitters to include (1d-vector)',0,'e');

% remove invalid numbers
zeroIDX = petind == 0; 
petind = petind(~zeroIDX);
greaterIDX = petind > numel(neurotransmitter); 
petind = petind(~greaterIDX);

neurotransmitterSel = [];
for i = 1:numel(petind)
    neurotransmitterSel{i}.id = neurotransmitter{petind(i)}.id;
    neurotransmitterSel{i}.listidx = neurotransmitter{petind(i)}.listidx;

end
    
end