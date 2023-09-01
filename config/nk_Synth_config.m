function [SYNTH, act, mess] = nk_Synth_config(SYNTH, mess, parentstr)
global NM

%% Set-up configuration interface

% Define variables
if ~exist('mess','var'), mess = []; end
na_str                      = '?';
synth_flag                  = 2;
k                           = 5;
distanceMeasure             = 'Euclidean';
standardize_data            = 2;
numSyntheticObservations    = round(size(NM.label,1)*0.1);
write2disk                  = 1;

% Generate synthetic data?
if isfield(SYNTH,'flag'),                       synth_flag = SYNTH.flag; else, SYNTH.flag = synth_flag; end
if isfield(SYNTH,'k'),                          k = SYNTH.k; else, SYNTH.k = k; end
if isfield(SYNTH,'distanceMeasure'),            distanceMeasure = SYNTH.distanceMeasure; else, SYNTH.distanceMeasure = distanceMeasure; end
if isfield(SYNTH,'standardize_data'),            standardize_data = SYNTH.standardize_data; else, SYNTH.standardize_data = standardize_data; end
if isfield(SYNTH,'numSyntheticObservations'),   numSyntheticObservations = SYNTH.numSyntheticObservations; else, SYNTH.numSyntheticObservations = numSyntheticObservations; end
if isfield(SYNTH,'write2disk'),                 write2disk = SYNTH.write2disk; else,  SYNTH.write2disk = write2disk;  end

mn_str = []; mn_act = [];

synthstr        = {'enabled','disabled'};
write2diskstr   = {'use data permanently by writing it to disk','recreate data in each model training run'};
stdstr          = {'standardize data', 'data already standardized/standardization not needed'};
mn_str = [mn_str sprintf('Generate synthetic data [ %s ]',synthstr{synth_flag})];       mn_act = [mn_act 1];
if synth_flag == 1
    mn_str = [mn_str sprintf('|Define no. of nearest neighbors for interpolation [ k=%g ]', k)];          mn_act = [mn_act 2];
    mn_str = [mn_str sprintf('|Choose distance measure [ %s ]', distanceMeasure)];      mn_act = [mn_act 3];
    mn_str = [mn_str sprintf('|Define whether data should be standardized before running distance analysis [ %s ]', stdstr{standardize_data})];      mn_act = [mn_act 4];
    mn_str = [mn_str sprintf('|Define how many synthetic observations you would like to be generated [ %g ]', numSyntheticObservations)];      mn_act = [mn_act 5];
    mn_str = [mn_str sprintf('|Choose how to manage synthetic data [ %s ]', write2diskstr{write2disk})];      mn_act = [mn_act 6];
end

nk_PrintLogo

if ~isempty(mess)
    for i=1:numel(mess)
        if isempty(mess(i).text), continue; end
        fprintf('\n');mess(i).text = regexprep(mess(i).text,'\','/');
        fprintf(mess(i).text); 
    end
    fprintf('\n')
    mess = [];
end

fprintf('\n'); mestr = 'Configure synthetc data generation options';  
navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); fprintf('You are here: %s >>> ',parentstr); 
act = char(nk_input(mestr,0,'mq', mn_str, mn_act));

switch act
    case 1
       if SYNTH.flag == 1, SYNTH.flag = 2; else, SYNTH.flag = 1; end

    case 2
        SYNTH.k = nk_input('Define no. of nearest neighbors for interpolation',0,'i', k);

    case 3
       sel = 'Euclidean|Manhattan|Cosine|Mahalanobis';
       defsel = strsplit(sel,'|');
       SYNTH.distanceMeasure = char(nk_input('Choose distance measure',0,'m', sel, defsel, find(contains(defsel, distanceMeasure))));

    case 4
        if SYNTH.standardize_data == 1, SYNTH.standardize_data = 2; else, SYNTH.standardize_data = 1; end

    case 5
        SYNTH.numSyntheticObservations = nk_input('How many synthetic observations should be generated',0,'i', numSyntheticObservations);

    case 6
        if SYNTH.write2disk == 1, SYNTH.write2disk = 2; else, SYNTH.write2disk = 1; end

end
