function [ GRAPHMETRICS, PX, act ] = graphMetrics_config(GRAPHMETRICS, PX, parentstr, defaultsfl)

MetricsYes = 1; 
MetricsList = [];


if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
   
    if isfield(GRAPHMETRICS,'metrics'), MetricsYes = GRAPHMETRICS.metrics; end
    if isfield(GRAPHMETRICS,'metricslist') && ~isempty(GRAPHMETRICS.metricslist), MetricsList = GRAPHMETRICS.metricslist; end

    if MetricsYes == 1, 
        METRICSSTR = 'yes'; 
    else 
        METRICSSTR = 'no'; 
    end
    
    if ~isempty(MetricsList)
        for i= 1:size(MetricsList,2)
            if i == 1
                METRICSLISTSTR = MetricsList{i}.id;
            else 
                METRICSLISTSTR = sprintf('%s, %s', METRICSLISTSTR, MetricsList{i}.id);
            end
        end
    else
        METRICSLISTSTR = 'no metrics selected';
    end

        
    menustr = ['Compute network metrics [' METRICSSTR ']|' ...
        'Choose network metrics [' METRICSLISTSTR ']'];
    menuact = 1:2;
    
    nk_PrintLogo
    mestr = 'Graph metrics'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            if MetricsYes == 1, MetricsYes = 2; elseif MetricsYes == 2, MetricsYes = 1; end
        case 2
            MetricsList = print_networkmetrics_quickselector();
    end
else
    act = 0;
end
GRAPHMETRICS.metric = MetricsYes;
if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end


    function netwmetricsSel = print_networkmetrics_quickselector()

netwmetrics{1}.id = '5HT1a';
netwmetrics{2}.id = '5HT1b';
netwmetrics{3}.id = '5HT21';
netwmetrics{4}.id = '5HT4';
netwmetrics{5}.id = '5HT6';
netwmetrics{6}.id = '5HTT';
netwmetrics{7}.id = 'A4B2';
netwmetrics{8}.id = 'CB1';


nk_PrintLogo
fprintf('\n\t'); fprintf('============================================= ');
fprintf('\n\t'); fprintf('***        Network metrics Selector        *** ');
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