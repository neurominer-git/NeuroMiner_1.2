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
GRAPHMETRICS.metricslist = MetricsList;
if exist('PX','var') && ~isempty(PX) && ...
        isfield (PX,'Px') && ...
        isfield(PX.Px,'Params') && ...
        ~isempty(PX.Px.Params)
    PX.opt = allcomb(PX.Px.Params,'matlab'); 
else
    PX.opt = [];
end


function netwmetricsSel = print_networkmetrics_quickselector()

netwmetrics{1}.id           = 'degree (local)';
netwmetrics{1}.listidx      =  1; 
netwmetrics{2}.id           = 'strength (local)';
netwmetrics{2}.listidx      =  2; 
netwmetrics{3}.id           = 'betweenness (local)';
netwmetrics{3}.listidx      =  3; 
netwmetrics{4}.id           = 'clustering coefficient (local)';
netwmetrics{4}.listidx      =  4; 
netwmetrics{5}.id           = 'clustering coefficient (global)';
netwmetrics{5}.listidx      =  5; 
netwmetrics{6}.id           = 'diameter';
netwmetrics{6}.listidx      =  6; 
netwmetrics{7}.id           = 'transitivity (global)';
netwmetrics{7}.listidx      =  7; 
netwmetrics{8}.id           = 'eigenvector';
netwmetrics{8}.listidx      =  8; 
netwmetrics{9}.id           = 'efficiency (local)';
netwmetrics{9}.listidx      =  9; 
netwmetrics{10}.id          = 'efficiency (global)';
netwmetrics{10}.listidx     =  10; 
netwmetrics{11}.id          = 'closeness';
netwmetrics{11}.listidx     =  11; 
netwmetrics{12}.id          = 'radius';
netwmetrics{12}.listidx     =  12; 
netwmetrics{13}.id          = 'pagerank';
netwmetrics{13}.listidx     =  13; 
netwmetrics{14}.id          = 'distance';
netwmetrics{14}.listidx     =  14; 

nk_PrintLogo
fprintf('\n\t'); fprintf('============================================= ');
fprintf('\n\t'); fprintf('***        Network metrics Selector        *** ');
fprintf('\n\t'); fprintf('============================================= ');
for i=1:numel(netwmetrics)
    fprintf('\n\t** [ %2g ]: ID : %s', i, netwmetrics{i}.id);
end
fprintf('\n');
metricsind = nk_input('Type sequence of neurotransmitters to include (1d-vector)',0,'e');

% remove invalid numbers
zeroIDX = metricsind == 0; 
metricsind = metricsind(~zeroIDX);
greaterIDX = metricsind > numel(netwmetrics); 
metricsind = metricsind(~greaterIDX);

netwmetricsSel = [];
for i = 1:numel(metricsind)
    netwmetricsSel{i}.id = netwmetrics{metricsind(i)}.id;
    netwmetricsSel{i}.listidx = netwmetrics{metricsind(i)}.listidx;
    
end
    
