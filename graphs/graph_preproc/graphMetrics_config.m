function [ GRAPHMETRICS, PX, act ] = graphMetrics_config(GRAPHMETRICS, PX, parentstr, defaultsfl)

MetricsYes = 1; 

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
   

    if isfield(GRAPHMETRICS,'metric'), MetricsYes = GRAPHMETRICS.metrics; end

    if MetricsYes == 1, METRICSSTR = 'yes'; else METRICSSTR = 'no'; end
        
    menustr = ['Compute network metrics [' METRICSSTR ']'];
    menuact = 1;
    
    nk_PrintLogo
    mestr = 'Graph metrics'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            if MetricsYes == 1, MetricsYes = 2; elseif MetricsYes == 2, MetricsYes = 1; end
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
