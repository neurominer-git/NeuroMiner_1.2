function [act, param] = nk_PSO_config(param, defaultsfl)

PSO.N                = 10;   % Number of particles
PSO.max_Iter         = 100;  % Number of iterations
PSO.c1               = 1;  % Cognitive factor
PSO.c2               = 1;  % Social factor
PSO.w                = 0.5; % Inertia weight
PSO.lb               = 0; % Lower bound for random feature initialization
PSO.ub               = 1; % upper bound for random feature initialization
PSO.thres            = 0.5; % threshold for feature selection

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end
if ~defaultsfl
    
    if isfield(param,'PSO') 
        if isfield(param.PSO,'N') && ~isempty( param.PSO.N ), PSO.c = param.PSO.N; end
        if isfield(param.PSO,'max_Iter') && ~isempty( param.PSO.max_Iter ), PSO.max_Iter = param.PSO.max_Iter; end
        if isfield(param.PSO,'c1') && ~isempty( param.PSO.c1), PSO.c1 = param.PSO.c1; end
        if isfield(param.PSO,'c2') && ~isempty( param.PSO.c2 ), PSO.c2 = param.PSO.c2; end
        if isfield(param.PSO,'w') && ~isempty( param.PSO.w ), PSO.w = param.PSO.w; end
        if isfield(param.PSO,'lb') && ~isempty( param.PSO.lb), PSO.lb = param.PSO.lb; end
        if isfield(param.PSO,'ub') && ~isempty( param.PSO.ub ), PSO.ub = param.PSO.ub; end
        if isfield(param.PSO,'thres') && ~isempty( param.PSO.thres ), PSO.thres = param.PSO.thres; end
    end
    nk_PrintLogo
    act = nk_input('Define Particle Swarm Optimization (PSO) parameters',0,'mq', ...
                    [sprintf('Number of particles [ %g ]|', PSO.N) ...
                     sprintf('Number of iterations [ %g ]|', PSO.max_Iter) ...
                     sprintf('Cognitive ratio [ %g ]|', PSO.c1) ...
                     sprintf('Social ratio [ %g ]|', PSO.c2) ... 
                     sprintf('Inertia weight [ %g ]|', PSO.w) ...
                     sprintf('Lower bound [ %g ]|', PSO.lb) ...
                     sprintf('Upper bound [ %g ]|', PSO.ub) ... 
                     sprintf('Threshold [ %g ]', PSO.thres)], 1:8);
    switch act
        case 1
            PSO.N = nk_input('Define no. of particles',0,'e',PSO.N);
        case 2
            PSO.max_Iter = nk_input('Define maximum no. of iterations',0,'e',PSO.max_Iter);
        case 3
            PSO.c1 = nk_input('Define cognitive factor',0,'e',PSO.c1);
        case 4
            PSO.c2 = nk_input('Define social factor',0,'e',PSO.c2);
        case 5
            PSO.w = nk_input('Define inertia weight',0,'e',PSO.w);
        case 6
            PSO.lb = nk_input('Define lower bound of feature weight',0,'e',PSO.lb);
        case 7
            PSO.ub = nk_input('Define upper bound of feature weight',0,'e',PSO.ub);
        case 8
            PSO.thres = nk_input('Define threshold for feature selection',0,'e',PSO.thres);
    end
else
    act = 0;
end
param.PSO = PSO; 
