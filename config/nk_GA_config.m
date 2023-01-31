function [act, param] = nk_GA_config(param, defaultsfl, parentstr)

GA.N                = 10;   % Number of chromosomes
GA.max_Iter         = 100;  % Number of iterations
GA.CR               = 0.5;  % Cross-over ratio
GA.MR               = 0.5;  % Mutation ratio

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end
if ~defaultsfl
    
    if isfield(param,'GA') 
        if isfield(param.GA,'N') && ~isempty( param.GA.N ), GA.N = param.GA.N; end
        if isfield(param.GA,'max_Iter') && ~isempty( param.GA.max_Iter ), GA.max_Iter = param.GA.max_Iter; end
        if isfield(param.GA,'CR') && ~isempty( param.GA.CR), GA.CR = param.GA.CR; end
        if isfield(param.GA,'MR') && ~isempty( param.GA.MR ), GA.MR = param.GA.MR; end
    end
    nk_PrintLogo
    mestr = 'Genetic algorithm setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);

    act = nk_input('Define Genetic Algorithm (GA) parameters',0,'mq', ...
                    [sprintf('Number of chromosomes [ %g ]|', GA.N) ...
                     sprintf('Number of iterations [ %g ]|', GA.max_Iter) ...
                     sprintf('Cross-over ratio [ %g ]|', GA.CR) ...
                     sprintf('Mutation ratio [ %g ]', GA.MR) ], 1:4);
    switch act
        case 1
            GA.N = nk_input('Define no. of chromosomes',0,'e',GA.N);
        case 2
            GA.max_Iter = nk_input('Define maximum no. of iterations',0,'e',GA.max_Iter);
        case 3
            GA.CR = nk_input('Define cross-over ratio',0,'e',GA.CR);
        case 4
            GA.MR = nk_input('Define mutation ratio',0,'e',GA.MR);
    end
else
    act = 0;
end
param.GA = GA; 
