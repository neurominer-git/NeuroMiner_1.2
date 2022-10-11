function [act, param] = nk_Wrapper_config(act, param, SVM, MODEFL, GRD, MULTI, defaultsfl, parentstr)
%
% Defaults
% --------
% Activate feature filtering
wrappertypes                = {'Greedy feature selection', ...
                               'Simulated Annealing', ...
                               'Genetic algorithm', ...
                               'Particle swarm optimization', ...
                               'Path finder algorithm'};
wrapperidx                  = 1:numel(wrappertypes);
wrappertag                  = {'GreedySearch', 'SA', 'GA', 'PSO', 'PFA'};
wrapperconffun              = {@nk_GreedySearch_config, @nk_SA_config, @nk_GA_config, @nk_PSO_config, @nk_PFA_config};

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end

if ~isempty(act) || ~defaultsfl
    
    if ~isfield(param,'Wrapper')
        if ~exist("MULTI", "var"), MULTI=[]; end
        if ~exist("parentstr", "var"), parentstr=[]; end
        [~, param] = nk_Wrapper_config(act, param, SVM, MODEFL, GRD, MULTI, true, parentstr);
    end
    if ~isfield(param.Wrapper,'GA')
        [ ~, param.Wrapper ] = nk_GA_config(param.Wrapper,1); % Define GA defaults
    end
    if ~isfield(param.Wrapper,'SA')
        [ ~, param.Wrapper ] = nk_SA_config(param.Wrapper,1); % Define SA defaults
    end
    if ~isfield(param.Wrapper,'PSO')
        [ ~, param.Wrapper ] = nk_PSO_config(param.Wrapper,1); % Define PSO defaults
    end
    if ~isfield(param.Wrapper,'PFA')
        [ ~, param.Wrapper ] = nk_PFA_config(param.Wrapper,1); % Define PFA defaults
    end
    Wrapper = param.Wrapper; 
    if ~isfield(param.Wrapper,'EnsembleStrategy')
        param.Wrapper = nk_EnsembleStrategy2_config(Wrapper, SVM, MODEFL, 1);
    end
    %if isfield(Wrapper,'optflag'), optflag = Wrapper.optflag; else, optflag = 1; end 
            
    nk_PrintLogo
    mestr = 'Wrapper-based model selection setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
    
    if ~Wrapper.flag
        
        act = nk_input(mestr,0,'mq','Activate wrapper methods', 1, 1);
        
    else
       
        d = nk_GetParamDescription2([],param,'FeatWrap');
        
        if isfield(param,'Filter') && param.Filter.flag && isfield(param.Filter,'SubSpaceFlag') && param.Filter.SubSpaceFlag
            e = nk_GetParamDescription2([],param.Wrapper,'EnsType');
            fprintf('\n================================================================')
            fprintf('\n* Pre-Wrapper filtering/subspace optimization method detected. *')
            fprintf('\n* Wrapper will be applied to pre-defined feature subspaces.    *')
            fprintf('\n* Further subspace optimization methods enabled.               *')
            fprintf('\n================================================================')

            SubStr = ['Evaluate features within subspaces [ ' e.EnsStrat ' ]|'];
            actind = [1 2];    
        else
            SubStr ='';
            actind = 1;
        end
        
        if prod(GRD.n_params)>1 || (isfield(param,'PreML') && ~isempty(param.PreML) && (numel(param.PreML.Params)>1 || numel(param.PreML.Params{1})>1))
            GrdOptStr = ['Learn wrapper model only at parameter optimum [ ' d.WrapperOptFlag ' ]|']; actind = [actind 3 ];
        else
            GrdOptStr = '';
        end
        
        act = nk_input(mestr,0,'mq',['Deactivate wrapper methods|'...
                                   SubStr ...
                                   GrdOptStr ...
                                   'Select wrapper algorithm [ ' wrappertypes{Wrapper.type} ' ]|'...
                                   'Configure ' wrappertypes{Wrapper.type} '[ ' print_struct_params(Wrapper.(wrappertag{Wrapper.type}), '; ') ' ]|' ...
                                   'CV1 data partitions for optimization [ ' d.WrapperDataMode ' ]|'...
                                   'Cross-CV1 Feature selection [ ' d.WrapperPFE ' ]'], [actind 4:7], 1);
            
    end
    
    switch act
        
        case 1
            Wrapper.flag = ~Wrapper.flag;    
            if ~Wrapper.flag && Wrapper.optflag == 1, Wrapper.optflag = 2; end
 
        case 2
            if isfield(param,'Filter') && isfield(param.Filter,'SubSpaceFlag')
                Wrapper.SubSpaceFlag = param.Filter.SubSpaceFlag;
            else
                Wrapper.SubSpaceFlag = 0;
            end
            Wrapper = nk_EnsembleStrategy2_config(Wrapper, SVM, MODEFL, [], navistr);
            if Wrapper.SubSpaceFlag
                SubSpaceSteppingFlag = nk_input('Define feature subspace stepping',0,'yes|no',[1,0],0);
                if SubSpaceSteppingFlag
                    Wrapper.SubSpaceStepping = nk_input('Divide maximum number of feature subspaces by',0, ...
                        'w1',Wrapper.SubSpaceStepping);
                end
            end
        case 3
            if Wrapper.optflag ==1,  Wrapper.optflag = 2; else, Wrapper.optflag = 1; end
            
        case 4
            Wrapper.type = nk_input('Choose wrapper algorithm',0,'m', strjoin(wrappertypes,'|'), wrapperidx, Wrapper.type);
        case 5
            switch Wrapper.type
                case 1
                    t_act = 1; while t_act>0, [ t_act , Wrapper ] = nk_GreedySearch_config(Wrapper, SVM, MULTI, 0, navistr); end
                otherwise
                    t_act = 1; while t_act>0, [ t_act , Wrapper ] = wrapperconffun{Wrapper.type}(Wrapper, 0, navistr); end
            end
        case 6
            Wrapper.datamode = nk_input('Samples for optimization',0,'m','CV1 training data|CV1 test data|CV1 training & test data',[1,2,3], Wrapper.datamode);
        case 7
            sact = 1; while sact>0, [Wrapper.PFE, sact] = nk_ProbalisticFea_config(Wrapper.PFE, navistr); end
    end
else
    Wrapper.flag                = 0;    % Perform wrapper-based feature selection
    Wrapper.optflag             = 2;    % Perform wrapper-based feature selection at parameter optimum
    Wrapper.type                = 1;
    Wrapper.datamode            = 2;
    Wrapper.CostFun             = 2;   % currently static

    % Subspace strategy
    Wrapper.SubSpaceFlag        = 0;
    Wrapper.SubSpaceStepping    = 0;
    [ ~, Wrapper ]              = nk_GreedySearch_config(Wrapper,SVM, MULTI, 1);
    [ ~, Wrapper ]              = nk_SA_config(Wrapper,1); % Define SA defaults
    [ ~, Wrapper ]              = nk_GA_config(Wrapper,1); % Define GA defaults
    [ ~, Wrapper ]              = nk_PSO_config(Wrapper,1); % Define PSO defaults
    [ ~, Wrapper ]              = nk_PFA_config(Wrapper,1); % Define PFA defaults
    Wrapper                     = nk_EnsembleStrategy2_config(Wrapper, SVM, MODEFL, 1);
    % Probabilistic Feature Extraction (PFE) across CV1 partitions 
    Wrapper.PFE                 = nk_ProbalisticFea_config([],1);
end
param.Wrapper = Wrapper;
