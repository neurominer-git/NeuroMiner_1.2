function [NM, act] = nk_SavingOptions_config(NM, defaultsfl, parentstr)

saving = 0;
switch NM.modeflag
    case 'classification'
        pref = 'ClassModel';
        if numel(NM.groupnames)>2
            grname = 'Multi';
        else
            grname = sprintf('%s-%s', NM.groupnames{1}, NM.groupnames{2});
            grname = regexprep(grname,' ','_');
        end
    case 'regression'
        pref = 'RegrModel';
        grname = regexprep(NM.groupnames{1},' ', '_');
end
matname = sprintf('%s_%s', pref, grname);
act = 0;

if ~defaultsfl
    if isfield(NM.TrainParam,'SAV')
        if isfield(NM.TrainParam.SAV,'savemodel'), saving = NM.TrainParam.SAV.savemodel; end
        if isfield(NM.TrainParam.SAV,'matname'), matname = NM.TrainParam.SAV.matname; end
    end
    if ~saving, SAVING_STR = 'no'; else SAVING_STR = 'yes'; end
    menustr = [ sprintf('Save models to disk [ %s ]|', SAVING_STR) ...
                sprintf('Define prefix of NM output files [ %s ]', matname)];
    menuact = [1:2];
    
    nk_PrintLogo
    mestr = 'Define saving options'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); fprintf('\nYou are here: %s >>> ',parentstr);
    act = nk_input(mestr, 0, 'mq', menustr, menuact);
    switch act
        case 1
            if saving, saving = 0; else, saving = 1; end 
        case 2
            matname = nk_input('Define prefix of NM output files',0,'s', matname); 
    end
end
NM.TrainParam.SAV.savemodel = saving;
NM.TrainParam.SAV.matname = matname;