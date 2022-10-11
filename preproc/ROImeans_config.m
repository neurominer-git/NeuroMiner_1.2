function [ROIMEANS, act ] = ROImeans_config(ROIMEANS, brainmask, parentstr, defaultsfl)


Atlas = [];

if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    if isfield(ROIMEANS,'atlas') && ~isempty(ROIMEANS.atlas)
        Atlas = ROIMEANS.atlas;
    end
    if ~isempty(Atlas)
        AtlasDef = 1;
        ATLASSTR = Atlas;
    else
        AtlasDef = 2;
        ATLASSTR = 'not defined';
    end

    menustr = ['Atlas                       [' ATLASSTR ']'];
    menuact = [1];
    nk_PrintLogo
    mestr = 'ROIMEANS Toolbox step setup'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>> ',parentstr);
    act = nk_input(mestr,0,'mq', menustr, menuact);

    switch act
        case 1
            hdrstr = 'Select atlas';
            startDir = what('ROIMEANS_v1.3');
            Atlas = nk_FileSelector(1, 'nifti', hdrstr, '.*\.nii$', [], sprintf('%s/atlas', startDir.path));
    end
else
    act = 0;
end

ROIMEANS.atlas = Atlas;
ROIMEANS.brainmask = brainmask;
end