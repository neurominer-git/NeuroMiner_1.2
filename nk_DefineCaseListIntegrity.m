function [passed, remap] = nk_DefineCaseListIntegrity(cases, tcases)

passed = 1;
nT = numel(tcases); nC = numel(cases);
remap = [];
if ~isempty(cases)

    fprintf('\n*** Case list integrity check ***')
    if ~isequal(cases,tcases)
        for i=1:nT
            if ~isequal(cases{i}, tcases{i}), a = '*'; else a = ''; end
            fprintf('\n'); cprintf('red','old: %s, new: %s\t%s ', cases{i}, tcases{i},a)
        end
        fprintf('\n=====================================================================\n');
        warning('Existing case identifiers are not identical to current ones.')
        fprintf('\n\nYou should know EXACTLY what you do right now !!!')
        abortflag = spm_input('Actions',0,'m','Ignore and Continue|Try to match new case identifiers to old ones, fill missing cases with NaN|Abort',0:2);
        switch abortflag 
            case 0
                passed = 1;
            case 1
                remap.vec1 = nan(nC,1);
                remap.vec2 = nan(nC,1); 
                ccases = char(cases);
                ctcases = char(tcases);
                if size(ctcases,2)>size(ccases,2)
                    remap.dir = 1;
                    for i=1:nC
                        Ix = find(contains(tcases, cases{i}));
                        if ~isempty(Ix)
                            remap.vec1(i) = Ix(1);
                            remap.vec2(i) = i;
                        end
                    end
                else
                    remap.dir = 2;
                    for i=1:nT
                        Ix = find(contains(cases, tcases{i}));
                        if ~isempty(Ix)
                            remap.vec2(i) = Ix(1);
                            remap.vec1(i) = i;
                        end
                    end
                end
                remap.vec1(isnan(remap.vec1))=[];
                remap.vec2(isnan(remap.vec2))=[];
                passed = 1;
            case 2
                passed = 0; 
        end
    else
        cprintf('green','\npassed.')
    end
end