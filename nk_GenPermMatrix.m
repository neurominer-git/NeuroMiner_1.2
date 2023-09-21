function indpermA = nk_GenPermMatrix(CV, SAV, inp)

% Generate or load permutation matrix from file.
if ~isfield(inp.PERM,'suffix') || isempty(inp.PERM.suffix)
    suffix = [ '_VISpermmat_ID' inp.id ];
else
    suffix = sprintf('%s%s',inp.PERM.suffix,inp.id);
end

permfile = fullfile(inp.rootdir,[SAV.matname suffix '.mat']);

if ~exist(permfile,'file')
    for h=1:inp.nclass
        fprintf('\nCreating parent permutation matrix with %g perms', inp.PERM.nperms(1))
        if inp.nclass>1
            if h==1, indpermA = cell(1, inp.nclass); end
            if numel(CV.class{1,1}{1}.groups)==1
                lb = [find(inp.labels == CV.class{1,1}{h}.groups(1)); find(inp.labels ~= CV.class{1,1}{h}.groups(1))];  
            else
                lb = [find(inp.labels == CV.class{1,1}{h}.groups(1)); find(inp.labels == CV.class{1,1}{h}.groups(2))]; 
            end
            indpermA{h} = zeros(size(inp.labels,1), inp.PERM.nperms(1));
            indpermAh = nk_VisXPermHelper('genpermlabel', size(lb,1), inp.PERM.nperms(1));
            for ii=1:inp.PERM.nperms(1)
                indpermA{h}(lb,ii) = lb(indpermAh(:,ii));
            end
        else
            indpermA = nk_VisXPermHelper('genpermlabel', size(inp.labels,1), inp.PERM.nperms(1));
        end
    end
    save(permfile,'indpermA');
else
    fprintf('\nLoading parent permutation matrix from file: \n%s', permfile);
    load(permfile,'indpermA');
    if iscell(indpermA)
        [~, nsA] = size(indpermA{1});
        if nsA ~= inp.PERM.nperms(1)
            error('Number of permutations in matrix (%g) does not match number of expected permutations (%g)', size(indpermA,2), inp.PERM.nperms(1))
        end
    else
        [msA, nsA] = size(indpermA);
        if nsA ~= inp.PERM.nperms(1)
            error('Number of permutations in matrix (%g) does not match number of expected permutations (%g)', size(indpermA,2), inp.PERM.nperms(1))
        end
        if msA ~= size(inp.labels,1) 
            error('Number of cases in permutation matrix (%g) does not match number of cases in NM file (%g)', size(indpermA,1), size(inp.labels,1))
        end
    end
end
