function indperm = nk_VisXPermHelper(act, N, nperms, L)

s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

switch act
    case 'genpermlabel'
        indperm = zeros(N, nperms);
        for perms = 1:nperms
            indperm(:,perms) = randperm(N); 
        end
    case 'genpermfeats'
        uL = unique(L); nuL = numel(uL);
        indperm = zeros(nuL, N, nperms);
        for i = 1:nuL
            for perms = 1:nperms
                indperm(i,:,perms) = randperm(N);
            end
        end  
end


