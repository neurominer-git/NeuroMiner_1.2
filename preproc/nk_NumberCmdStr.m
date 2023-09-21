function [PREPROC, cmdseqold, cmdseqnew] = nk_NumberCmdStr(PREPROC)

nA = numel(PREPROC.ACTPARAM);
cmdseqold = cell(nA,1);
cmdseqnew = cell(nA,1);

for i=1:nA
    cmdseqold{i} = PREPROC.ACTPARAM{i}.cmd;
end

uCmd = unique(cmdseqold);
nuCmd = numel(uCmd);

for i=1:nuCmd
    idx = find(strcmp(cmdseqold,uCmd{i}));
    for j=1:numel(idx)
        cmdseqnew{idx(j)} = sprintf('%s_%g',cmdseqold{idx(j)},j);
    end
end

for i=1:nA
    PREPROC.ACTPARAM{i}.cmd = cmdseqnew{i};
end
