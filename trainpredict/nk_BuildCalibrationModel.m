function model = nk_BuildCalibrationModel(SVM, MODEFL, model, Y, label)

if strcmp(MODEFL,'classification') && isfield(SVM,'BBQ') && SVM.BBQ.flag
    cmd = ['nk_GetTestPerf_' SVM.prog]; 
    [~, PTR] = feval(cmd, [], Y, label, model);
    PTR = exp(PTR)./(1+exp(PTR)); 
    if numel(unique(PTR))==1, return; end
    L=label;L(label==-1 | label==2)=0;
    m.md = model;
    m.BBQ = buildBBQ(PTR,L,[]);
    model = m;
end