function m = nk_SVEN(SVM, X, Y, Params)
global TRAINFUNC

[~,p] = size(X);
m.t = Params(1); m.lambda = Params(2);
Xnew = [bsxfun(@minus,X,Y./m.t) bsxfun(@plus,X,Y./m.t)]';
Ynew = [ones(p,1);-ones(p,1)];
m.C = 1/(2*m.lambda); cParams = m.C;
if ~isempty(TRAINFUNC)
	model = TRAINFUNC( Ynew, Xnew, ModelOnly, cParams);
	m.w = nk_GetPrimalW(model);
else
	CMDSTR = nk_DefineCmdStr(SVM, 'classification');
    %% Build parameter string
    cmdstr = CMDSTR.simplemodel;
    cParams = num2str(cParams','%1.10f');
    for i = 1:numel(CMDSTR.ParamStr)
        cmdstr = [ ' -' CMDSTR.ParamStr{i} ' ' strtrim(cParams(i,:)) cmdstr ];
    end
	model = train_liblin244(Ynew, sparse(Xnew), cmdstr);
    m.w = model.w; m.w = m.w';
end
m.alpha = m.C * max(1- Ynew.*(Xnew * m.w),0);
m.beta = m.t *(m.alpha(1:p) - m.alpha(p+1:2*p)) /sum(m.alpha);