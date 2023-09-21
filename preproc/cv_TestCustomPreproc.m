function Y = cv_TestCustomPreproc(Y, params)
for i = 1:size(params,2)
    Y = Y+params(i);
end
end