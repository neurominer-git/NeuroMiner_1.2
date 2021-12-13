function [op, str, sortdir, optstart, minmaxstr, minmaxfl] = nk_ReturnEvalOperator(GridParam)

switch GridParam
    case {1, 2, 5, 6, 7, 10, 13, 14, 16, 17, 15, 19}
        op = 'ge'; str='above'; sortdir = 'descend'; minmaxstr = 'max'; minmaxfl = 1;
        switch GridParam
            case 19
                optstart = -100;
            otherwise
                optstart = 0;
        end
    otherwise
        op = 'le'; str='below'; sortdir = 'ascend'; optstart = Inf; minmaxstr = 'min'; minmaxfl = 2;
end
   

