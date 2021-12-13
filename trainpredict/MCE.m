function [ res ] = MCE( labels, predictions )   
    ordered = sortrows([predictions,labels]);
    N = size(ordered,1);
    rest = mod(N,10);
    for i = 1 : 10
        if (i <= rest)
            group = ordered((i-1) * ceil(N / 10) + 1 : i * ceil(N / 10),:);
        else
            group = ordered(rest + (i-1) * fix(N / 10) + 1 : rest + i * fix(N / 10),:);
        end
        n = size(group,1);
        observed = mean(group(:,2));
        expected = mean(group(:,1));
        S(i) = abs(expected-observed);
    end
    res = max(S);
end

