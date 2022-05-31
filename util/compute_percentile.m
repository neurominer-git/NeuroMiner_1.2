function centile = compute_percentile(data,value)
    
    data = data(:)';
    value = value(:);
    
    nless = sum(data < value, 2);
    nequal = sum(data == value, 2);
    centile = 100 * (nless + 0.5.*nequal) / length(data);
    
end