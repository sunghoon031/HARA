function [val_sorted, freq_sorted] = SortNonFamilyNeighbors(x, isFamily)

    % Sort the (non-family) nodes in the order of descending frequency.
    
    if (isempty(x))
        val_sorted = [];
        freq_sorted = [];
        return;
    end
    
        
    val = [];
    freq = [];
    for i = 1:length(x)
        xi = x(i);
        if (isFamily(xi))
            continue;
        end
        
        idx = find(val==xi);
        if (isempty(idx))
            val(end+1) = xi;
            freq(end+1) = 1;
        else
            freq(idx) = freq(idx) + 1;
        end
    end
    
    if (isempty(freq))
        val_sorted = [];
        freq_sorted = [];
        return;
    end
    
    [freq_sorted, sort_idx] = sort(freq, 'descend');
    val_sorted = val(sort_idx);

end