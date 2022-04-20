function out = UniformSampling(data, r)

    % Sample r numbers from data uniformly and deterministically.
    
    n = length(data);
    if (r >= n)
        out = data;
        return;
    end
    quotient = floor(n/r);
    remainder = n-r*quotient;
    if (remainder==0)
        out = quotient:quotient:n;
        return;
    end
    out = nan(1, r);
    i = 0;
    c = 0;
    out_idx = 0;
    while (out_idx < r)
        i = i + 1;
        if (i > length(data))
            i = 1;
            quotient = floor(length(data)/(r-out_idx));
        end
        if (c == 0)
            out_idx = out_idx + 1;
            out(out_idx) = data(i);
            data(i) = [];
            i = i-1;
        end
        c = c + 1;
        if (c == quotient+1)
            c = 0;
        end
    end
end