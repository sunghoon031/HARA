function R = ChordalL1Mean(R_input, b_outlier_rejection, n_iterations, thr_convergence)
    % https://github.com/sunghoon031/RobustSingleRotationAveraging
    
    % 1. Initialize
    n_samples = length(R_input);
    
    vectors_total = zeros(9,n_samples);
    for i = 1:n_samples
        vectors_total(:,i)= R_input{i}(:);
    end
    
  
    s = median(vectors_total,2);
    
    
    % 2. Optimize
    for j = 1:n_iterations
        if (sum(sum(abs(vectors_total-s))==0) ~= 0)
            s = s+rand(size(s,1),1)*0.001;
        end

        v_norms = zeros(1,n_samples);
        for i = 1:n_samples
            v =  vectors_total(:,i)-s;
            v_norm = norm(v);
            v_norms(i) = v_norm;
        end

        % Compute the inlier threshold (if we reject outliers).
        thr = inf;
        if (b_outlier_rejection)
            sorted_v_norms = sort(v_norms);
            v_norm_firstQ = sorted_v_norms(ceil(n_samples/4));

            if (n_samples <= 50)
                thr = max(v_norm_firstQ, 1.356);
                % 2*sqrt(2)*sin(1/2) is approximately 1.356
            else
                thr = max(v_norm_firstQ, 0.7);
                % 2*sqrt(2)*sin(0.5/2) is approximately 0.7
            end
        end

        step_num = 0;
        step_den = 0;

        for i = 1:n_samples
            v_norm = v_norms(i);
            if (v_norm > thr)
                continue;
            end
            step_num = step_num + vectors_total(:,i)/v_norm;
            step_den = step_den + 1/v_norm;
        end


        s_prev = s;
        s = step_num/step_den;

        update_medvec = s-s_prev;
        if (norm(update_medvec) < thr_convergence)
            break;
        end

    end
    
    R = ProjectOntoSO3(reshape(s, [3 3]));

end