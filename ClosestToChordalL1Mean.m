function R = ClosestToChordalL1Mean(R_input, b_outlier_rejection, n_iterations, thr_convergence)
    
    R_ = ChordalL1Mean(R_input, b_outlier_rejection, n_iterations, thr_convergence);
    
    min_err = inf;
    for i = 1:length(R_input)
        delta_R = R_-R_input{i};
        err_sq = sum(sum(delta_R.^2));
        
        if (err_sq < min_err)
            min_err = err_sq;
            R = R_input{i};
        end
    end

end