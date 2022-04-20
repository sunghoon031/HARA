function [R_geo2, mean_error, rms_error] = AlignRotationL2(R_true, R_est)
    
    nViews = size(R_true,2);
    errors = zeros(1,nViews);
    R_transform = cell(1, nViews);
    for i = 1:nViews
         R_transform{i} = R_est{i}'*R_true{i};
    end

    R_sum = zeros(3,3);
    for i = 1:nViews
        R_sum = R_sum + R_transform{i};
    end
    R_geo2 = ProjectOntoSO3(R_sum);
  
    for j = 1:10
        v = zeros(3,1);
        for i = 1:nViews
            v =  v + LogMap(R_transform{i}*R_geo2');
        end
        v = v/nViews;
        R_delta = ExpMap(v);
        R_geo2 = R_delta*R_geo2;
    end
    
    for i = 1:nViews
        error = abs(acosd((trace(R_true{i}*(R_est{i}*R_geo2)')-1)/2));
        errors(i) = error;
    end
    
    mean_error = mean(errors);
    rms_error = sqrt(mean(errors.^2));
end