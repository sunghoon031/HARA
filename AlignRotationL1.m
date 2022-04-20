function [R_geo1, errors, mean_error, rms_error] = AlignRotationL1(R_true, R_est)

    
    nViews = size(R_true,2);
    
    errors = zeros(1,nViews);
    R_transform = cell(1, nViews);
    
    for i = 1:nViews
        R_transform{i} = R_est{i}'*R_true{i};
    end
  
    vectors_total = zeros(9,nViews);
    for i = 1:nViews
        vectors_total(:,i)= R_transform{i}(:);
    end
    med_vectors_total = median(vectors_total,2);
    [U,~,V] = svd(reshape(med_vectors_total, [3 3]));
    R_med = U*V.';
    if (det(R_med) < 0)
        V(:,3) = -V(:,3);
        R_med = U*V.';
    end

    R_geo1 = R_med;
    for j = 1:10
        step_num = 0;
        step_den = 0;
        for i = 1:nViews
            v =  LogMap(R_transform{i}*R_geo1');
            v_norm = norm(v);
            step_num = step_num + v/v_norm;
            step_den = step_den + 1/v_norm;
        end
        delta = step_num/step_den;
        delta_angle = norm(delta);

        delta_axis = delta/delta_angle;
        so3_delta = SkewSymmetricMatrix(delta_axis);
        R_delta = eye(3)+so3_delta*sin(delta_angle)+so3_delta^2*(1-cos(delta_angle));
        R_geo1 = R_delta*R_geo1;
    end
    
    

    for i = 1:nViews
        error = abs(acosd((trace(R_true{i}*(R_est{i}*R_geo1)')-1)/2));
        errors(i) = error;
    end
    mean_error = mean(errors);
    rms_error = sqrt(mean(errors.^2));
end