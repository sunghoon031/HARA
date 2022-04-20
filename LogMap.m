function out = LogMap(in)
    if (in(1,1) == 1 && in(2,2) == 1 && in(3,3) == 1)
        out = [0;0;0];
        return;
    end
    
    cos_theta = (trace(in)-1)/2;
    sin_theta = sqrt(1-cos_theta^2);
    theta = acos(cos_theta);
    ln_R = theta/(2*sin_theta)*(in-in');
    out = [ln_R(3,2);ln_R(1,3);ln_R(2,1)];
end