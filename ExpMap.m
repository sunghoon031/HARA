function out = ExpMap(in)
    angle = norm(in);
    if (angle == 0)
        out = eye(3);
        return;
    end
    axis = in/angle;
    so3 = SkewSymmetricMatrix(axis);
    R = eye(3)+so3*sin(angle)+so3^2*(1-cos(angle));
    out = R;
end