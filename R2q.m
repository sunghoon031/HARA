function q = R2q(R)
    % https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    q = zeros(4,1);
    tr = trace(R);
    if (tr > 0)
        s = sqrt(1+tr)*2; %s = 4*qw
        q(1) = 0.25*s;
        q(2) = (R(3,2)-R(2,3))/s;
        q(3) = (R(1,3)-R(3,1))/s;
        q(4) = (R(2,1)-R(1,2))/s;
    elseif (R(1,1) > R(2,2) && R(1,1) > R(3,3))
        s = sqrt(1+R(1,1)-R(2,2)-R(3,3))*2; %s = 4*qx
        q(1) = (R(3,2)-R(2,3))/s;
        q(2) = 0.25*s;
        q(3) = (R(1,2)+R(2,1))/s;
        q(4) = (R(1,3)+R(3,1))/s;
    elseif (R(2,2) > R(3,3))
        s = sqrt(1+R(2,2)-R(1,1)-R(3,3))*2; %s = 4*qy
        q(1) = (R(1,3)-R(3,1))/s;
        q(2) = (R(1,2)+R(2,1))/s;
        q(3) = 0.25*s;
        q(4) = (R(2,3)+R(3,2))/s;
    else
        s = sqrt(1+R(3,3)-R(1,1)-R(2,2))*2; %s = 4*qz
        q(1) = (R(2,1)-R(1,2))/s;
        q(2) = (R(1,3)+R(3,1))/s;
        q(3) = (R(2,3)+R(3,2))/s;
        q(4) = 0.25*s;
    end
end