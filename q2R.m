function R = q2R(q)
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    R = zeros(3,3);
    R(1,1) = 1 - 2*qy^2 - 2*qz^2;
    R(1,2) = 2*qx*qy - 2*qz*qw;
    R(1,3) = 2*qx*qz + 2*qy*qw;
    R(2,1) = 2*qx*qy + 2*qz*qw;
    R(2,2) = 1 - 2*qx^2 - 2*qz^2;
    R(2,3) = 2*qy*qz - 2*qx*qw;
    R(3,1) = 2*qx*qz - 2*qy*qw;
    R(3,2) = 2*qy*qz + 2*qx*qw;
    R(3,3) = 1 - 2*qx^2 - 2*qy^2;
end