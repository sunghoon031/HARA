function out = RandomRotation(angle_deg)
    axis = rand(3,1)-0.5;
    axis = axis/norm(axis);
    angle = angle_deg/180*pi;
    rotvec = angle*axis;
    out = ExpMap(rotvec);
end