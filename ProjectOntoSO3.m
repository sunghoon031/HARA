function R = ProjectOntoSO3(M)   
    [U,~,V] = svd(M);
    R = U*V.';
    if (det(R) < 0)
        V(:,3) = -V(:,3);
        R = U*V.';
    end
end