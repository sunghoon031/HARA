function [R_gt, edge_IDs, RR] = CreateSyntheticData(nViews, connectivity, outlier_ratio, inlier_noise_deg)


    % Input:
    % - nViews = Number of absolute rotations (views) we simulate.
    % - connectivity = Proportion of the established edges among all
    %                  pairs of views.
    % - outlier_ratio = Outlier ratio in the edges.
    % - inlier_noise_deg = Noise level in the inlier edges.
    %
    % Output:
    % - R_gt = Ground-truth absolute rotations.
    % - edge_IDs = Edges between views.
    %              If an i-th edge is established between view j and k, 
    %              edge_IDs(:,i) = [j;k].
    % - RR = Relative rotations for all edges. 
    %        If an i-th edge is established between view j and k, 
    %        RR(:,:,i) = Estimate of (R_j)*(R_k)^T.
    
    
    R_gt = cell(1, nViews);
    for i = 1:nViews
        R_gt{i} = RandomRotation(rand(1)*360);
    end

    nEdgesThr = round(connectivity*(nViews*(nViews-1)/2));

    edge_IDs = [];
    RR = [];

    interval = 1;
    c = 0;
    while (true)

        for i = 1:nViews
            j = i + interval;
            if (j > nViews)
                j = j - nViews;
            end

            c = c + 1;

            edge_IDs(:,c) = [i;j];
            RR(:,:,c) = R_gt{i}*R_gt{j}'*RandomRotation(normrnd(0, inlier_noise_deg));


            if (size(edge_IDs,2) > nEdgesThr)
                break;
            end
        end

        if (size(edge_IDs,2) > nEdgesThr)
            break;
        end

        if (interval == 1)
            nSafeEdges = c;
        end

        interval = interval + 1;

        %disp(['Progress = ' num2str(size(edge_IDs,2)/nEdgesThr*100) '%'])
    end

    nEdges = size(edge_IDs,2);


    % Turn some of the edges into outliers
    nOutliersThr = round(outlier_ratio*nEdges);
    outlierIDs = [];

    while (true)

        i = randi(nEdges);

        if (i <= nSafeEdges)
            continue;
        end

        if (ismember(i, outlierIDs))
            continue;
        end
        
        outlierIDs(end+1) = i;
        RR(:,:,i) = RandomRotation(rand(1)*360);
        
        if (length(outlierIDs) >= nOutliersThr)
            break;
        end
    end

    % Shuffle
    ii = randperm(nEdges);
    edge_IDs = edge_IDs(:,ii);
    RR = RR(:,:,ii);

end