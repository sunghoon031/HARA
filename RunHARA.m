function [R_est, time_initialization, time_optimization] = ...
    RunHARA(...
    nViews, ...
    edge_IDs, ...
    RR, ...
    nSamples, ...
    s_init, ...
    thrs_proportion, ...
    triplet_thr, ...
    median_thr, ...
    err_sq_thr)
    

    %% Definitions:
    % Input:
    % - nViews = Number of absolute rotations (views) we want to estimate.
    % - edge_IDs = Edges between views.
    %              If an i-th edge is established between view j and k, 
    %              edge_IDs(:,i) = [j;k].
    % - RR = Relative rotations for all edges. 
    %        If an i-th edge is established between view j and k, 
    %        RR(:,:,i) = Estimate of (R_j)*(R_k)^T.
    % - nSamples = Number of triplets we sample per edge (default = 10).
    % - s_init = Initial support threshold (default = 10).
    % - thrs_proportion = Related to the percentiles of the sampled loop
    %                     errors. We use it to set the loop thresholds.
    %                     If thrs_proportion = [0.1 0.2 0.3], it means that
    %                     we set the loop thresholds (thrs) to the 10th,
    %                     20th, and 30th percentile of the collected
    %                     errors.
    %                     (default = [0.1, 0.2, 0.3])
    % - triplet_thr = When compute the loop thresholds from the sampled
    %                 loop errors, we ignore the errors above this 
    %                 threshiold. This is because we assume that any loop 
    %                 error above this level is likely to be caused by one
    %                 or more outliers in the triplet.
    %                 (default = 1.0)
    % - median_thr = If the median of all sampled loop errors is greater
    %                than this threshold, we skip edge filtering.
    %                (default = 1.0)
    % - err_sq_thr = Corresponds to tau^2 (tau is explained in the paper).
    %                If the edge does not conform to our initial solution,
    %                we filter this edge. This is the threshold we use for
    %                this check. 
    %                (default = 1.0).
    %
    % Output:
    % - R_est = Estimated absolute rotations.
    % - time_initialization = Time taken to obtain the initial solution.
    % - time_optimization = Time taken to filter edges and perform the
    %                       local refinement.
    
    %% Step 0: Preprocessing
    
    t_init = tic;
    
    nEdges = size(edge_IDs,2);
    
    edge_info = cell(1, nViews); % Lists the views connected to each view.
    for i = 1:nEdges
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        edge_info{j} = [edge_info{j}, k];
        edge_info{k} = [edge_info{k}, j];  
    end

    nEdgesForEachView = zeros(1, nViews);
    for i = 1:nViews
        nEdgesForEachView(i) = length(edge_info{i});
    end
        
    A = zeros(nViews,nViews);
    G = eye(3*nViews);
    for i = 1:nEdges
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        R_jk = RR(:,:,i);
        G(3*j-2:3*j, 3*k-2:3*k) = R_jk;
        G(3*k-2:3*k, 3*j-2:3*j) = R_jk';
        A(j,k) = 1;
        A(k,j) = 1;
    end

    % Sample loop errors:
    errs = [];
    errs_total = [];
    for i = 1:nEdges
        j = edge_IDs(1,i);
        k = edge_IDs(2,i);
        R_jk = G(3*j-2:3*j, 3*k-2:3*k);

        ls = intersect(edge_info{j}, edge_info{k});
        ls = UniformSampling(ls, nSamples);

        for l = ls
            R_jl = G(3*j-2:3*j, 3*l-2:3*l);
            R_kl = G(3*k-2:3*k, 3*l-2:3*l);

            delta_R = R_jk - (R_jl*R_kl');
            err_sq = sum(sum(delta_R.^2));

            errs_total(end+1) = err_sq;

            if (err_sq < triplet_thr)
                errs(end+1) = err_sq;
            end
        end
    end

    med_err = median(errs_total);

    errs = sort(errs);
    nErrs = length(errs);
    thrs = errs(round(nErrs*thrs_proportion));

    
    %% Step 1: Initialization

    thr_idx = 1;
    s = s_init;

    isFamily = zeros(1, nViews);
    R_init = cell(1,nViews);

    % Find the most connected node:
    [~, id_to_check] = max(nEdgesForEachView);

    % Add it to family.
    newFamily = id_to_check;
    newFamily_nEdges = nEdgesForEachView(id_to_check);

    isFamily(id_to_check) = 1; 
    R_init{id_to_check} = eye(3);


    supportedNeighborsTable = zeros(nViews, length(thrs), s_init);
    % Example: For each view, it's the transpose of the following table:
    %
    %                                        | thr1 | thr2 | thr3 | ...
    % #nbors supported by 1 other neighbor   |  12  |  14  |  10  | ...
    % #nbors supported by 2 other neighbors  |  20  |   31 |  24  | ...
    %                 ...
    % #nbors supported by s_init neighbors   |  50  |  100 |  80  | ...
    
    votes = zeros(1,nViews);
    has_voted = zeros(1,nViews);

    while(min(isFamily)==0)
        
        % If the neighbors of the first one in newFamily are consistent
        % with each other, add them to the family.
        
        while(~isempty(newFamily))
            id_to_check = newFamily(1);
            newFamily(1) = [];
            newFamily_nEdges(1) = [];

            % Set the non-family neighbors' rotations:
            
            neighbors = edge_info{id_to_check};
            
            all_neighbors_are_family = true;
            for j = neighbors
                if (~isFamily(j))
                    all_neighbors_are_family = false;
                    R_ji = G(3*j-2:3*j, 3*id_to_check-2:3*id_to_check);
                    R_init{j} = R_ji*R_init{id_to_check};
                end
            end

            if (all_neighbors_are_family)
                supportedNeighborsTable(id_to_check, :, :) = zeros(length(thrs), s_init);
                continue;
            end

            consistent_neighbors = cell(1, length(thrs));
            % Example:
            % If (1, 2), (1, 3), (2, 3) have error less than thrs(1), 
            % and (1, 2), (1, 3), (2, 3), (2, 4), (3, 4) have error less
            % than thrs(2), then we have
            % consistent_neighbors{1} = [1, 2, 1, 3, 2, 3]
            % consistent_neighbors{2} = [1, 2, 1, 3, 2, 3, 2, 4, 3, 4]


            % Check consistency between them:
            for jj = 1:length(neighbors)-1
                j = neighbors(jj);
                for kk = jj + 1:length(neighbors)
                    k = neighbors(kk);
                    if (A(j,k)==0)
                        continue;
                    end
                    if (isFamily(j) && isFamily(k))
                        continue;
                    end

                    R_jk = G(3*j-2:3*j, 3*k-2:3*k);
                    delta_R = R_init{j}*R_init{k}'-R_jk;
                    err_sq = sum(sum(delta_R.^2));

                    for ii = 1:length(thrs)
                        if (err_sq < thrs(ii))
                            consistent_neighbors{ii} = [consistent_neighbors{ii}, j, k];
                        end
                    end
                end
            end
            
            % Next, we find out how many times each neighbor appears in
            % consistent_neighbors at the current loop threshold.
            % If a neighbor appears s or more times, we add it to family.
            
            [consistent_neighbors_sorted, freqs_sorted] = SortNonFamilyNeighbors(consistent_neighbors{thr_idx}, isFamily);
            

            family_candidates = [];
            for i = 1:length(consistent_neighbors_sorted)
                if (freqs_sorted(i) >= s)
                    if (~isFamily(consistent_neighbors_sorted(i)))
                        family_candidates(end+1) = consistent_neighbors_sorted(i);
                    end
                else
                    break;
                end
            end
                    

            for i = 1:length(family_candidates)
                j = family_candidates(i);

                isFamily(j) = 1;
                [newFamily, newFamily_nEdges] = AddToNewFamily(...
                    j, ...
                    nEdgesForEachView(j), ...
                    newFamily, ...
                    newFamily_nEdges);

                % Reset thr and s (since new family member is added).
                s = s_init;
                thr_idx = 1;
            end   
            
            
            % Now that consistent_neighbors have been updted, update 
            % supportedNeighborsTable too:
            
            for ii = 1:length(thrs)
                [~, freq_sorted] = SortNonFamilyNeighbors(consistent_neighbors{ii}, isFamily);

                for i = 1:s_init
                    supportedNeighborsTable(id_to_check, ii, i) = sum(freq_sorted >= i);
                end
            end
    
        end
        % newFamily is now empty.

        
        % Either (1) move on to the next base node, (2) increase the loop
        % threshold, or (3) decrease the support threshold, s:
        
        [max_support, id_min] = max(supportedNeighborsTable(:, thr_idx, s));
                
        if (max_support > 0)
            newFamily = id_min;
            newFamily_nEdges = nEdgesForEachView(id_min); 
        else
            if (thr_idx < length(thrs))
                thr_idx = thr_idx + 1;
            else
                s = s-1;
                thr_idx = 1;
            end
        end
        
        if (s > 0)
            continue;
        end

        % Here, s = 0. Do SRA (single rotation averaging)!
        % Let every family member vote for nonfamily neighbors.
        % The one with the most votes gets added to family through SRA.
        votes(isFamily==1) = 0;
        voters = 1:nViews;
        voters = voters(~has_voted & isFamily);
        for i = voters
            for j = edge_info{i}
                if (~isFamily(j))
                    votes(j) = votes(j) + 1;
                end
            end
        end
        has_voted(voters) = 1;


        % Add the one with the most votes to the family.
        [maxVotes, id_to_check] = max(votes);
        if (maxVotes == 0)
            break;
        end
        isFamily(id_to_check) = 1; 
        newFamily = id_to_check;
        newFamily_nEdges = nEdgesForEachView(id_to_check);
        

        % Determine the rotation of the new family member.
        c = 0;
        R_i = cell(1,1);
        for j = edge_info{id_to_check}
            if (isFamily(j))
                c = c + 1;
                R_ij = G(3*id_to_check-2:3*id_to_check, 3*j-2:3*j);
                R_i{c} = R_ij*R_init{j};
            end
        end
        R_init{id_to_check} = ClosestToChordalL1Mean(R_i, true, 10, 0.001);

        % Reset thr and s (since new family member is added).
        s = s_init;
        thr_idx = 1;
    end


    time_initialization = toc(t_init);


    %% Step 2: Edge filtering

    t_opti = tic;
    
    if (med_err > median_thr) 
        edge_IDs_ = edge_IDs;
        RR_ = RR;
    else
        RR_ = zeros(3,3);
        edge_IDs_ = zeros(2,1);

        cc = 0;
        for i = 1:nEdges
            j = edge_IDs(1,i);
            k = edge_IDs(2,i);

            R_jk = RR(:,:,i);
            R_jk_est = R_init{j}*R_init{k}';

            delta_R = R_jk_est-R_jk;
            err_sq = sum(sum(delta_R.^2));

            if (err_sq < err_sq_thr)
                cc = cc + 1;
                RR_(:,:,cc) = RR(:,:,i);
                edge_IDs_(:,cc) = [j;k];
            end
        end
        %disp(['Edges removed = ', num2str((nEdges-cc)/nEdges*100), '%'])
    end

    %% Step 3: Local optimization
    
    nIterations = 250;
    [R_est, ~, ~] = LocalOptimization(R_init, edge_IDs_, RR_, 'L0.5', nIterations);

    time_optimization = toc(t_opti);
end