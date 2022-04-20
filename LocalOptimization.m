function [R_est, time_iterations, iterations] = LocalOptimization(R_init, edge_IDs, RR, mode, nIterations)
    
    % For this part, I used Chatterjee's code as a reference:
    % https://ee.iisc.ac.in/labs/cvl/research/rotaveraging/
    

    someVerySmallNumber = 1e-6;
    nEdges = size(edge_IDs, 2);
    nViews = length(R_init);
    
    % Fix the ambiguity (gauge freedom) by not updating the first view.
    nEntriesExceptAnchor = sum(sum(edge_IDs ~= 1));
    row_A = zeros(1, nEntriesExceptAnchor); col_A = row_A; val_A = row_A;
    
    c = 0;
    for i = 1:nEdges
        j = edge_IDs(1,i)-1;
        k = edge_IDs(2,i)-1;
        
        if (j > 0)
            c = c + 1;
            row_A(c) = i; col_A(c) = j; val_A(c) = 1; %A(i, j) = 1 
        end
        
        if (k > 0)
            c = c + 1;
            row_A(c) = i; col_A(c) = k; val_A(c) = -1; %A(i, k) = -1 
        end
    end

    A = sparse(row_A, col_A, val_A); % [nEdges x (nViews-1)]
    
    w = ones(nEdges, 1);
    v = zeros(nViews, 3); % Update in terms of rotation vector.
    
    qjk_qk = zeros(4,nEdges);
    invOfqj_qjk_qk = zeros(4,nEdges);
    qjk_all = zeros(4,nEdges);
    for i = 1:nEdges
        qjk_all(:,i) = R2q(RR(:,:,i));
    end
    q_all = zeros(4,nViews);
    for i = 1:nViews
        q_all(:,i) = R2q(R_init{i});
    end
    q_all_updated = zeros(4,nViews);
    delta_q_all = zeros(4,nViews);
    
    js = edge_IDs(1,:);
    ks = edge_IDs(2,:);
    
    tIterations = tic;
    for it = 1:nIterations
       
        qjk_qk(1,:) = ... %scalar term of (q_jk)(q_k)
            qjk_all(1,:).*q_all(1,ks) - sum(qjk_all(2:4,:).*q_all(2:4,ks),1); 
        
        qjk_qk(2:4,:) = ... % vector term of (q_jk)(q_k)
           qjk_all(1,:).*q_all(2:4,ks) + q_all(1,ks).*qjk_all(2:4,:) ...
           + [qjk_all(3,:).*q_all(4,ks) - qjk_all(4,:).*q_all(3,ks);...
              qjk_all(4,:).*q_all(2,ks) - qjk_all(2,:).*q_all(4,ks);...
              qjk_all(2,:).*q_all(3,ks) - qjk_all(3,:).*q_all(2,ks)];
       
        invOfqj_qjk_qk(1,:) = ... %scalar term of inv(q_j)(q_jk)(q_k)
            -q_all(1,js).*qjk_qk(1,:) - sum(q_all(2:4,js).*qjk_qk(2:4,:),1); 
        
        invOfqj_qjk_qk(2:4,:) = ... %vector term of inv(q_j)(q_jk)(q_k)
           -q_all(1,js).*qjk_qk(2:4,:) + qjk_qk(1,:).*q_all(2:4,js) ...
           + [q_all(3,js).*qjk_qk(4,:) - q_all(4,js).*qjk_qk(3,:);...
              q_all(4,js).*qjk_qk(2,:) - q_all(2,js).*qjk_qk(4,:);...
              q_all(2,js).*qjk_qk(3,:) - q_all(3,js).*qjk_qk(2,:)]; 

        vij_norm = sqrt(sum(invOfqj_qjk_qk(2:4,:).^2, 1));
        vij_theta = 2*atan2(vij_norm, invOfqj_qjk_qk(1,:));
        
        ids_theta_smaller_than_minus_pi = vij_theta < -pi;
        vij_theta(ids_theta_smaller_than_minus_pi) = vij_theta(ids_theta_smaller_than_minus_pi) + 2*pi;
        ids_theta_larger_than_pi = vij_theta > pi;
        vij_theta(ids_theta_larger_than_pi) = vij_theta(ids_theta_larger_than_pi) - 2*pi;
            
        % The three lines below prevent division by near-zero.
        ids_theta_too_small = vij_norm < someVerySmallNumber;
        vij_norm(ids_theta_too_small) = 1;
        vij_theta(ids_theta_too_small) = 0;
        
        B = (vij_theta./vij_norm).*invOfqj_qjk_qk(2:4,:);
        B = B';
                
        WB = w.*B; % [nEdges  x 3]
        W = sparse(1:length(w), 1:length(w), w, length(w), length(w));
        WA =W*A;  % [nEdges x (nViews-1)]
        
        v(2:end,:) = WA\WB; % [(nViews-1) x 3]
        
        residuals_sq = A*v(2:end,:)-B;
        residuals_sq = sum(residuals_sq.^2, 2);
        
        
        
        if (strcmp(mode, 'L1'))
%             w = min(1e4, residuals.^(-1/2)); % sqrt of weight from L1 norm
            w = min(1e4, residuals_sq.^(-1/4)); % sqrt of weight from L1 norm
        elseif (strcmp(mode, 'L0.5'))
%             w = min(1e4, residuals.^(-3/4)); % sqrt of weight from L0.5 norm
            w = min(1e4, residuals_sq.^(-3/8)); % sqrt of weight from L0.5 norm
        end
        
        
        v = v'; %[nViews x 3]
        v_theta = sqrt(sum(v.*v, 1));
        sin_half_theta = sin(v_theta/2);
        cos_half_theta = cos(v_theta/2);
        
        mean_theta = mean(v_theta(2:end));
        
        % The four lines below prevent division by near-zero.
        ids_theta_too_small = v_theta < someVerySmallNumber;
        sin_half_theta(ids_theta_too_small) = 0;
        cos_half_theta(ids_theta_too_small) = 1;
        v_theta(ids_theta_too_small) = 1;
        
        delta_q_all(1,:) = cos_half_theta;
        delta_q_all(2:4,:) = (sin_half_theta./v_theta).*v;
        v = v'; %[3 x nViews]
        
                
        q_all_updated(1,:) = ... %scalar term of q*delta_q
            q_all(1,:).*delta_q_all(1,:) - sum(q_all(2:4,:).*delta_q_all(2:4,:),1); 
        
        q_all_updated(2:4,:) = ... % vector term of q*delta_q
           q_all(1,:).*delta_q_all(2:4,:) + delta_q_all(1,:).*q_all(2:4,:) ...
           + [q_all(3,:).*delta_q_all(4,:) - q_all(4,:).*delta_q_all(3,:);...
              q_all(4,:).*delta_q_all(2,:) - q_all(2,:).*delta_q_all(4,:);...
              q_all(2,:).*delta_q_all(3,:) - q_all(3,:).*delta_q_all(2,:)];
        
        q_all = q_all_updated;  

        if (mean_theta < 1e-3) % Same as Chatterjee's
            break;
        end
    end
    time_iterations = toc(tIterations);
    iterations = it;
    
    R_est = cell(1,nViews);
    for i = 1:nViews
        R_est{i} = q2R(q_all(:,i));
    end
end
