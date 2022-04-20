function [newFamily, newFamily_nEdges] = AddToNewFamily(id_to_add, nEdges_of_id_to_add, newFamily, newFamily_nEdges)
    
    % Add to newFamily such that the one with fewer edges comes at the end.
    
    if (isempty(newFamily))
        newFamily = id_to_add;
        newFamily_nEdges = nEdges_of_id_to_add;
        return;
    end
    
    if (nEdges_of_id_to_add <= newFamily_nEdges(end))
        newFamily(end+1) = id_to_add;
        newFamily_nEdges(end+1) = nEdges_of_id_to_add;
        return;
    end

    newFamily_ = [newFamily 0];
    newFamily_nEdges_ = [newFamily_nEdges 0];
    
    for i = 1: length(newFamily)
        nE = newFamily_nEdges(i);
        if (nEdges_of_id_to_add >= nE)
            newFamily_(i) = id_to_add;
            newFamily_nEdges_(i) = nEdges_of_id_to_add;
            newFamily_(i+1:end) = newFamily(i:end);
            newFamily_nEdges_(i+1:end) = newFamily_nEdges(i:end);
            break;
        end
    end
    
    newFamily = newFamily_;
    newFamily_nEdges = newFamily_nEdges_;
end