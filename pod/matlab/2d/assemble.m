%----------------------------------------------------------------------
% simple test for model order reduction in 2D
%             Pavel Kůs 
%                2014
%----------------------------------------------------------------------

function [A, f] = assemble(coef_mat, load_mat, elems_per_block)
    n_blocks = length(coef_mat);
    h = 1/(n_blocks * elems_per_block);
    dofs_per_row = n_blocks * elems_per_block - 1;
    ndofs = (dofs_per_row)^2;
    fprintf('assembling problem with %d dofs\n', ndofs)

    local_stiffness = 1/6 * [4, -1, -1, -2;
                             -1, 4, -2, -1;
                             -1, -2, 4, -1;
                             -2, -1, -1, 4];                         
    local_rhs = h*h/4;
    
    ndof_neighbours = 9;
    
    rows = ones(1, ndof_neighbours*ndofs);
    cols = ones(1, ndof_neighbours*ndofs);
    vals = zeros(1, ndof_neighbours*ndofs);
    
    offsets = zeros(2*dofs_per_row + 3, 1);
    offsets(1) = 0;
    offsets(2) = 1;
    offsets(3) = 2;
    offsets(dofs_per_row + 1) = 3;
    offsets(dofs_per_row + 2) = 4;
    offsets(dofs_per_row + 3) = 5;
    offsets(2 * dofs_per_row + 1) = 6;
    offsets(2 * dofs_per_row + 2) = 7;
    offsets(2 * dofs_per_row + 3) = 8;
    
    f = zeros(ndofs, 1);
    element_idx = 1;
    %i rade, j sloupec
    for i_block = 1:n_blocks
        for i_in_block = 1:elems_per_block
            i_total = elems_per_block * (i_block - 1) + i_in_block;
            for j_block = 1:n_blocks
                for j_in_block = 1:elems_per_block 
                    j_total = elems_per_block * (j_block - 1) + j_in_block;
                    assert(element_idx == (elems_per_block*n_blocks)*(i_total-1) + j_total);
                    coef = coef_mat(i_block, j_block);
                    load = load_mat(i_block, j_block);
                    
                    %okolo daneho elementu v tomto poradi:
                    %   1 2
                    %   3 4
                    global_indices = [-1,-1,-1,-1];
                    if(i_total > 1 && j_total > 1) 
                        global_indices(1) = (elems_per_block * n_blocks - 1) * (i_total - 2) + j_total - 1;
                    end
                    if(i_total > 1 && j_total < n_blocks*elems_per_block)
                        global_indices(2) = (elems_per_block * n_blocks - 1) * (i_total - 2) + j_total;
                    end
                    if(i_total < n_blocks*elems_per_block && j_total > 1)
                        global_indices(3) = (elems_per_block * n_blocks - 1) * (i_total - 1) + j_total - 1;
                    end
                    if(i_total < n_blocks*elems_per_block && j_total < n_blocks*elems_per_block)
                        global_indices(4) = (elems_per_block * n_blocks - 1) * (i_total - 1) + j_total;
                    end
                    
                    for i = 1:4
                        if global_indices(i) >= 0
                            row = global_indices(i);
                            offset_start = (row - 1) * ndof_neighbours + 1;
                            for j = 1:4
                                 if global_indices(j) >= 0
                                     col = global_indices(j);
                                     offset = offset_start + offsets(col - row + dofs_per_row + 2);
                                     rows(offset) = row;
                                     cols(offset) = col;
                                     vals(offset) = vals(offset) + coef * local_stiffness(i,j);
                                     
                                 end
                            end
                           
                            f(global_indices(i)) = f(global_indices(i)) + load * local_rhs;
                        end
                        
                    end
                    
                    element_idx = element_idx + 1;
                end
            end
        end
    end

    A = sparse(rows, cols, vals);
end