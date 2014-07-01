function extract_modes(what_changes)
    n_blocks = 4;
    n_per_block = 10;
    
    dofs = (n_blocks * n_per_block - 1)^2;
%    Q = zeros(dofs, n_blocks^2);
    
    idx = 1;
    for i = 1:n_blocks
        for j = 1:n_blocks
            for i2 = 1:n_blocks
                for j2 = 1:n_blocks
                    coef_mat = ones(n_blocks);
                    load_mat = ones(n_blocks);
                    changing_mat = ones(n_blocks);
                    changing_mat(i,j) = 5;
                    changing_mat(i2,j2) = 5;
                    if strcmp(what_changes, 'load')
                        load_mat = changing_mat;
                    else
                        coef_mat = changing_mat;
                    end
                    [A,f] = assemble(coef_mat, load_mat, n_per_block);
                    x = A\f;            
                
                    Q(:, idx) = x;            
                    idx = idx + 1;
                end
            end
        end
    end
    
    c = Q*Q';
    
    use_modes = 50;
    show_modes = round(1.2*use_modes);

    [vec, lambda] = eigs(c, show_modes);
    lambda = diag(lambda);
    for i = 1:show_modes
        subplot(1, 2, 1)
        semilogy(1:show_modes, lambda, '*r', i, lambda(i), '*b');
        subplot(1, 2, 2)
        show_solution(vec(:, i));
%        waitforbuttonpress;
            
    end
       
    modes = vec(:, 1:use_modes);
    
    coef_mat = ones(n_blocks);
    load_mat = ones(n_blocks);
    changing_mat = ones(n_blocks);

    changing_mat(2,2) = 5;
    changing_mat(2,3) = 5;
    changing_mat(3,3) = 3;
    changing_mat(4,3) = 2;
    changing_mat(4,4) = 6;
    if strcmp(what_changes, 'load')
        load_mat = changing_mat;
    else
        coef_mat = changing_mat;
    end
    
    [A,f] = assemble(coef_mat, load_mat, n_per_block);
    
    x1 = A\f;
    reduced_matrix = modes'*A*modes;
    reduced_rhs = modes'*f;
    xi2 = reduced_matrix\reduced_rhs;
    x2 = modes * xi2;
    fprintf('original dofs %d, reduced dofs %d, difference %d\n', length(A), length(reduced_matrix), norm(x1-x2)/norm(x1));
    subplot(1, 2, 1)
    show_solution(x1);
    subplot(1, 2, 2)
    show_solution(x2);
    
end

