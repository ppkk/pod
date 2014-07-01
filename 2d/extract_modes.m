function extract_modes
    n_blocks = 4;
    n_per_block = 20;
    
    dofs = (n_blocks * n_per_block - 1)^2;
%    Q = zeros(dofs, n_blocks^2);
    
    load_mat = ones(n_blocks);
    idx = 1;
    for i = 1:n_blocks
        for j = 1:n_blocks
            for j2 = 1:n_blocks
                coef_mat = ones(n_blocks);
                coef_mat(i,j) = 5;
                coef_mat(i,j2) = 5;
                [A,f] = assemble(coef_mat, load_mat, n_per_block);
                x = A\f;            
            
                Q(:, idx) = x;            
                idx = idx + 1;
            end
        end
    end
    
    c = Q*Q';
    
    num_modes = 30;
    [vec, lambda] = eigs(c, num_modes);
    lambda = diag(lambda);
    for i = 1:num_modes
        subplot(1, 2, 1)
        semilogy(1:num_modes, lambda, '*r', i, lambda(i), '*b');
        subplot(1, 2, 2)
        show_solution(vec(:, i));
        waitforbuttonpress;
            
    end
end

