%----------------------------------------------------------------------
% simple test for model order reduction in 2D
%             Pavel Kůs 
%                2014
%----------------------------------------------------------------------

% the problem domain is divided into 4X4 blocks
% in each block, different coefficient (which affects matrix) or load 
% (which affects rhs) is prescribed
% parameter 'coef' or 'load' determines which of the two apply

% the problem is solved by solving a series of simmilar problems with
% different distribution of load or coef and usign MOR

% for load it works exactly (if sufficient number of eigenvalues is used), 
% for coef not that good 
% the reason for the first one is probably the linear dependence of
% solution on rhs (but not on matrix!)

function extract_modes(what_changes)
    n_blocks = 4;
    n_per_block = 10;
    
    dofs = (n_blocks * n_per_block - 1)^2;
%    Q = zeros(dofs, n_blocks^2);
    
    %create a set of problems used as basis for MOR 
    idx = 1;
    for i = 1:n_blocks
        for j = 1:n_blocks
%            uncoment this (and correspondent end) and modify line with changing_mat to get larger basis
%            for i2 = 1:n_blocks
                for j2 = 1:n_blocks
                    coef_mat = ones(n_blocks);
                    load_mat = ones(n_blocks);
                    changing_mat = ones(n_blocks);
                    changing_mat(i,j) = 5;
                    changing_mat(i,j2) = 5;
%                    changing_mat(i2,j2) = 5;
                    if strcmp(what_changes, 'load')
                        load_mat = changing_mat;
                    elseif strcmp(what_changes, 'coef')
                        coef_mat = changing_mat;
                    else
                        assert(false);
                        
                    end
                    [A,f] = assemble(coef_mat, load_mat, n_per_block);
                    x = A\f;            
                
                    Q(:, idx) = x;            
                    idx = idx + 1;
                end
  %          end
        end
    end
    
    c = Q*Q';
    
    use_modes = 30;
    show_modes = round(1.2*use_modes);

    [vec, lambda] = eigs(c, show_modes);
    lambda = diag(lambda);
    lambda = max(lambda, 10e-15);
    for i = 1:show_modes
        subplot(1, 2, 1)
        semilogy(1:show_modes, lambda, '*r', i, lambda(i), '*b');
        subplot(1, 2, 2)
        show_solution(vec(:, i));
        waitforbuttonpress;
            
    end
       
    modes = vec(:, 1:use_modes);
    
    coef_mat = ones(n_blocks);
    load_mat = ones(n_blocks);
    changing_mat = ones(n_blocks);

    %this is a distribution of coef or load which I want to solve
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

