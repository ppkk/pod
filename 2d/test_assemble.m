function test_assemble
    n_blocks = 4;
    n_per_block = 60;
    coef_mat = ones(n_blocks);
    load_mat = ones(n_blocks);
    %coef_mat(2,2) = 5;
    %coef_mat(3,4) = 5;
    load_mat = [1,1,1,1;
                1,5,1,1;
                1,5,1,1;
                1,1,5,1];
            
    
    rows = n_blocks * n_per_block - 1;
    
    [A,f] = assemble(coef_mat, load_mat, n_per_block);
    fprintf('assembled\n')
%    imagesc(A);
%    full(A)
    %plot(f);
    %spy(A)
    x = A\f;
    show_solution(x);
end

