function decompostion()
    n_space = 100;
    n_time = 300;
    Q = solve(n_space, n_time, ones(n_space, 1), eye(n_space));
    subplot(1,3,1)
    plot(Q)
    c = Q*Q';
    [eigvec, lambda] = eig(c);
    diag(lambda)
    
    num_eigvec_keep = 6;
    important_eigvec = zeros(n_space, num_eigvec_keep+1);
    important_eigvec(:, 1) = 0.1*ones(n_space, 1);
    important_eigvec(:, 2:num_eigvec_keep+1) = eigvec(:, length(eigvec)-num_eigvec_keep + 1:length(eigvec));
    norm(important_eigvec(:, 2))

    subplot(1,3,2)
    %plot(important_eigvec)
    semilogy(diag(lambda) / lambda(100,100), 'r*')
    grid on
    initial = zeros(num_eigvec_keep + 1, 1);
    initial(1) = 10;
    reduced = solve(n_space, n_time, initial, important_eigvec)
    
    solutions = important_eigvec * reduced;
    subplot(1,3,3)
    plot(solutions)
end