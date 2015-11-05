function decompostion()
    n_space = 100;
    n_time = 300;
    Q = solve(n_space, n_time, ones(n_space, 1), eye(n_space), false);
    c = Q*Q';
    [eigvec, lambda] = eig(c);
    diag(lambda)
    
    num_eigvec_keep = 3;
    important_eigvec = zeros(n_space, num_eigvec_keep+1);
    important_eigvec(:, 1) = 0.1*ones(n_space, 1);
    important_eigvec(:, 2:num_eigvec_keep+1) = eigvec(:, length(eigvec)-num_eigvec_keep + 1:length(eigvec));
    norm(important_eigvec(:, 2))

    subplot(1,2,1)
    lambda = diag(lambda);
    lambda = max(lambda, 10e-15);
    semilogy(lambda, 'r*')
    grid on
    subplot(1,2,2)
    plot(important_eigvec)
    waitforbuttonpress;

    initial = zeros(num_eigvec_keep + 1, 1);
    initial(1) = 10;
    reduced = solve(n_space, n_time, initial, important_eigvec, false);
    
    subplot(1,2,1)
    plot(Q)
    solutions = important_eigvec * reduced;
    subplot(1,2,2)
    plot(solutions)

    waitforbuttonpress;
    Q = solve(n_space, n_time, ones(n_space, 1), eye(n_space), true);
    subplot(1,2,1)
    plot(Q)
    initial = zeros(num_eigvec_keep + 1, 1);
    initial(1) = 10;
    reduced = solve(n_space, n_time, initial, important_eigvec, true);
    solutions = important_eigvec * reduced;
    subplot(1,2,2)
    plot(solutions)

end