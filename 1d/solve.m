function solutions = solve(n_space, n_time, initial, B)
    lambda = 0.01;
    time_total = 30;
    dt = time_total / n_time;
    K = stiffness(n_space);
    M = mass(n_space);
    A = lambda * K + M/dt;
    
    size_B = size(B);
    solutions = zeros(size_B(2), n_time);
    solutions(:, 1) = initial;
    rhs = zeros(n_space, 1);
    time = 0;
    K = B'*A*B;
    M = B'*(M/dt)*B;
    for(m=2:n_time)
        time = time + dt;
        rhs(1) = lambda*derivative(time);
        q = B'*rhs;
        solutions(:, m) = K\ ( M * solutions(:, m-1) + q);
    end
    
    solutions = solutions(:, [1,10*(1:30)]);
end

function q = derivative(time)
    if(time <= 10)
        q = 1;
    else
        q = 0;
    end
end

function K = stiffness(n)
    vec = zeros(n,1);
    vec(1) = 2;
    vec(2) = -1;
    K = toeplitz(vec);
    K(1,1) = 1;
    K(n,n) = 1;
    K = K*n;
end

function M = mass(n)
    vec = zeros(n,1);
    vec(1) = 4;
    vec(2) = 1;
    M = toeplitz(vec);
    M(1,1) = 2;
    M(n,n) = 2;
    M = M/(6*n);

end