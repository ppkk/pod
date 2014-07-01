function show_solution(x)
    rows = sqrt(length(x));
    x_mat = vec2mat(x, rows, rows);
    
    %fprintf('maximal value %d\n', max(max(x_mat)));
    
    x_mat = smooth(x_mat);
    x_mat = smooth(x_mat);
    x_mat = smooth(x_mat);
    imagesc(x_mat);   
end

function x = smooth(a)
    len = length(a);
   
    vec = zeros(1, 2*len-1);
    vec(1) = 1;
    vec(2) = 0.5;
    smoother = toeplitz(vec);
    select = 2*(1:len)-1;
    smoother = smoother(select, :);
    x = smoother' * a * smoother;
    
end