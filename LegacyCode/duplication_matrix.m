% See Magnus and Neudecker (1988), Matrix differential calculus with applications in statistics and econometrics.
% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Created: 8 May 1995
function d = duplication_matrix(n)

d = zeros (n * n, n * (n + 1) / 2);
count = 0;
for j = 1 : n
    d ((j - 1) * n + j, count + j) = 1;
    for i = (j + 1) : n
        d ((j - 1) * n + i, count + i) = 1;
        d ((i - 1) * n + j, count + i) = 1;
    end
    count = count + n - j;
end

end