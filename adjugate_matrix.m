function [B] = adjugate_matrix(A)
%   function that calculates the adjugate of a matrix

B = transpose(A);
[n, m] = size(B);

for i = 1:n
    for j = 1:m
        B(i,j) = conj(B(i,j));
    end
end

end
