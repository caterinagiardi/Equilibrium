function [result] = interaction_function(M, input_x, n)
%   return the interaction function array with results calculated in x
%   implementation of the interaction function

result = zeros(n+1, n+1);

%   g for index beteween 1 and n

syms x y m;
g(x, y, m) = m * x * y;

%   g for the other indexes

g_newnode(x, y) = x * y;

for i = 1:n
    for j = 1:n
        result(i,j) =  g(input_x(i), input_x(j), M(i,j));
    end
    result(i, n+1) = g_newnode(input_x(i), input_x(n+1));
    result(n+1, i) = g_newnode(input_x(i), input_x(n+1));
end
result(n+1, n+1) = g_newnode(input_x(n+1), input_x(n+1));

end