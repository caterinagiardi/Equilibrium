function [result] = growth_function(input_x, input_l, input_k, n)
%   returns the array with the function caluclated for each node status (x)
%   implements the growth function

result = zeros(n+1, 1);

%   f for index beteween 1 and n

syms x l k;
f(x, l, k) = l * x * (1 - (x / k));

%   f for n+1 index

f_newnode(x, k) = x * (1 - (x / k));

for i = 1:n
    result(i, 1) = f(input_x(i), input_l(i), input_k(i));
    result(n+1, 1) = f_newnode(input_x(i), input_k(i));
end
