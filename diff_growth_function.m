function [result] = diff_growth_function(input_x, input_l, input_k, n)
%   return the array with the derivative of the f function calculated for
%   each node status (x) 
%   implementation of the derivative function of f.
result = zeros(n+1, 1);

%   f for index beteween 1 and n

syms x l k;
f(x, l, k) = l * x * (1 - (x / k));

%   f for n+1 index

f_newnode(x, k) = x * (1 - (x / k));

diff_f = diff(f, x);
diff_f_newnode = diff(f_newnode, x);

for i = 1:n
    result(i, 1) = diff_f(input_x(i), input_l(i), input_k(i));
    result(n+1, 1) = diff_f_newnode(input_x(i), input_k(i));
end

end

