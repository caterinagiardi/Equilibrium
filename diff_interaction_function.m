function [result] = diff_interaction_function(input_x, M, n)
%   return the array with the derivative respect to to the first and the   
%   second of the g function calculated for each node status (x)  
%   implementation of the partial derivative function of g

result = zeros(n+1, n+1, 2);

%   g for index beteween 1 and n

syms x y m;
g(x, y, m) = m * x * y;

%   g for the other indexes

g_newnode(x, y) = x * y;


diff_g_wrt_x = diff(g, x);
diff_g_wrt_y = diff(g, y);

diff_g_newnode_wrt_x = diff(g_newnode, x);
diff_g_newnode_wrt_y = diff(g_newnode, y);

for i = 1:n
    for j = 1:n
        result(i, j, 1) =  diff_g_wrt_x(input_x(i), input_x(j), M(i,j));
        result(i, j, 2) =  diff_g_wrt_y(input_x(i), input_x(j), M(i,j));
    end
    result(i, n+1, 1) = diff_g_newnode_wrt_x(input_x(i), input_x(n+1));
    result(n+1, i, 1) = result(i, n+1, 1);

    result(i, n+1, 2) = diff_g_newnode_wrt_y(input_x(i), input_x(n+1));
    result(n+1, i, 2) = result(i, n+1, 2);
end
end


