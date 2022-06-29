function [outputArg1] = equilibrium(M, k, l, n, target)
%EQUILIBRIUM Summary of this function goes here
%   n different species
%   x(i) population density of species i 
%   f(i) growth function for species i
%   g(i, j) interaction functions between species i and j
%   m(i, j) interaction coefficients: j is the prey and i the predator

target_f = zeros(n+1, 1);
target_g = zeros(n+1, 1);
G = zeros(n);
q = zeros(n,1);
p = sym(zeros(n,1));
I = eye(n); 

%   STEP 1: find alpha (together with beta they make up the interaction 
%   coefficient)

alpha = zeros(n+1, 1);

%   first we need the growth function and interaction function calculated 
%   in the target

target_f = growth_function(target, l, k, n);
target_g = interaction_function(M, target, n);

%   now we can proceed with calculating alpha
%   sum contains the sum over rows of the g matrix minus the g(i,j) with 
%   i = j. 

% NB: we need target_g(i,n+1) != 0

for i = 1:n+1
    alpha(i) = (-1 / (target_g(i, n+1))) .* (target_f(i) + sum(target_g(i,:)) - target_g(i, i));
end

%   we also need the derivative of the growth function(f) and the
%   interaction function(g)

diff_target_f = diff_growth_function(target, l, k, n);
diff_target_g = diff_interaction_function(target, M, n);

%    _       _
%   | G  |  q |
%   |____|____|
%   | p  | r  |
%    -       -
%
%   we calculate the matrix G (nxn) : the Jacobian of the system 
%   defined over the original n-node network, as well as p, q and r

sym beta;
beta = sym('beta',[1, n+1]);
beta_array = sym2cell(beta);

r(beta) = beta(n+1) * diff_target_f(n+1);

for i = 1:n

    G(i,i) = diff_target_f(i) + alpha(i) * diff_target_g(i,n+1,1);
    q(i) = alpha(i) * diff_target_g(i, n+1, 2);
    p(i) = beta(i) * diff_target_g(n+1, i, 1);
    r(beta_array{:}) = r(beta_array{:}) + beta(i)*diff_target_g(n+1, i, 1);

    for j = 1:n
        if j ~= i

            G(i,j) = diff_target_g(i,j,2);

            % per implementare la sommatoria scorro su tutti i j della i
            % corrente

            G(i,i) = G(i,i) + diff_target_g(i,j,1);
            
        end
    end
end

syms lambda;

a_G = det(lambda * I - G);

A_G = adjoint(lambda * I - G);

a_J(beta_array{:}) = (lambda - r(beta_array{:})) * a_G - transpose(p) * A_G * q;

% b has a dimension of n+1, it is distributed in 2 variables (b_1_to_n and
% b_newnode)
% d has a dimension of n

B = zeros(n+1, n+1);

sigma = rand(1)*10;
target_a_J = (lambda + sigma)*(lambda + 100 * sigma)^n;
a = fliplr(coeffs(target_a_J, 'All'));
a = a(1:n+1);

% calcolo B e d

d = fliplr(coeffs((lambda * a_G), 'All'));
d = d(1:n+1);

size_d = size(d);
size_a = size(a);

% B structure: b(j,i) is the i-degree coefficients of the polynomials b(j)
% so i goes from 1, that is 0-degree, to n+1, that is the n-degree
% instead j goes from 1, that is the first polynomial, to n+1, that is the
% last one.

% First, we calculate b_j(lambda)

b = sym(zeros(n+1));
s = sym(zeros(n+1));

for j = 1:n
    for i = 1:n
        s(j) = s(j) + (A_G(j,i) * q(i));
    end
    b(j) = - diff_target_g(n+1, j, 1) * a_G - diff_target_g(n+1, j, 2) * s(j);
end
b(n+1) = - diff_target_f(n+1) * a_G;

% Then, we fill the B matrix (it will not contain lambda because it takes
% the coefficient on the polynomial of lambda)

% Declaration of v
v = zeros(n+1);

for j = 1:n+1
   B(:, j) = fliplr(coeffs(b(j), 'All'));
   if j ~= n+1
        v(j) = target_g(n+1,j);
   else
       v(j) = target_f(j);
   end
end

% B matrix is n+1 x n+1 and it contains all the couefficients of the
% polynomials b. signed_B is instead n+1 x n. We then calculate
% signed_beta_star which contains the first n betas. To calculate the last
% beta_newnode we rely on equation (26).

size_B = size(B);
signed_B = B(:,1:n) - B(:,n+1)*v(1:n)/v(n+1);
size_signed_B = size(signed_B);
beta_star = zeros(n);
notStabilized = true;

while(notStabilized)
    notStabilized = false;
    sigma = sigma*0.9;

    % se signed_B è full rank, cioè il numero di colonne è uguale al rango
    if rank(signed_B) == size(signed_B, 2)
        signed_beta_star = (transpose(signed_B)*signed_B) \ (transpose(signed_B) * transpose(a - d));
        beta_newnode =  - v(1:n) * signed_beta_star / v(n+1);
        signed_beta_star = cat(1, signed_beta_star, beta_newnode)
    end
 
    % ricalcolo r e p con i nuovi beta

    r = signed_beta_star(n+1) * diff_target_f(n+1);

    % p(beta, i) = beta(i) * diff_target_g(n+1, i, 1);
    
    for i = 1:n
        p(i) = signed_beta_star(i) * diff_target_g(n+1, i, 1);
        r = r + signed_beta_star(i) * diff_target_g(n+1, i, 1);
    end

    hat_a_J = (lambda - r) * a_G - transpose(p) * A_G * q
    s_hat_a_J = simplify(hat_a_J)
    %tmp = coeffs(hat_a_J, lambda, 'All')
    prva = solve(lambda^2-lambda-6)
    eigenvalues = solve(hat_a_J)
    for i = 1:size(eigenvalues, 2)
        if(((real(eigenvalues(i))) > 0) && (real(eigenvalues(i)) ~= 0))
            notStabilized = true;
        else 
            notStabilized = false;
        end
    end

end

outputArg1 = a_J(lambda ,beta_array{:});

end

