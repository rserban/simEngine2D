function [] = Assignment_6_1()
% Sample program for solving the system of nonlinear equations:
%    3x + sin(xy) - 4 = 0
%    x - exp(cos(y)) = 0 
% using Newton's method, with a tolerance 1e-8, using as an initial guess
%    x0 = 1.5
%    y0 = 1.5
%

% Exact solution:
q_an = [1; pi/2];

% Specify maximum number of N-R iterations and stopping tolerance
max_it = 20;
tol = 1e-8;

% Initial guess
q = [1.5; 1.5];

fprintf('Initial guess q = [%20.16f %20.16f]  ||error|| = %20.16f\n', ...
        q, norm(q - q_an));

% Perform at most max_it N-R iterations
for it = 1:max_it
    % Calculate residual and Jacobian at current iterate
    res = calcPhi(q);
    jac = calcJac(q);
    
    % Calculate the correction
    corr = jac\res;
    norm_corr = norm(corr);
    q = q - corr;   
    
    % Output information at current iteration
    fprintf('Iteration: %i  q = [%20.16f %20.16f]  ||error|| = %20.16f',...
            it, q,  norm(q - q_an));
    fprintf('   ||correction|| = %e\n', norm_corr);
        
    % Stop if the norm of the correction is below the specified tolerance.
    if norm_corr <= tol
        break;
    end
end


end


function Phi = calcPhi(q)
Phi = [3 * q(1) + sin(q(1)*q(2)) - 4; 
       q(1) - exp(cos(q(2)))];
end

function Jac = calcJac(q)
x = q(1);
y = q(2);

Jac = [3 + y*cos(x*y) , x*cos(x*y)
       1              , sin(y)*exp(cos(y))];
end