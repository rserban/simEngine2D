function [] = Assignment_6_1()

max_it = 20;
tol = 1e-8;

% Initial guess
q = [1.5; 1.5];

% Perform at most max_it N-R iterations
for it = 1:max_it
    res = calcPhi(q);
    jac = calcJac(q);
    corr = jac\res;
    q = q - corr;
    norm_corr = norm(corr);
    
    fprintf('Iteration: %i  ||correction|| = %e  q = [%f  %f]\n', it, norm_corr, q);
    
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