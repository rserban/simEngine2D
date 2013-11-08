function data = dynamics(sys, tstart, tend, dt, nOut, tol, maxIt)
%DYNAMICS - ...
%     data = dynamics(sys, tstart, tend, dt, nOut)

% Use default tolerance and maximum number of iterations if not specified.
if nargin < 6 || isempty(tol)
    tol = 1e-8;
end
if tol < eps
    warning('Specified tolerance too small. Reset at %g', eps);
    tol = eps;
end
if nargin < 7 || isempty(maxIt)
    maxIt = 100;
end

% Calculate output time grid.
tout = linspace(tstart,tend,nOut+1);

% Preallocate output (one additional output at tstart)
data.t = tout;
data.q = zeros(sys.n,nOut+1);
data.qd = zeros(sys.n,nOut+1);
data.qdd = zeros(sys.n,nOut+1);
data.lam = zeros(sys.m,nOut+1);

%% Initialization.
t = tstart;

% Extract initial conditions. Check position and velocity constraints.
[q,qd] = sys.setIC;
Phi = sys.evalPhi(t, q);
Phi_q = sys.evalPhi_q(t, q);
Nu = sys.evalNu(t, q);

nrm = norm(Phi, inf);
if nrm > tol
    warning('Max. position constraint violation at initial time: %g', nrm);
end
nrm = norm(Phi_q*qd-Nu, inf);
if nrm > tol
    warning('Max. velocity constraint violation at initial time: %g', nrm);
end

% Calculate initial values for accelerations and Lagrange multipliers.
M = sys.getM;
Q = sys.evalQ(t, q, qd);
Gamma = sys.evalGamma(t, q, qd);
Mbar = [M , Phi_q' ; Phi_q , zeros(sys.m)];
x = Mbar\[Q;Gamma];
qdd = x(1:sys.n);
lam = x(sys.n+1:end);

% Set parameters in Newmark
beta = 0.3025;
gam = 0.6;

% Store data at tstart
nextOut = 1;
data.q(:,nextOut) = q;
data.qd(:,nextOut) = qd;
data.qdd(:,nextOut) = qdd;
data.lam(:,nextOut) = lam;
nextOut = nextOut + 1;

%% Loop until we reach the final time.
hw = waitbar(0,'Initializing waitbar...');

while true
    % Increment time, adjusting to hit exactly the output times.
    t_str = sprintf('t = %.4f', t);
    waitbar(t/tend,hw,t_str)
    
    t_next = t + dt;
    if abs(t_next - tout(nextOut)) < dt/10
        h = tout(nextOut) - t;
        t = tout(nextOut);
    else
        h = dt;
        t = t + dt;
    end
        
    %t = t + h;
    inv_beta = 1/(beta * h^2);
    
    q_prev = q;
    qd_prev = qd;
    qdd_prev = qdd;
    
    % Iterate to find solution (i.e. 'qdd' and 'lam') at new time.
    % We use as initial guess the solution from the previous time.
    for iter = 1:maxIt
        % Update position and velocity (using Newmark)
        q = q_prev + h*qd_prev + 0.5*h^2*((1-2*beta)*qdd_prev + 2*beta*qdd);
        qd = qd_prev + h*((1-gam)*qdd_prev + gam*qdd);
        
        % Evaluate residual and (simplified) Jacobian.
        Phi = sys.evalPhi(t, q);
        Q = sys.evalQ(t, q, qd);
        Phi_q = sys.evalPhi_q(t, q);
        
        R = [M*qdd + Phi_q' * lam - Q;
            inv_beta * Phi];
        J = [M , Phi_q' ; Phi_q , zeros(sys.m)];
        
        % Solve for corrections and update solution.
        del = J\R;
        qdd = qdd - del(1:sys.n);
        lam = lam - del(sys.n+1:end);
        
        % Check for convergence (exclude Lagrange multipliers).
        if norm(del(1:sys.n), inf) <= tol
            break;
        end
    end
    
    
    %fprintf('t = %g  h = %g  iter = %i  nrm  %g\n', t, h, iter, norm(del,inf));
    
    
    
    % If this is an output time, store data.
    if t >= tout(nextOut)
        %fprintf('output at t = %g\n', t);
        data.q(:,nextOut) = q_prev + h*qd_prev + 0.5*h^2*((1-2*beta)*qdd_prev + 2*beta*qdd);
        data.qd(:,nextOut) = qd_prev + h*((1-gam)*qdd_prev + gam*qdd);
        data.qdd(:,nextOut) = qdd;
        data.lam(:,nextOut) = lam;
        
        nextOut = nextOut+1;
    end
    
    % Stop now if we reached the final time.
    if t >= tend
        break;
    end
    
end

close(hw);

end

