function data = kinematics(sys, tstart, tend, dt, nOut, tol, maxIt)
%KINEMATICS - perform kinematic analysis for the specified mechanism
%     data = KINEMATICS(SYS, TSTART, TEND, DT, NOUT) performs kinematic
%     analysis for the mechanism represented by SYS over the time interval
%     [TSTART, TEND], using a step-size DT and storing results at NOUT
%     equally-spaced time points.
%     data = KINEMATICS(SYS, TSTART, TEND, DT, NOUT, TOL, MAXIT) also
%     specifies the tolerance and maximum number of iterations for the
%     Newton method used for position analysis. Default values are
%     TOL = 1e-6 and MAXIT = 100.
%
%     SYS is an object of type MBsys, constructed from an ADM file. See
%     MBsys for details.
%
%     The simulation results are returned in the structure DATA which has
%     the following fields:
%        t   - output time points: a (1 x NT) row vector
%        q   - body generalized coordinates: a (3*NB x NT) matrix
%        qd  - body generalized velocities: a (3*NB x NT) matrix
%        qdd - body generalized accelerations: a (3*NB x NT) matrix
%     where NB is the number of bodies in the system and M is the number of
%     constraint equations.
%     The DATA structure can be passed to the varioius reporting methods of
%     the SYS object to display results.
%
%     See also MBSYS, DYNAMICS

% Use default tolerance and maximum number of iterations if not specified.
if nargin < 6 || isempty(tol)
    tol = 1e-6;
end
if tol < eps
    warning('Specified tolerance too small. Reset at %g', eps);
    tol = eps;
end
if nargin < 7 || isempty(maxIt)
    maxIt = 100;
end

% Calculate output frequency and output time grid.
tout = linspace(tstart,tend,nOut+1);

% Preallocate output (one additional output at tstart)
data.t = tout;
data.q = zeros(sys.n,nOut+1);
data.qd = zeros(sys.n,nOut+1);
data.qdd = zeros(sys.n,nOut+1);
data.stats = zeros(1,nOut+1);

% Initialize time, set initial guess at start time, and evaluate Jacobian.
t = tstart;
[q,~] = sys.getIC;
jac = sys.evalPhi_q(t, q);
nextOut = 1;

% Loop until we reach the final time.
while true
    
    % Newton iterations for position analysis.
    % We use a modified Newton method, with the Jacobian evaluated at
    % the initial guess for the current time (which is the solution at
    % the previous step).
    % The iteration stops when the correction is less than the specified
    % tolerance or when the maximum number of iteration was reached.
    
    for iter = 1:maxIt
        res = sys.evalPhi(t, q);
        delta = jac \ res;
        q = q - delta;
        if norm(delta) <= tol
            break;
        end
    end
    %fprintf('  solved at %g\n', t);
    
    % Evaluate Jacobian at converged states.
    jac = sys.evalPhi_q(t, q);
    
    % Collect data if we reached the next output time.
    if t >= tout(nextOut)
        % Perform velocity and acceleration analysis.
        Nu = sys.evalNu(t, q);
        qd = jac \ Nu;
        Gamma = sys.evalGamma(t, q, qd);
        qdd = jac \ Gamma;
        
        % Collect output data.
        data.q(:, nextOut) = q;
        data.qd(:, nextOut) = qd;
        data.qdd(:, nextOut) = qdd;
        data.stats(nextOut) = iter;
        %fprintf('output %i at %g\n', nextOut, t);
        
        nextOut = nextOut + 1;
    end
    
    % Stop now if we reached the final time.
    if t >= tend
        break;
    end
    
    % Advance time, making sure we hit the next output time.
    t = t + min(dt, tout(nextOut)-t);
    
end

end

