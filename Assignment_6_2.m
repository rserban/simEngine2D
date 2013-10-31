function [] = Assignment_6_2()

% Read ADM file
sys = MBsys('Models/simplePend.adm');
sys.print();

% At t = 0
fprintf('\n\nAt initial time:\n');
[q0,qd0] = sys.setIC;
Phi = sys.evalPhi(0, q0);
Jac = sys.evalPhi_q(0, q0);
nu = sys.evalNu(0, q0);
gamma = sys.evalGamma(0, q0, qd0);

fprintf('Phi   = [%10.6f\n         %10.6f\n         %10.6f]\n', Phi);
fprintf('Phi_q = [%10.6f  %10.6f  %10.6f]\n', Jac(1,:));
fprintf('        [%10.6f  %10.6f  %10.6f]\n', Jac(2,:));
fprintf('        [%10.6f  %10.6f  %10.6f]\n', Jac(3,:));
fprintf('nu    = [%10.6f\n         %10.6f\n         %10.6f]\n', nu);
fprintf('gamma = [%10.6f\n         %10.6f\n         %10.6f]\n', gamma);


% Loop on a time grid...
dt = 0.01;
nt = 100;

fprintf('---------------------------\n');
fprintf('At t = 0.37\n');
for it = 1:nt
    t = dt*it;
    % Consistent q at current time 
    phi = pi/2 + 2*pi*t;
    x = 1 - cos(phi);
    y = 1 - sin(phi);
    q = [x;y;phi];
    % Evaluate Phi, Jac, Nu at current time
    Phi = sys.evalPhi(t, q);
    Jac = sys.evalPhi_q(t, q);
    nu = sys.evalNu(t, q);
    qd = Jac\nu;
    gamma = sys.evalGamma(0, q, qd);
    qdd = Jac\gamma;
    if it == 37
        fprintf('q     = [%10.6f\n         %10.6f\n         %10.6f]\n', q);
        fprintf('qd    = [%10.6f\n         %10.6f\n         %10.6f]\n', qd);
        fprintf('qdd   = [%10.6f\n         %10.6f\n         %10.6f]\n', qdd);
        fprintf('Phi   = [%10.6f\n         %10.6f\n         %10.6f]\n', Phi);
        fprintf('Phi_q = [%10.6f  %10.6f  %10.6f]\n', Jac(1,:));
        fprintf('        [%10.6f  %10.6f  %10.6f]\n', Jac(2,:));
        fprintf('        [%10.6f  %10.6f  %10.6f]\n', Jac(3,:));
        fprintf('nu    = [%10.6f\n         %10.6f\n         %10.6f]\n', nu);
        fprintf('gamma = [%10.6f\n         %10.6f\n         %10.6f]\n', gamma);
    end
end
