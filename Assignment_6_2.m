function [] = Assignment_6_2()

% Read ADM file
sys = MBsys('Models/simplePend.adm');
sys.print();

% At t = 0
fprintf('At initial time:\n');
[q0,qd0] = sys.setIC;
Phi = sys.evalPhi(0, q0)
Jac = sys.evalPhi_q(0, q0)
nu = sys.evalNu(0, q0)
gamma = sys.evalGamma(0, q0, qd0)


% Loop on a time grid...
dt = 0.01;
nt = 100;

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
    if it == 37
       disp(Phi)
       disp(Jac)
       disp(nu)
    end
end
