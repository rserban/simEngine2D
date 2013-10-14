% Read the sliderCrankKIN model
sys = MBsys('Models/sliderCrankKIN.adm');

% Print out information about this model
sys.print();

% Go through all constaints in the system and evaluate the following
% quantities for each one of them: Phi, Jac, nu, gamma
% These quantities are evaluated at t=0, q=q0 and qd = 0
qd0 = zeros(3,1);
flags = [1,1,1,1];
bodies = sys.bodies;

for ic = 1:sys.nC
    this = sys.constraints(ic);
    fprintf('%2i  %-15s\n', ic, this.getType());

    qi = bodies(this.bodyI).q0;

    if isempty(this.bodyJ)
        this.evalAndPrint(0, qi, qd0); 
    else
        qj = bodies(this.bodyJ).q0;
        this.evalAndPrint(0, qi, qj, qd0, qd0); 
    end
end