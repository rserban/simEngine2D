function F = pistonForce(t, qi, qdi)
x3 = qi(1);
x3d = qdi(1);

F = [0;0];
if x3d > 0
    if x3 > 5.5
        F(1) = 0;
    elseif x3 > 5
        F(1) = -110000 * ( 1 - sin(2*pi*(x3-5.25)) );
    elseif x3 >= 1.5
        F(1) = 62857 - 282857/(6-x3);
    else
        F(1) = 0;
    end
end
