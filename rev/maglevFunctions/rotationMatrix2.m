%{
    Look at rotationMatrix.m for explanation.
    This is the same rotation matrix without YAW/Mz(gamma) used because of
    the simmetry around z-axis
%}

function Mr2 = rotationMatrix2(alpha,beta)
    Mr2 = [
            cos(beta),  sin(alpha)*sin(beta),   cos(alpha)*sin(beta),   0;
            0,          cos(alpha),             -sin(alpha),            0;
            -sin(beta), sin(alpha)*cos(beta),   cos(alpha)*cos(beta),   0;
            0,          0,                      0,                      1
          ];
end