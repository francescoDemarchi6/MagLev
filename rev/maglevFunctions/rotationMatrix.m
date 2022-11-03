%{
    https://en.wikipedia.org/wiki/Rotation_matrix
    
    Mr =  Mz(gamma) * My(alpha) * Mx(beta) = YAW * PITCH * ROLL where:
        
        Mz(gamma) = [ 1,  0,            0;
                      0,  cos(gamma),   -sin(gamma);
                      0   sin(gamma)    cos(gamma); ];

        My(alpha) = [ cos(alpha),  -sin(alpha),    0;
                      sin(alpha),  cos(alpha),     0;
                      0            0               1 ];

        Mx(beta) = [ cos(beta),     0,    sin(beta);
                     0,             1,    0;
                     -sin(beta),    0,    cos(beta) ];    
    
    4x4 matrix where:
        top-left 3x3 sub-matrix is Mr
        top-left 1x3 sub-matrix represents a translation matrix (0_1x3 => no translation),
        last row 4x1 is always [0 0 0 1]
%}

function Mr = rotationMatrix(alpha,beta,gamma)
    Mr = [
            cos(beta)*cos(gamma),   sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma),    cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma),    0;
            cos(beta)*sin(gamma),   sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),    cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),    0;
            -sin(beta),             sin(alpha)*cos(beta),                                       cos(alpha)*cos(beta),                                       0;
            0,                      0,                                                          0,                                                          1
         ];
end