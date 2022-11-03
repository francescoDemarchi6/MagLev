%{
    Same fun implemented in teensy code. PLEASE NOTE: array indices in
    matlab starts with 1 so instead of ux(0,1,2) we have ux(1,2,3)
%}

function u = PIDfun(y,PIDparams)
    persistent ex ey ez ux uy uz;
    
    if isempty(ux)
        ux = zeros(1,2); uy = ux; uz = ux;
        ex = zeros(1,3); ey = ex; ez = ex;
    end

    KPx = PIDparams.KPx; KIx = PIDparams.KIx; KDx = PIDparams.KDx;
    KPy = PIDparams.KPy; KIy = PIDparams.KIy; KDy = PIDparams.KDy;
    KPz = PIDparams.KPz; KIz = PIDparams.KIz; KDz = PIDparams.KDz;
    
    %setpoints? as for the potentiometers?
    rx = PIDparams.rx; ry = PIDparams.ry; rz = PIDparams.rz;
   
    % Inputs x, y and z readings
    yx = y(1);
    yy = y(5);
    yz = y(9);

    % Measurement error
    ex(1) = rx - yx;
    ey(1) = ry - yy;
    ez(1) = rz - yz;

    % Compute gain
    ux(1) = ux(1) + KPx*(ex(1) - ex(2)) + KIx*ex(1) + KDx*(ex(1) - 2*ex(2) + ex(3));
    uy(1) = uy(1) + KPy*(ey(1) - ey(2)) + KIy*ey(1) + KDy*(ey(1) - 2*ey(2) + ey(3));
    uz(1) = uz(1) + KPz*(ez(1) - ez(2)) + KIz*ez(1) + KDz*(ez(1) - 2*ez(2) + ez(3));

    % Update values
    circshift(ex,1);
    circshift(ey,1);
    circshift(ez,1);

    circshift(ux,1);
    circshift(uy,1);
    circshift(uz,1);

    ux(1) = min(max(-255, ux(1)), 255);
    uy(1) = min(max(-255, uy(1)), 255);
    uz(1) = min(max(-100, uz(1)), 100);

    %% TURNX/Y
    sumx1 = +ux(1) +uz(1);
    sumx2 = -ux(1) +uz(1);

    sumy1 = +uy(1) +uz(1);
    sumy2 = -uy(1) +uz(1);

    sumx1 = min(max(-255, sumx1), 255);
    sumx2 = min(max(-255, sumx2), 255);

    sumy1 = min(max(-255, sumy1), 255);
    sumy2 = min(max(-255, sumy2), 255);
    
    %% Gain to current in solenoids (multiplied by 0.5/255):
    sumx1 = sumx1* 0.5/255;
    sumx2 = sumx2* 0.5/255;
    sumy1 = sumy1* 0.5/255;
    sumy2 = sumy2* 0.5/255;

    u = [sumx1 sumx2 sumy1 sumy2];
end