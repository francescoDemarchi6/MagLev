%{ 
    Linearizes the magnetic levitation system in every point of a cube and 
    calculates the condition number 
%}

clear; close all;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat');

approximationType = input("approxType [0/1]> ");
eccentricity = input("Eccentricity [0:1]> ");
eccentricity = max (eccentricity, 0); eccentricity = min (eccentricity, 1); 

%% Searching parameters 
steps = 150; % odd numbers will intersecate the equilibria in (0, 0), resulting in +infinity
L = .10;
bsqr = 1 / (1 - eccentricity^2);

if(approximationType == 0)
    eq = results.zeq.zeq_fst;
    params.magnets.I = results.neo_vs_neo.curr_fst;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_fst;
else
    eq = results.zeq.zeq_acc;
    params.magnets.I = results.neo_vs_neo.curr_acc;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_acc;
end

%% Area to research
Ptx = linspace(-L/2, L/2, steps);
Pty = linspace(-L/2*bsqr, L/2*bsqr, steps); % assuming the ellipse's major axis is y
Zs = linspace(eq-L*9/20, eq+L*1/10, steps);
Cns = zeros(length(Ptx), length(Pty), length(Zs),1); 

%% parameters for the linearization of system
x0 = zeros(12,1); x0(3) = eq;
sys = maglevSystem(x0, params, approximationType, eccentricity);

uLp = zeros(params.solenoids.N,1); % input (it doesn't affect the condition number, but it's needed to linearize)

delta = 1e-4;
dimX = 12;
dimU = params.solenoids.N; % 4

h = waitbar(0);
max = 3*steps;

B = zeros(dimX,dimU); % 12x4

%% Plane XY
%calculationg k in this way can only be done if the linearization has
%already been completed
%[~, k] = min(Cns(ceil(length(Ptx)/2), ceil(length(Pty)/2), :)); % k is the Zs index of the plane with the lowest cn in the middle
%calculating k in this way gives an approximation of the nearest k to the
%equilibrium
k = round(length(Zs)*9/11); % chosen height of the plane to plot: around the equilibrium

fprintf("k: %f  Zs(k): %f\n", k, Zs(k));
planeXY = zeros(length(Ptx), length(Pty)); % contains the condition numbers of the chosen XY plane 
for i = 1:length(Ptx) %x
    for j = 1:length(Pty) %y
        % linearizing
        x0(1) = Ptx(i); x0(2) = Pty(j); x0(3) = Zs(k); % linearizing
        for l = 1:dimU % linearizing around this point
            B(:,l) = (sys.f(x0,uLp+(l==1:dimU)'*delta) ...
                     -sys.f(x0,uLp-(l==1:dimU)'*delta)) ...
                     /(2*delta);               
        end
        temp = B(7:11, :); % coping only the non 0 rows
        Cns(i, j, k) = cond(pinv(temp)); % saving the condition number of the Moore-Penrose pseudinverse
        % saving the CN
        planeXY(i, j) = Cns(i, j, k);
        if planeXY(i, j) > 5e3 % for some reason max() gives an error
            if Cns(i, j, k) > 1e14 && ~(Ptx(i) == 0 && Pty(j) == 0)
                    fprintf("possible equilibrium at x: %f y: %f z: %f, or i: %f j: %f k: %f\n", Ptx(i), Pty(j), Zs(k), i, j, k);
            end
            planeXY(i, j) = 5e3; % saturating the results for better plotting
        end
    end
    waitbar(i/max)
end

%% Plane XZ
j = ceil(length(Pty)/2); % y = 0
fprintf("j: %f  Pty(j): %f\n", j, Pty(j));
planeXZ = zeros(length(Ptx), length(Zs)); % contains the condition numbers of the chosen XZ plane 
for i = 1:length(Ptx) %x
    for k = 1:length(Zs) %z
        % linearizing
        x0(1) = Ptx(i); x0(2) = Pty(j); x0(3) = Zs(k); % linearizing
        for l = 1:dimU % linearizing around this point
            B(:,l) = (sys.f(x0,uLp+(l==1:dimU)'*delta) ...
                     -sys.f(x0,uLp-(l==1:dimU)'*delta)) ...
                     /(2*delta);
        end
        temp = B(7:11, :); % coping only the non 0 rows
        Cns(i, j, k) = cond(pinv(temp)); % saving the condition number of the Moore-Penrose pseudinverse
        % saving the CN
        planeXZ(i, k) = Cns(i, j, k);
        if planeXZ(i, k) > 5e3 % for some reason max() gives an error
            if Cns(i, j, k) > 1e14 && Ptx(i) ~= 0
                    fprintf("possible equilibrium at x: %f y: %f z: %f, or i: %f j: %f k: %f\n", Ptx(i), Pty(j), Zs(k), i, j, k);
            end
            planeXZ(i, k) = 5e3;
        end
    end
    waitbar((steps+i)/max)
end

%% Plane YZ
i = ceil(length(Ptx)/2); % x = 0
fprintf("i: %f  Ptx(i): %f\n", i, Ptx(i));
planeYZ = zeros(length(Pty), length(Zs)); % contains the condition numbers of the chosen YZ plane 
for j = 1:length(Pty) %y
    for k = 1:length(Zs) %z
        % linearizing
        x0(1) = Ptx(i); x0(2) = Pty(j); x0(3) = Zs(k); % linearizing
        for l = 1:dimU % linearizing around this point
            B(:,l) = (sys.f(x0,uLp+(l==1:dimU)'*delta) ...
                     -sys.f(x0,uLp-(l==1:dimU)'*delta)) ...
                     /(2*delta);
        end
        temp = B(7:11, :); % coping only the non 0 rows
        Cns(i, j, k) = cond(pinv(temp)); % saving the condition number of the Moore-Penrose pseudinverse  
        % saving the CN
        planeYZ(j, k) = Cns(i, j, k);
        if planeYZ(j, k) > 5e3 % for some reason max() gives an error            
            if Cns(i, j, k) > 1e14 && Pty(j) ~= 0
                    fprintf("possible equilibrium at x: %f y: %f z: %f, or i: %f j: %f k: %f\n", Ptx(i), Pty(j), Zs(k), i, j, k);
            end
            planeYZ(j, k) = 5e3;
        end
    end
    waitbar((2*steps+i)/max)
end

close(h);

%% Plotting XY
figure('Name', 'XY Plane'); 
xlabel('X'); ylabel('Y'); hold on;
grid minor; axis equal;
[X, Y] = meshgrid(Ptx, Pty);
contour3(X, Y, planeXY, 'ShowText','on');
c.LineWidth = 3;
hold off;

figure('Name', 'XY Plane v2');
grid minor; axis equal;
[X, Y] = meshgrid(Ptx, Pty);
surf(Ptx, Pty, planeXY);
xlabel('X'); ylabel('Y'); zlabel('Cn')

%% Plotting XZ
figure('Name', 'XZ Plane');
xlabel('X'); ylabel('Z'); hold on;
grid minor; axis equal;
[X, Z] = meshgrid(Ptx, Zs);
contour3(X, Z, planeXZ.', 'ShowText','on');
c.LineWidth = 6;
hold off;

figure('Name', 'XZ Plane v2'); 
grid minor; axis equal;
[X, Z] = meshgrid(Ptx, Zs);
surf(Ptx, Zs, planeXZ.');
xlabel('X'); ylabel('Z'); zlabel('Cn')

%% Plotting YZ
figure('Name', 'YZ Plane');
xlabel('Y'); ylabel('Z'); hold on;
grid minor; axis equal;
[Y, Z] = meshgrid(Pty, Zs);
contour3(Y, Z, planeYZ.', 'ShowText','on');
c.LineWidth = 6;
hold off;

figure('Name', 'YZ Plane v2');
grid minor; axis equal;
[Y, Z] = meshgrid(Pty, Zs);
surf(Pty, Zs, planeYZ.');
xlabel('Y'); ylabel('Z'); zlabel('Cn')
