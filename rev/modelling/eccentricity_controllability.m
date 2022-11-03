%{
    measures the volume of the controllable area*
    for different eccentricities of the system (both magnets and solenoids)

    use the paramiters / figures for a single OR varying eccentricity, do
    not use both

    *the exact controllable area is difficult to measure, so the 
    smallest circumbscribing rectangoloid is calculated instead
%}

clear; close all;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat');

approximationType = input("approxType [0/1]> ");

%% Searching parameters
steps = 512;
L = .10;
% parameters for a varying eccentricity
min_eccentricity = 0; max_eccentricity = .5;
eccentricity_steps = max_eccentricity * 200; % plotting every .005
% parameters for a single eccentricity to check if the program
% works (faster to compare with solenoids_controlling_power
% and z_graph)0
%min_eccentricity = .25; max_eccentricity = min_eccentricity;
%eccentricity_steps = 1;

% load correct parameters
if(approximationType == 0)
    eq = results.zeq.zeq_fst;
    params.magnets.I = results.neo_vs_neo.curr_fst;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_fst;
else
    eq = results.zeq.zeq_acc;
    params.magnets.I = results.neo_vs_neo.curr_acc;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_acc;
end

%% Derived parameters
% Doubles the space for two halfs respectively for positive error
% (solenoids = +-) and for negative (solenoids = -+)
Xs = linspace(-L/2, L/2, 2*steps);
Ys = linspace(-L, L, 2*steps);
Zs = linspace(.04, .07, steps);

Es = linspace(min_eccentricity, max_eccentricity, eccentricity_steps); % eccentricities evaluated
Volumes = zeros(eccentricity_steps, 7); % will contain the volumes and edges of the rectangoloids
eqs = zeros(length(Es), 1);

Fxs = zeros(1,length(Xs));
Fys = zeros(1,length(Ys));
Fzs = zeros(1,length(Zs));

h = waitbar(0);

%% Start searching
x0 = zeros(12,1);
sys = maglevSystem(x0, params, approximationType, 0);

for ec = 1:length(Es)    
    Zs = linspace(.04, .07, 2*steps); % these parameters need to be reset to find eq_z 
    Fzs = zeros(1,length(Zs));
    % searching for the new equilibrium to reset Zs
    for k = 1:length(Zs) 
            temp = sys.f([0,0,Zs(k),zeros(1,9)]',zeros(params.solenoids.N,1));
            Fzs(k) = temp(9);        
    end
    [~,idx]=min(abs(Fzs)); % new equilibrium found
    eqs(ec) = Zs(idx);
    Zs = linspace(eqs(ec)-L/2, eqs(ec)+L/2, 2*steps); % resetting Zs around the new equilibria
    
    % calculating points along X axis where Fx is closest to 0
    x0 = zeros(12,1); x0(3) = eqs(ec);
    sys = maglevSystem(x0, params, approximationType, Es(ec));

    for i = 1:length(Xs)/2
        x0(1) = Xs(i);
        temp = sys.f(x0,[.5 0 -.5 0]');
        Fxs(i) = temp(7);
    end
    for i = length(Xs)/2+1:length(Xs)
        x0(1) = Xs(i);
        temp = sys.f(x0,[.5 0 -.5 0]');
        Fxs(i) = temp(7);
    end

    % calculating points along Y axis where Fy is closest to 0
    x0 = zeros(12,1); x0(3) = eqs(ec);
    sys = maglevSystem(x0, params, approximationType, Es(ec));

    for j = 1:length(Ys)/2
        x0(2) = Ys(j);
        temp = sys.f(x0,[0 .5 0 -.5]');
        Fys(j) = temp(8);
    end
    for j = length(Ys)/2+1:length(Ys)
        x0(2) = Ys(j);
        temp = sys.f(x0,[0 .5 0 -.5]');
        Fys(j) = temp(8);
    end

    % calculating points along Z axis where Fz is closest to 0
    x0 = zeros(12,1); x0(3) = eqs(ec);
    sys = maglevSystem(x0, params, approximationType, Es(ec));

    for k = 1:length(Zs)/2
        x0(3) = Zs(k);
        temp = sys.f(x0,[.5 .5 .5 .5]');
        Fzs(k) = temp(9);
    end
    for k = length(Zs)/2+1:length(Zs)
        x0(3) = Zs(k);
        temp = sys.f(x0,[-.5 -.5 -.5 -.5]');
        Fzs(k) = temp(9);
    end

    % finding the indexes of the coordinates of the points where the components of F are smallest
    idxs = zeros(1,6);
    [~,idxs(1)] = min(abs(Fxs(1:end/2))); % smallest force for x < 0
    [~,idxs(2)] = min(abs(Fxs(end/2+1:end))); idxs(2) = idxs(2) + length(Fxs)/2; % smallest force for x > 0
    [~,idxs(3)] = min(abs(Fys(1:end/2))); % smallest force for y < 0
    [~,idxs(4)] = min(abs(Fys(end/2+1:end))); idxs(4) = idxs(4) + length(Fys)/2; % smallest force for y > 0
    [~,idxs(5)] = min(abs(Fzs(1:end*3/8))); % smallest force for z << eq
    [~,idxs(6)] = min(abs(Fzs(end/2+1:end))); idxs(6) = idxs(6) + length(Zs)/2; % smallest force for z > eq

    % dimensions of the rectangoloid
    [Volumes(ec, 1), Volumes(ec, 2)] = deal(Xs(idxs(1)),Xs(idxs(2))-Xs(idxs(1))); % vertex x, edge xl
    [Volumes(ec, 3), Volumes(ec, 4)] = deal(Ys(idxs(3)),Ys(idxs(4))-Ys(idxs(3))); % vertex y, edge yl
    [Volumes(ec, 5), Volumes(ec, 6)] = deal(Zs(idxs(5)),Zs(idxs(6))-Zs(idxs(5))); % vertex z, edge zl
    Volumes(ec, 7) = Volumes(ec, 2) * Volumes(ec, 4) * Volumes(ec, 6) * 1e6; % saving the volume in cm^3
    waitbar(ec/length(Es))
end

[~, idx] = max(Volumes(:,7));
fprintf("The eccentricity that grants the biggest controllable volume is: %f\n", Es(idx));
fprintf("Volume is: %f\n", Volumes(idx, 7));

close(h);

%% Plotting

% figures for a varying eccentricity
figure('Name', 'EQz / ecc');
plot(Es, eqs);
xlabel('Eccentricity'); ylabel('Height of the equilibrium')

figure('Name', 'Controllable Volume');
plot(Es, Volumes(:,7));
xlabel('Eccentricity'); ylabel('Volume')

%{
% figures for a single eccentricity
figure(1);
plot3(Xs,zeros(1,length(Xs)),Fxs); hold on;
plot3(zeros(1,length(Ys)),Ys,Fys);
xlabel('X'); ylabel('Y'); zlabel('Fx/Fy')
grid minor; xline(0); yline(0);
% rescaling X-Y axes to be equal, leaving Z axis as unvaried
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
plot3(Xs(idxs(1, 1:2)),0,0,'ko');
plot3(0,Ys(idxs(1, 3:4)),0,'ko'); hold off;

figure(2);
plot(Zs,Fzs); hold on;
xlabel('Z'); ylabel('Fz')
grid minor; yline(0);
Z = Zs(idxs(1,5:6));
plot(Z,0,'ko'); hold off;
%}
% plots the biggest (or only) rectangoloid
figure('Name', 'Controllability Region');
x0 = zeros(12,1); x0(3) = eq;
sys = maglevSystem(x0, params, approximationType, Es(idx)); 
hold on; draw(sys, 'fancy'); grid on; axis equal; view([45,30]);
plotcube([Volumes(idx, 2) Volumes(idx, 4) Volumes(idx, 6)], ...
         [Volumes(idx, 1) Volumes(idx, 3) Volumes(idx, 5)], .1, [1,0,0]); hold off;
     
     
     