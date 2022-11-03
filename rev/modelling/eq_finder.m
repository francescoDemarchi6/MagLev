%{
    calculates the height of the equilibrium in (x=0,y=0)
    for different eccentricities of the ellipse

    an ellipses is described by the equation
    x^2/a^2 + y^2/b^2 = 1
    sqrt(1-a^2/b^2) is the eccentricity of the ellipse
%}

clear; close all;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat')

approximationType = input("approxType [0/1]> ");

%% Searching parameters 
min_height = .03;
max_height = .07;
height_steps = 300;
min_eccentricity = 0;
max_eccentricity = .6;
eccentricity_steps = max_eccentricity * 200; % plotting every .005

if(approximationType == 0)
    eq = results.zeq.zeq_fst;
    params.magnets.I = results.neo_vs_neo.curr_fst;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_fst;
else
    eq = results.zeq.zeq_acc;
    params.magnets.I = results.neo_vs_neo.curr_acc;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_acc;
end

x0 = zeros(12,1); x0(3) = eq;

%% Finding the height of the equilibrium for each eccentricity
Zs = linspace(min_height, max_height, height_steps); % heights at which the force will be calculated
Es = linspace(min_eccentricity, max_eccentricity, eccentricity_steps); % eccentricities evaluated
Bs = zeros(size(Es)); % will contain the length of the Y axis of the ellipse
Fzs = zeros(size(Zs)); % vector that will contain the force at each height
Heights = zeros(size(Es)); % vector with the results

h = waitbar(0);
for i = 1:length(Es)
    Bs(i) = params.magnets.R*(1/(1-Es(i)^2));
    sys = maglevSystem(x0, params, approximationType, Es(i)); % modelling a system of eccentricity Es(i)
    for k = 1:length(Zs)
        temp = sys.f([0,0,Zs(k),zeros(1,9)]',zeros(params.solenoids.N,1));
        Fzs(k) = temp(9);        
    end
    [~,idx]=min(abs(Fzs));
    Fzs = zeros(size(Zs)); % resetting Fzs for precaution
    Heights(i) = Zs(idx); % saving the height of the equilibrium
    waitbar(i/length(Es))
end
close(h);

%% Plotting
figure('Name', 'Height of the equilibrium / Eccentricity');
grid minor; axis equal;
plot(Es, Heights)
xlabel('Eccentricity'); ylabel('Height');

figure('Name', 'Height of the equilibrium / Length of the Y axis');
grid minor; axis equal;
plot(Bs, Heights)
xlabel('Length'); ylabel('Height');