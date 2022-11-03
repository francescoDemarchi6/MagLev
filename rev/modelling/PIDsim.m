%{
    Simulates the system controlled by a PID controller equivalent to the one
    implemented on the Teensy microcontroller.
    Requires PIDfun.m to function
%}
clear; close all;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat');

approximationType = input("approxType [0/1]> ");

if(approximationType == 0)
    eq = results.zeq.zeq_fst;
    params.magnets.I = results.neo_vs_neo.curr_fst;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_fst;
else
    eq = results.zeq.zeq_acc;
    params.magnets.I = results.neo_vs_neo.curr_acc;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_acc;
end

%% Declaring PID parameters
PIDparams.KPx = 1.5; %2.6
PIDparams.KIx = 1.5;   %1
PIDparams.KDx = 0;

PIDparams.KPy = 1.5; %2.6
PIDparams.KIy = 1.5;   %1
PIDparams.KDy = 0;

PIDparams.KPz = 1.5; %-1.2
PIDparams.KIz = 1.5; %-0.8
PIDparams.KDz = 0;

x0 = zeros(12,1); x0(3) = eq;
sys = maglevSystem(x0, params, approximationType);
r = sys.h(x0, zeros(4,1)); % Getting reference sensor reading
PIDparams.rx = r(1); PIDparams.ry = r(5); PIDparams.rz = r(9);

%% Simulating the system
tspan = linspace(0,0.1,10);
[t,x] = ode45(@(t,x) odefunc(t,x,sys, PIDparams),tspan,x0);

%% Plotting the results
figure(2)
hold on;
plot(t,x(:,1))
plot(t,x(:,2))
plot(t,x(:,3))
legend('x','y','z');
hold off;
grid;
ylim([0 0.1])
xlim([0 2])

figure(3)
hold on;
plot(t,x(:,4))
plot(t,x(:,5))
plot(t,x(:,6))
legend('psi','theta','phi');
hold off;
grid;
ylim([-0.5 0.5])
xlim([0 2])
%%

function dxdt = odefunc(t, x, sys, PIDparams)
    persistent y
    if t == 0
        u = zeros(4,1);
    else
        u = PIDfun(y,PIDparams);
    end
    y = round(sys.h(x,u),4);
    temp = sys.f(x,u);
    dxdt = temp;  
end
