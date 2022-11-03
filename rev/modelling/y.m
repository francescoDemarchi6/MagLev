%{
    simulates y output of the system when u=0 and we move the levmag
    checks if the sensors can estimate the position of the magnet
%}

clear; close all;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat');

approximationType = input("approxType [0/1]> ");
eccentricity = input("Eccentricity [0:1]> ");
eccentricity = max (eccentricity, 0); eccentricity = min (eccentricity, 1); 

if(approximationType == 0)
    eq = results.zeq.zeq_fst;
    params.magnets.I = results.neo_vs_neo.curr_fst;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_fst;
else
    eq = results.zeq.zeq_acc;
    params.magnets.I = results.neo_vs_neo.curr_acc;
    params.levitatingmagnet.I = results.neo_vs_lev.curr_acc;
end

%% EDIT PARAMS
steps = 300; %divisible by 6
L = .01; %double it to obtain the cube side surrounding the path

%{
    create a W oblique movement trace

      o       ^
      |       |
 2L   |   ^   |   x-y plane [meaning of "W"]  
      | /   \ |
       /     \

          2L

    on the x-plane it goes (for every of the 4 segment): down-up-down-up
    [meaning of "oblique"]

%}
px = [-L*ones(1,steps/3) linspace(-L,0,steps/6) linspace(0,L,steps/6) L*ones(1,steps/3)];
py = [linspace(-L,L,steps/3) linspace(L,0,steps/6) linspace(0,L,steps/6) linspace(L,-L,steps/3)];
pz = [linspace(L,-L,steps/3) linspace(-L,0,steps/6) linspace(0,L,steps/6) linspace(L,-L,steps/3)];
pz = pz+eq;

figure(1);
plot3(px,py,pz); hold on;
plot3(px(1), py(1), pz(1), 'o');
plot3(px(end), py(end), pz(end), '^');
xlabel('x'); ylabel('y'); zlabel('z');

x0 = zeros(12,1); x0(3) = eq;
sys = maglevSystem(x0, params, approximationType, eccentricity);
draw(sys,'fancy'); hold off;

yeq = sys.h(x0,zeros(params.solenoids.N,1));
bx_eq = yeq(1);
by_eq = yeq(5);
bz_eq = yeq(9);

% empty arrays
bx = zeros(1,steps); by = bx; bz = bx;

for i=1:steps
    x0(1) = px(i); x0(2) = py(i); x0(3) = pz(i);
    temp = sys.h(x0,zeros(params.solenoids.N,1));
    %sys.h is a 9 values-variable that stores Bx,By,Bz for every of the 3
    %sensors. In the real-life x-sensor can estimate only Bx, y only By and
    %z only Bz. We proceed to trash the other values
    bx(i) = temp(1);
    by(i) = temp(5);
    bz(i) = temp(9);
end

figure(2);
plot3(bx,by,bz); hold on;
plot3(bx(1), by(1), bz(1), 'o');
plot3(bx(end), by(end), bz(end), '^');
xlabel('x'); ylabel('y'); zlabel('z');

dbx = bx-bx_eq;
dby = by-by_eq;
dbz = bz-bz_eq;

figure(3);
plot3(dbx,dby,dbz); hold on;
plot3(dbx(1), dby(1), dbz(1), 'o');
plot3(dbx(end), dby(end), dbz(end), '^');
xlabel('x'); ylabel('y'); zlabel('z');