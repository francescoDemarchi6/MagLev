%{
    show the force felt by the levmag in a 3d cube of points
%}

clear; close all;
addpath('../maglevFunctions');
load('params.mat');
load('results.mat');

approximationType = input("approxType [0/1]> ");
eccentricity = input("Eccentricity [0:1]> ");
eccentricity = max (eccentricity, 0); eccentricity = min (eccentricity, 1); 

%% Searching parameters
steps = 9; % odd number for planes along the axis
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

%% Derived parameters
Ptx = linspace(-L/2, L/2, steps);
Pty = linspace(-L/2*bsqr, L/2*bsqr, steps); % assuming the ellipse's major axis is y
Zs = linspace(eq-L*9/20, eq+L*1/10, steps);
Fzs = zeros(length(Ptx), length(Pty), length(Zs),3); %3D matrix of a 3-page vals

x0 = zeros(12,1); x0(3) = eq;
u = [0; 0; 0; 0]; % x+, y+, x-, y-
sys = maglevSystem(x0, params, approximationType, eccentricity);

% YZ Plane
i = ceil(length(Ptx)/2); % x = 0
for j = 1:length(Pty) %y
    for k = 1:length(Zs) %z
        x0(1) = Ptx(i); x0(2) = Pty(j); x0(3) = Zs(k);      %   YZ plane
        temp = sys.f(x0, u); 
        % normalize to make a good plot
        vec = temp(7:9);
        Fzs(j,i,k,:) = vec./(norm(vec)); %x,y,z forces       
    end
end

%% XZ Plane
j = ceil(length(Pty)/2); % y = 0
for i = 1:length(Ptx) %x
    for k = 1:length(Zs) %z
        x0(1) = Ptx(i); x0(2) = Pty(j); x0(3) = Zs(k);      %   XZ plane
        temp = sys.f(x0, u); 
        % normalize to make a good plot
        vec = temp(7:9);
        Fzs(j,i,k,:) = vec./(norm(vec)); %x,y,z forces       
    end
end


%% Plotter
figure(1);
[X,Y,Z] = meshgrid(Ptx, Pty, Zs);
[U,V,W] = deal(Fzs(:,:,:,1), Fzs(:,:,:,2), Fzs(:,:,:,3));

idx = Fzs(:,:,:,3) >= 0;
quiver3(X(idx),Y(idx),Z(idx), U(idx),V(idx),W(idx), .6, 'r'); hold on;
quiver3(X(~idx),Y(~idx),Z(~idx), U(~idx),V(~idx),W(~idx), .6, 'b'); hold on;
grid on; axis equal; view([45,15]);
xlabel('x'); ylabel('y'); zlabel('z');

load params;
%params.solenoids.ri = 0;
%params.solenoids.ro = 0; % solenoids are covering the forces, better to not plot them
x0 = zeros(12,1); x0(3) = eq;
sys = maglevSystem(x0, params, approximationType, eccentricity);
draw(sys, 'fancy'); hold off;