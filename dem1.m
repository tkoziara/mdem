
% one particle

density = 1E3;

center = [0, 0, 0.05];
      
radius = [0.05];
     
linear = [0, 0, 0];
      
angular = [0, 20, 0];
    
gravity = [0, 0, -10];

plane = [0, 0, 0, 0, 0, 1];

spring = 1E3;

damper = 0.1;

friction = [0.5, 0.5, 0.1, 0];

% get interface
dem = mdem();

% initialise some data
[rotation, inertia, inverse, ...
   mass, invm, plane, contact] = dem.init (center, radius, density, plane);

% step multipler
mult = 1.0;

% set up integration step
step = mult*0.1*dem.crit (mass, spring, damper);

% number of time steps
m = 300/mult;

% number of particles
n = size(center,1);

% save initial conditions
rot0 = rotation;
cen0 = center;
ang0 = angular;
lin0 = linear;

% allocate animation space
pos = zeros(m,n,3);
rot = zeros(m,n*3,3);

% allocate plots space
ang1 = zeros(m);
ang2 = zeros(m);
time = zeros(m);

% initialise progress bar
progressbar;

% run classical DEM simulation
for i = 1:m
  [center, rotation, linear, ...
          angular, contact] = dem.go1 (step, center, rotation, ...
          linear, angular, radius,inertia, inverse, mass, invm, ...
          gravity, spring, damper, friction, plane, contact);
  
  pos(i,:,:) = center;
  rot(i,:,:) = rotation;
  %ang1(i) = angular(2);
  ang2(i) = linear(3);
  time(i) = step*i;
  progressbar (i/m);
end

% run alternative DEM simulation
% rotation = rot0;
% center = cen0;
% angular = ang0;
% linear = lin0;
% for i = 1:m
%   [center, rotation, linear, ...
%           angular, contact] = dem.go2 (step, center, rotation, ...
%           linear, angular, radius,inertia, inverse, mass, invm, ...
%           gravity, spring, damper, friction, plane, contact);
%   
%   ang2(i) = angular(2);
%   progressbar (i/m);
% end

% show classical DEM animation
mgui (pos, rot, radius, [3,3,3], 3);

% compare angular velocity histories
plot ([time, time], [ang1, ang2]);
