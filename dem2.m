
% many particles

n = 10;

density = 1E3;

r = 0.05;

for i = 1:n
  center(i,:) = [0.0, 0, 2*r*i-r];
  radius(i) = r;
  linear(i,:) = [0, 0, 0];   
  angular(i,:) = [0, 0, 0];
end

angular(2,:) = [0, 20, 0];
    
gravity = [0, 0, -10];

plane = [0, 0, 0, 0, 0, 1;
        -1, 0, 0, 1, 0, 0;
         1, 0, 0,-1, 0, 0;
         0,-1, 0, 0, 1, 0;
         0, 1, 0, 0,-1, 0];

spring = 1E4;

damper = 1.0;

friction = [0.5, 0.5, 0, 0.0];

dem = mdem();

[rotation, inertia, inverse, ...
   mass, invm, plane, contact] = dem.init (center, radius, density, plane);

step = 0.1*dem.crit (mass, spring, damper);

m = 1500;
pos = zeros(m,n,3);
rot = zeros(m,n*3,3);
progressbar;

for i = 1:m

  [center, rotation, linear, ...
          angular, contact] = dem.go1 (step, center, rotation, ...
          linear, angular, radius,inertia, inverse, mass, invm, ...
          gravity, spring, damper, friction, plane, contact);
  
  pos(i,:,:) = center;
  rot(i,:,:) = rotation;
  progressbar (i/m);
end

mgui (pos, rot, radius, [5,5,5]);