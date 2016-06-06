% Copyright (c) 2016, Tomasz Koziara
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function out = mdem()
  out.init = @init;
  out.crit = @crit;
  out.go1 = @go1;
  out.go2 = @go2;
end

% initialise contact and rotation data; compute inertia properties
function [rotation, inertia, inverse, ...
        mass, invm, pout, contact] = init (center, radius, density, plane)
  n = size (center, 1);
  m = size (plane, 1);
  contact = {[],cell(n,n),[],cell(n,m)};
  rotation = zeros (n*3,3);
  mass = density*(4./3.)*pi*radius.*radius.*radius;
  invm = 1./mass;
  inertia = zeros (n*3,3);
  inverse = zeros (n*3,3);
  % set up inertia
  for i = 1:n
    rotation (3*(i-1)+1:3*i,:) = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    x = 0.4*mass(i)*radius(i)*radius(i);
    y = 1.0/x;
    inertia(3*(i-1)+1:3*i,:) = [x, 0, 0; 0, x, 0; 0, 0, x];
    inverse(3*(i-1)+1:3*i,:) = [y, 0, 0; 0, y, 0; 0, 0, y];
  end
  % normalise planes
  for i = 1:m
    plane (i,4:6) = plane(i,4:6)/norm(plane(i,4:6));
  end
  pout = plane;
end

function [W] = skew (O)
  W = [0, -O(3), O(2);
       O(3), 0, -O(1);
      -O(2), O(1), 0];
end
               
function R = expmap (O)
  dp = dot(O,O);
  ln = sqrt(dp);
  W = skew (O);
  if ln < 1E-15
    R = eye(3,3);
  else
    R = eye(3,3) + (sin(ln)/ln)*W + ((1-cos(ln))/dp)*W*W;
  end
end

% initial time integration half-step
function [center, rotation, force] = half1 (step, center, rotation, ...
                         linear, angular, inertia, inverse, mass, gravity)
  n = size(center,1);
  force = zeros (n,6);
  for i = 1:n
      idx = 3*(i-1)+1:3*i;
      O = angular(i,:)';
      J = inertia(idx,:);
      I = inverse(idx,:);
      R1 = expmap (0.5*step*O);
      R = R1*rotation(idx,:);
      T = [0, 0, 0]'; % could be applied torques = R'* t (spatial)
      O = I*(R1'*J*O + 0.5*step*T);
      force(i,1:3) = T - skew(O)*J*O;
      force(i,4:6) = mass(i)*gravity;
      rotation (idx,:) = R;
  end
  center = center + 0.5*step*linear;
end

% final time integration half-step
function [center, rotation, linear, angular] = half2 (step, center, ...
                          rotation, linear, angular, inverse, invm, force)
  n = size(center,1);
  for i = 1:n
      idx = 3*(i-1)+1:3*i;
      I = inverse(idx,:);
      angular (i,:) = angular(i,:) + (step*I*force(i,1:3)')';
      linear (i,:) = linear(i,:) + step*invm(i)*force(i,4:6);
      R1 = expmap (0.5*step*angular(i,:));
      rotation (idx,:) = R1*rotation(idx,:);
  end
  center = center + 0.5*step*linear;
end

% contact detection update
function [contact] = condet (center, radius, plane, contact)
  kd = KDTreeSearcher (center);
  n = size (center,1);
  m = size (plane,1);
  cpnt1 = contact{1};
  data1 = contact{2};
  % release previous sphere to sphere contacts
  for k = 1:size(cpnt1,1)
      i = cpnt1(k,1);
      j = cpnt1(k,2);
      d = norm(center(i,:)-center(j,:));
      depth = radius(i)+radius(j)-d;
      if depth < 0
          data1{i,j} = []; %invalidate history
      end
  end
  cpnt2 = contact{3};
  data2 = contact{4};
  % release previous sphere ro plane contacts
  for k = 1:size(cpnt2,1)
      i = cpnt2(k,1);
      j = cpnt2(k,2);
      depth = radius(i) - (center(i,:)-plane(j,1:3))*plane(j,4:6)';
      if depth < 0
          data2{i,j} = []; %invalidate history
      end
  end
  % detected contacts between sphere pairs
  rmax = max(radius);
  cpnt1 = [];
  for i = 1:n
      [idx, dst] = rangesearch (kd, center(i,:), radius(i)+rmax);
      j = idx{1};
      d = dst{1};
      for k = 1 : size(j,2)
          if i < j(k)
              depth = radius(i)+radius(j(k))-d(k);
              if depth > 0
                normal = center(i,:)-center(j(k),:);
                normal = normal / norm(normal);
                point = 0.5*(center(i,:)+center(j(k),:));
                item = data1{i,j(k)}; % upper diagonal
                if size(item) == [0,0]
                  data1{i,j(k)} = [0,0,0,0,0];
                end
                cpnt1 = [cpnt1; [i, j(k), point, normal, depth, [0,0,0]]];
              end
          end
      end
  end
  % detect contacts between spheres and planes
  cpnt2= [];
  for j = 1:m
      a = plane(j,1:3);
      normal = plane(j,4:6);
      for i = 1:n
          x = center(i,:);
          r = radius(i);
          y = x-a;
          depth = r - y*normal';
          if depth > 0
              point = x - r*normal;
              item = data2{i,j};
              if size(item) == [0,0]
                data2{i,j} = [0,0,0,0,0];
              end
              cpnt2 = [cpnt2; [i, j, point, normal, depth, [0,0,0]]];
          end
      end
  end
  contact = {cpnt1, data1, cpnt2, data2};
end

function [base] = localbase (n)
  e = [[1., 0., 0.]; [0., 1., 0.]];
  a = e(1,2)*n(3)-e(1,3)*n(2);
  b = e(1,3)*n(1)-e(1,1)*n(3);
  c = e(1,1)*n(2)-e(1,2)*n(1);
  len = sqrt(a*a+b*b+c*c);
  if (len < 1E-6)
    a = e(2,2)*n(3)-e(2,3)*n(2);
    b = e(2,3)*n(1)-e(2,1)*n(3);
    c = e(2,1)*n(2)-e(2,2)*n(1);
    len = sqrt(a*a+b*b+c*c);
  end
  a = a/len;
  b = b/len;
  c = c/len;
  d = b*n(3)-c*n(2);
  e = c*n(1)-a*n(3);
  f = a*n(2)-b*n(1);
  base = [[a,b,c]',[d,e,f]',n'];
end

% incremental elastop-plastic frictional force
function [RT,S1,S2] = fric1 (step, S1, S2, spring, eta, friction, RN, U)
  springT = (2./7.)*spring; % J. Shafer, J Phys. France, 6 (1996) 5-20
  etaT = 0.5*eta; % FIXME (what value is most commonly used?)
  RT = zeros(2,1);
  S1 = S1 + step*U(1);
  S2 = S2 + step*U(2);
  RT(1) = (springT*S1 + etaT*U(1)); % FIXME (why this sign works?)
  RT(2) = (springT*S2 + etaT*U(2));
  len = norm(RT);
  lim = friction(1)*RN;
  if len > lim
    lim = friction(2)*RN/len;
    RT = lim * RT;
  end
end

% incremental elasto-plastic rolling and rilling torque
function [t,T1,T2,T3] = rolldril1 (step, T1, T2, T3, E, Ii, Ij, mi, mj, ...
                     Ri, Rj, Oi, Oj, ri, rj, spring, damper, friction, RN)
  if rj > 0.0 % sphere to sphere contact
    oi = Ri*Oi;
    oj = Rj*Oj;
    O = E'*(oj - oi); % rolling (tangential) and drilling (normal) relative angular velocity
    r = ri*rj/(ri + rj);
    I = 1.0/(1.0/(trace(Ii)/3+mi*ri*ri) + 1.0/(trace(Ij)/3+mj*rj*rj)); % 20 in [1]
  else % sphere to plane contact
    oi = Ri*Oi;
    O = -E'*oi;
    r = ri;
    I = trace(Ii)/3+mi*ri*ri;
  end

  mur = friction(3);
  mud = friction(4);
  %k = 0.5*spring*r*r; % (7a) in [1] J. Ai et al. / Powder Technology 206 (2011) 269-282
  %k = 3*spring*mur*mur*r*r; % (7b) in [1]
  k = r*RN; % (7c) in [1]
  e = 2*sqrt(I*k); % (18) in [1] with eta_t = 1.0 (critically damped)
  T1 = T1 + step*k*O(1);
  T2 = T2 + step*k*O(2);
  T3 = T3 + step*k*O(3);    
  lim = mur*r*RN;
  len = (T1*T1+T2*T2);
  T = zeros(3,1);
  if len > lim
      eps = lim/len;
      T1 = eps*T1;
      T2 = eps*T2;
      T(1) = T1;
      T(2) = T2;
  else
      T(1) = T1 + e*O(1);
      T(2) = T2 + e*O(2);
  end
  lim = mud*r*RN;
  len = abs(T3);
  if len > lim
      eps = lim/len;
      T3 = eps*T3;
      T(3) = T3;
  else
      T(3) = T3 + e*O(3);
  end
  t = E*T;
end

% classical DEM contact resolution
function [force, contact] = resolve1 (step, center, rotation, linear, ...
  angular, inertia, mass, radius, spring, damper, friction, contact, force)
  cpnt1 = contact{1};
  data1 = contact{2};
  cpnt2 = contact{3};
  data2 = contact{4};
  n = size(cpnt1,1);
  m = size(cpnt2,1);
  % resolve sphere to sphere contacts
  for k = 1:n
      i = cpnt1(k,1);
      j = cpnt1(k,2);
      point = cpnt1(k,3:5);
      normal = cpnt1(k,6:8);
      depth = cpnt1(k,9);
      hist = data1{i,j};
      S1 = hist(1);
      S2 = hist(2);
      T1 = hist(3);
      T2 = hist(4);
      T3 = hist(5);      
      E = localbase (normal);
      Ri = rotation(3*(i-1)+1:3*i,:);
      Rj = rotation(3*(j-1)+1:3*j,:);
      Ii = inertia(3*(i-1)+1:3*i,:);
      Ij = inertia(3*(j-1)+1:3*j,:);
      mi = mass(i);
      mj = mass(j);
      ri = radius(i);
      rj = radius(j);
      Oi = angular(i,:)';
      Oj = angular(j,:)';
      Hi = E'*[skew(center(i,:)-point)*Ri,eye(3,3)];
      Hj = E'*[skew(center(j,:)-point)*Rj,eye(3,3)];
      Ui = Hi*[angular(i,:),linear(i,:)]';
      Uj = Hj*[angular(j,:),linear(j,:)]';
      U = Uj-Ui;
      mij = mi*mj/(mi+mj);
      eta = damper*2*sqrt(spring*mij);
      RN = depth*spring + eta*U(3);
      [RT,S1,S2] = fric1 (step, S1, S2, spring, eta, friction, RN, U);
      [t,T1,T2,T3] = rolldril1 (step, T1, T2, T3, E, Ii, Ii, mi, mi, ...
                     Ri, Ri, Oi, Oi, ri, 0, spring, damper, friction, RN);
      R = [RT; RN];
      data1{i,j} = [S1, S2, T1, T2, T3];
      force(i,:) = force(i,:) + (Hi'*R)' + [(Ri'*t)',0,0,0];
      force(j,:) = force(j,:) - (Hj'*R)' - [(Rj'*t)',0,0,0];
      cpnt1(k,10:12) = R; % output force
  end
  % resolve sphere to plane contacts
  for k = 1:m
      i = cpnt2(k,1);
      j = cpnt2(k,2);
      point = cpnt2(k,3:5);
      normal = cpnt2(k,6:8);
      depth = cpnt2(k,9);
      hist = data2{i,j};
      S1 = hist(1);
      S2 = hist(2);
      T1 = hist(3);
      T2 = hist(4);
      T3 = hist(5);
      E = localbase (normal);
      Ri = rotation(3*(i-1)+1:3*i,:);
      Ii = inertia(3*(i-1)+1:3*i,:);
      mi = mass(i);
      ri = radius(i);
      Oi = angular(i,:)';
      Hi = E'*[skew(center(i,:)-point)*Ri,eye(3,3)];      
      Ui = Hi*[angular(i,:),linear(i,:)]';      
      U = -Ui;
      eta = damper*2*sqrt(spring*mi);
      RN = depth*spring + eta*U(3);
      [RT,S1,S2] = fric1 (step, S1, S2, spring, eta, friction, RN, U);
      [t,T1,T2,T3] = rolldril1 (step, T1, T2, T3, E, Ii, Ii, mi, mi, ...
                     Ri, Ri, Oi, Oi, ri, 0, spring, damper, friction, RN);
      R = [RT; RN];
      data2{i,j} = [S1, S2, T1, T2, T3];
      force(i,:) = force(i,:) + (Hi'*R)' + [(Ri'*t)',0,0,0];
      cpnt2(k,10:12) = R; % output force
  end
  contact = {cpnt1,data1,cpnt2,data2};
end

% alternative frictional force
function [RT] = fric2 (step, W, B, friction, RN)  
  RT = W(1:2,1:2)\(B(1:2)+step*W(1:2,3)*RN)/step; % FIXME (why this sign works?)
  len = norm(RT);
  lim = friction(1)*RN;
  if len > lim
    lim = friction(2)*RN/len;
    RT = lim*RT;
  end
end

% alternative rolling and rilling torque
function [t] = rolldril2 (step, E, Ii, Ij, Ri, Rj, Oi, Oj, ...
                          fi, fj, ri, rj, friction, RN)
  if rj > 0.0 % sphere to sphere contact
    oi = Ri*Oi;
    oj = Rj*Oj;
    r = ri*rj/(ri + rj);
    j = Ri*Ii*Ri' + Rj*Ij*Rj';
    bi = step*Ri*Ii*fi;
    bj = step*Rj*Ij*fj;
    tstick = j\[oj-oi+bj-bi]/step;
  else % sphere to plane contact
    oi = Ri*Oi;
    O = -E'*oi;
    r = ri;
    j = Ri*Ii*Ri';
    bi = step*Ri*Ii*fi;
    tstick = j\[-oi-bi]/step;
  end
  
  T = E'*tstick;
  lim = friction(3)*r*RN;
  len = norm(T(1:2));
  if len > lim
      eps = lim/len;
      T(1) = eps*T(1);
      T(2) = eps*T(2);
  end 
  lim = friction(4)*r*RN;
  len = abs(T(3));
  if len > lim
      eps = lim/len;
      T(3) = eps*T(3);
  end
  t = E*T;
end

% alternative contact resolution
function [force, contact] = resolve2 (step, center, rotation, linear, ...
  angular, inverse, mass, invm, radius, spring, damper, friction, contact, force)
  cpnt1 = contact{1};
  data1 = contact{2};
  cpnt2 = contact{3};
  data2 = contact{4};
  n = size(cpnt1,1);
  m = size(cpnt2,1);
  reac = zeros(size(force));
  % resolve sphere to sphere contacts
  for k = 1:n
      i = cpnt1(k,1);
      j = cpnt1(k,2);
      point = cpnt1(k,3:5);
      normal = cpnt1(k,6:8);
      depth = cpnt1(k,9);
      E = localbase (normal);
      Oi = angular(i,:)';
      Oj = angular(j,:)';
      ri = radius(i);
      rj = radius(i);
      Ri = rotation(3*(i-1)+1:3*i,:);
      Rj = rotation(3*(j-1)+1:3*j,:);
      Ii = inverse(3*(i-1)+1:3*i,:);
      Ij = inverse(3*(j-1)+1:3*j,:);
      Mi = blkdiag(Ii,invm(i)*eye(3,3));
      Mj = blkdiag(Ij,invm(j)*eye(3,3));
      Hi = E'*[skew(center(i,:)-point)*Ri,eye(3,3)];
      Hj = E'*[skew(center(j,:)-point)*Rj,eye(3,3)];
      W = Hi*Mi*Hi' + Hj*Mi*Hj';
      Ui = Hi*[angular(i,:),linear(i,:)]';
      Uj = Hj*[angular(j,:),linear(j,:)]';
      U = Uj-Ui;
      B = U + step*(Hj*[Ij*force(j,1:3)';invm(j)*force(j,4:6)'] - ...
                    Hi*[Ii*force(i,1:3)';invm(i)*force(i,4:6)']);
      eta = damper*2*sqrt(spring/W(3,3));
      %RN = (spring*depth+spring*0.5*step*B(3)+eta*B(3))/(1 + W(3,3)*(spring*0.5*step*step+eta*step));
      RN = spring*depth + eta*U(3);
      RT = fric2 (step, W, B, friction, RN);
      R = [RT; RN];
      fi = force(i,1:3)' + Hi(:,1:3)'*R;
      fj = force(j,1:3)' - Hj(:,1:3)'*R;
      t = rolldril2 (step, E, Ii, Ij, Ri, Rj, Oi, Oj, fi, fj, ri, rj, friction, RN);
      reac(i,:) = reac(i,:) + (Hi'*R)' + [(Ri*t)',0,0,0];
      reac(j,:) = reac(j,:) - (Hj'*R)' - [(Rj*t)',0,0,0];
      cpnt1(k,10:12) = R; % output force
  end
  % resolve sphere to plane contacts
  for k = 1:m
      i = cpnt2(k,1);
      j = cpnt2(k,2);
      point = cpnt2(k,3:5);
      normal = cpnt2(k,6:8);
      depth = cpnt2(k,9);
      E = localbase (normal);
      Oi = angular(i,:)';
      ri = radius(i);
      Ri = rotation(3*(i-1)+1:3*i,:);
      Ii = inverse(3*(i-1)+1:3*i,:);
      Mi = blkdiag(Ii,invm(i)*eye(3,3));
      Hi = E'*[skew(center(i,:)-point)*Ri,eye(3,3)];
      W = Hi*Mi*Hi';
      Ui = Hi*[angular(i,:),linear(i,:)]';
      U = -Ui;
      B = U - step*Hi*[Ii*force(i,1:3)';invm(i)*force(i,4:6)'];
      eta = damper*2*sqrt(spring/W(3,3));
      %RN = (spring*depth+spring*0.5*step*B(3)+eta*B(3))/(1 + W(3,3)*(spring*0.5*step*step+eta*step));
      RN = spring*depth + eta*U(3);
      RT = fric2 (step, W, B, friction, RN);
      R = [RT; RN];
      fi = force(i,1:3)' + Hi(:,1:3)'*R;
      t = rolldril2 (step, E, Ii, Ii, Ri, Ri, Oi, Oi, fi, fi, ri, 0, friction, RN);
      reac(i,:) = reac(i,:) + (Hi'*R)' + [(Ri*t)',0,0,0];
      cpnt2(k,10:12) = R; % output force
  end
  force = force + reac;
  contact = {cpnt1,data1,cpnt2,data2};
end

% critical time step computation
function [step] = crit (mass, spring, damper)
  minm = min (mass);
  omega = sqrt(spring/minm);
  step  = (2.0/omega)*(sqrt(1.0+damper*damper) - damper);
end

% classical DEM integration step
function [center, rotation, linear, angular, contact] = go1 (step, ...
          center, rotation, linear, angular, radius, inertia, inverse, ...
          mass, invm, gravity, spring, damper, friction, plane, contact)

  [center, rotation, force] = half1 (step, center, rotation, linear, ...
                               angular, inertia, inverse, mass, gravity);
                           
  [contact] = condet (center, radius, plane, contact);
  
  [force, contact] = resolve1 (step, center, rotation, linear, angular, ...
          inertia, mass, radius, spring, damper, friction, contact, force);
                  
  [center, rotation, linear, angular] = half2 (step, center, rotation, ...
                                  linear, angular, inverse, invm, force);
end

% alternative integration step
function [center, rotation, linear, angular, contact] = go2 (step, ...
          center, rotation, linear, angular, radius, inertia, inverse, ...
          mass, invm, gravity, spring, damper, friction, plane, contact)

  [center, rotation, force] = half1 (step, center, rotation, linear, ...
                               angular, inertia, inverse, mass, gravity);
                           
  [contact] = condet (center, radius, plane, contact);
  
  [force, contact] = resolve2 (step, center, rotation, linear, angular, ...
     inverse, mass, invm, radius, spring, damper, friction, contact, force);
                  
  [center, rotation, linear, angular] = half2 (step, center, rotation, ...
                                  linear, angular, inverse, invm, force);
end
