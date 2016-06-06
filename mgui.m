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

function mgui (pos, rot, rad, camera, angle)

n = max(size(rad));
m = size(pos,1);

rmax = max(rad);
xmin = min(min(pos(:,:,1))) - rmax;
xmax = max(max(pos(:,:,1))) + rmax;
ymin = min(min(pos(:,:,2))) - rmax;
ymax = max(max(pos(:,:,2))) + rmax;
zmin = min(min(pos(:,:,3))) - rmax;
zmax = max(max(pos(:,:,3))) + rmax;

view(3);
%title('DEM simulation')
%xlabel('x'); ylabel('y'); zlabel('z') 

set (gcf, 'Renderer', 'opengl');

[x, y, z] = sphere;

% Draw spheres
for j = 1:m
  clf;
  if exist ('angle', 'var')
    angle_mode = 'manual';
    angle_value = angle;
  else
    angle_mode = 'auto';
    angle_value = 1;
  end
  ax = axes('XLim',[xmin xmax],'YLim',[ymin ymax],'ZLim',[zmin zmax],...
    'Visible', 'on', 'XTick', [], 'YTick', [], 'ZTick', [], ...
    'CameraPosition', camera, 'CameraViewAngle', angle_value, ...
    'CameraViewAngleMode', angle_mode);
  daspect([1,1,1]);
  h1 = surface (x, y, z, 'EdgeColor', 'none');
  
  for i = 1:n
    t2 = hgtransform ('Parent',ax);
    h2 = copyobj (h1,t2);
    T = makehgtform('translate', reshape(pos(j,i,:),3,1));
    S = makehgtform('scale', rad(i));
    R = reshape(rot(j,3*(i-1)+1:3*i,:),3,3);
    T(1:3,1:3) = R;
    set (t2, 'Matrix', T*S);
  end
  
  set(h1,'visible','off');  
  drawnow;
end

end
