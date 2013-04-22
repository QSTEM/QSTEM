% File generate rotation angle:
% function [rotAngles,Mrot] = rotateZoneAxis(zone,M,refZone)
function [rotAngles,Mrot] = rotateZoneAxis(zone,M,refZone)

if nargin < 3
    refZone = [0; 0; 1];
end
if nargin < 2
    M = eye(3);
end
if isempty(M)
    M = eye(3);    
end
    
    
% The rotations are in the following sequence: X-Y-Z
% The rotation about the z-axis will be last.
% The rotation angles must tilt the vector 'zone' to refZone
[Ny,Nx] = size(zone);
if (Ny == 1), zone = zone.'; end
% zone = M*zone;
zone2 = zone/sqrt(sum(zone.^2));

[Ny,Nx] = size(refZone);
if (Ny == 1), refZone = refZone.'; end
% refZone = M*refZone;
refZone2 = refZone/sqrt(sum(refZone.^2));



rotAngles = 20*rand(1,3);

[rotAngles,chi2] = fminsearch(@(rotAngles) rotationCostFunc(rotAngles,M,zone2,refZone2), rotAngles);
% [rotAngles,chi2] = ga(@(rotAngles) rotationCostFunc(rotAngles,zone),3);
% rotAngles(3) = 0;
[chi2,Mrot] = rotationCostFunc(rotAngles,M,zone2,refZone2);
fprintf('chi2 = %f\n',chi2);
% Mrot*M

function [chi2,Mrot] = rotationCostFunc(rotAngles,M,zone,refZone)
% Generate the rotation matrix with which one can rotate any atom position:
Mrot = zeros(3);
phi_x = rotAngles(1)*pi/180;
phi_y = rotAngles(2)*pi/180;
phi_z = rotAngles(3)*pi/180;

if (1)
    cx = cos(phi_x);    sx = sin(phi_x);
    cy = cos(phi_y);    sy = sin(phi_y);
    cz = cos(phi_z);    sz = sin(phi_z);

    Mx = [1 0 0;0 cx -sx;0 sx cx];
    My = [cy 0 sy;0 1 0; -sy 0 cy];
    Mz = [cz -sz 0;sz cz 0; 0 0 1];
    Mrot = Mz*My*Mx;
else
    Mrot(1,1) = cos(phi_z)*cos(phi_y);
    Mrot(1,2) = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
    Mrot(1,3) = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);
    
    Mrot(2,1) = sin(phi_z)*cos(phi_y);
    Mrot(2,2) = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
    Mrot(2,3) = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);
    
    Mrot(3,1) = -sin(phi_y);
    Mrot(3,2) = cos(phi_y)*sin(phi_x);
    Mrot(3,3) = cos(phi_y)*cos(phi_x);
end

v = Mrot*M*zone;
v = v/sqrt(sum(v.^2));
% [v refZone]
chi2 = sum((v-refZone).^2);
% check also the direction of the x-vector:
% It should point to the right.
if (refZone(1) == 0)
    vx = Mrot*M*[1;0;0];
    chi2 = chi2+vx(2)^2;
end
