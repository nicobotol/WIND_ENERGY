function [cl, cd] = interpolate_cl_cd(alpha)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%interpolate the values to the different thicknesses
for k=1:6; % k indicate the airfoil
clthick(k) = aoa(:,k),cl(:, angle_of_attack attack);
cdthick(k) = aoa(:,k),cd(:, angle_of_attack attack);
end

% then interpolate to the actual thickness
% thick_prof =(100,60,48,36,30.1,24.1), i indicates the element nr .
clift = thick_prof clthick (:),thick(
cdrag = thick_prof cdthick (:),thick(

end