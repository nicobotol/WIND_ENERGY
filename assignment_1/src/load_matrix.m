function [cl_mat, cd_mat] = load_matrix()
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

filenames = ["airfoil_data\FFA-W3-241", "airfoil_data\FFA-W3-301", "airfoil_data\FFA-W3-360", "airfoil_data\FFA-W3-480", "airfoil_data\FFA-W3-600", "airfoil_data\cylinder"];

for i = 1:6
  mat = readmatrix(filenames(i));
  cl_mat(:,i) = mat(:,2);
  cd_mat(:,i) = mat(:,3);
end
end