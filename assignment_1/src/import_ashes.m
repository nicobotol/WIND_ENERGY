function [pn_ashes, pt_ashes, r_ashes] = import_ashes(sensor_file_name)
% Import data from the ashes code

sensor_file_item = size(sensor_file_name, 2);

pn_ashes = zeros(sensor_file_item, 40);
pt_ashes = zeros(sensor_file_item, 40);
r_ashes = zeros(1, 40);

for i=1:sensor_file_item % loop over different files
  file_name = sensor_file_name(i);  % extract the filename from the vector
  loaded_mat = load(file_name);  % load the file

  pn = loaded_mat.data(:, 40*17:40*17+39); % normal load
  pt = loaded_mat.data(:, 40*19:40*19+39); % tangential load

  pn_ashes(i,:) = mean(pn, 1);
  pt_ashes(i,:) = mean(pt, 1);
end

bladepos = load("ashes_data\bladepos.mat");
r_ashes = bladepos.blade_pos(1,:) + 2.8; % load the radial position

end