% This script rotates and translates premade ii categories centered at the 
% origin (0,0) to make 2 ii and 2 rb categories

% Obtain screen size and pixel info, then create an appropriate circlespace
% matrix
% visual_angle_in_degrees = 5; 
% dist_to_screen_cm = 46;
% 
% [screen_size_mm(1), screen_size_mm(2)] = Screen('DisplaySize', 0);
% screen_width_cm = screen_size_mm(1)/10;
% screen_height_cm = screen_size_mm(2)/10;
% [screen_width_pixels,screen_height_pixels] = Screen('WindowSize', 0);
% 
% screen_pixel_density = screen_width_pixels/screen_width_cm;
% stim_width = 2*dist_to_screen_cm*tan((visual_angle_in_degrees/2)*(pi/180));
% stim_pixels = stim_width*screen_pixel_density;
% stim_radius = ceil(stim_pixels/2);

clear all
close all
clc

stim_radius = 512;

load zero_center_positive_ii_stim.dat
stim = zero_center_positive_ii_stim;

figure
plot2dstim(stim, [-50 50 -50 50], 0)

exportfig(gcf, [pwd '/output/stim_figs/stim_origin.eps'], 'Color', 'cmyk');

% create positive sloping ii cats
ii_positive(:,1) = stim(:,1);
ii_positive(:,2) = stim(:,2) + 50;
ii_positive(:,3) = stim(:,3) + 50;

% create negative sloping ii cats
theta = pi/2;
rot = [cos(theta) sin(theta);
       -sin(theta) cos(theta);];
   
ii_negative = [stim(:,2) stim(:,3)]*rot;

ii_negative(:,1) = ii_negative(:,1) + 50;
ii_negative(:,2) = ii_negative(:,2) + 50;

ii_negative = [stim(:,1) ii_negative];

% create vertical rb categories
theta = pi/4;
rot = [cos(theta) sin(theta);
       -sin(theta) cos(theta);];
   
rb_vertical = [stim(:,2) stim(:,3)]*rot;

rb_vertical(:,1) = rb_vertical(:,1) + 50;
rb_vertical(:,2) = rb_vertical(:,2) + 50;

rb_vertical = [stim(:,1) rb_vertical];

% create horizontal rb categories
theta = -pi/4;
rot = [cos(theta) sin(theta);
       -sin(theta) cos(theta);];
   
rb_horizontal = [stim(:,2) stim(:,3)]*rot;

rb_horizontal(:,1) = rb_horizontal(:,1) + 50;
rb_horizontal(:,2) = rb_horizontal(:,2) + 50;

rb_horizontal = [stim(:,1) rb_horizontal];

% Plot the stim in 0 to 100 space
xyaxes = [-10 100 -10 100];

figure
subplot(2,2,1)
plot2dstim(ii_positive, xyaxes, 0)
title('ii positive')
xlabel('spatial frequency')
ylabel('orientation')
subplot(2,2,2)
plot2dstim(ii_negative, xyaxes, 0)
xlabel('spatial frequency')
ylabel('orientation')
title('ii negative')
subplot(2,2,3)
plot2dstim(rb_vertical, xyaxes, 0)
title('rb vertical')
xlabel('spatial frequency');
ylabel('orientation')
subplot(2,2,4)
plot2dstim(rb_horizontal, xyaxes, 0)
title('rb horizontal')
xlabel('spatial frequency')
ylabel('orientation');

exportfig(gcf, [pwd '/output/stim_figs/stim_abstract.eps'], 'Color', 'cmyk');

% write the categories to files for later reference
f1 = fopen([pwd '/output/stim_files/ii_positive.dat'], 'w');
f2 = fopen([pwd '/output/stim_files/ii_negative.dat'], 'w');
f3 = fopen([pwd '/output/stim_files/rb_vertical.dat'], 'w');
f4 = fopen([pwd '/output/stim_files/rb_horizontal.dat'], 'w');

fprintf(f1,'%i %7.4f %1.4f \n', ii_positive');
fprintf(f2,'%i %7.4f %1.4f \n', ii_negative');
fprintf(f3,'%i %7.4f %1.4f \n', rb_vertical');
fprintf(f4,'%i %7.4f %1.4f \n', rb_horizontal');

fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);

% transform these categories into physically meaningful dimensions

final_ii_positive(:,1) = ii_positive(:,1);
final_ii_positive(:,2) = 1 + (ii_positive(:,2)/30); %x - axis (cpd)
final_ii_positive(:,3) = (ii_positive(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

final_ii_negative(:,1) = ii_negative(:,1);
final_ii_negative(:,2) = 1 + (ii_negative(:,2)/30); %x - axis (cpd)
final_ii_negative(:,3) = (ii_negative(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

final_rb_vertical(:,1) = rb_vertical(:,1);
final_rb_vertical(:,2) = 1 + (rb_vertical(:,2)/30); %x - axis (cpd)
final_rb_vertical(:,3) = (rb_vertical(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

final_rb_horizontal(:,1) = rb_horizontal(:,1);
final_rb_horizontal(:,2) = 1 + (rb_horizontal(:,2)/30); %x - axis (cpd)
final_rb_horizontal(:,3) = (rb_horizontal(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

% Plot the stim in the transformed space
xyaxes = [1 4.58 0 pi/2+pi/9];

figure
subplot(2,2,1)
plot2dstim(final_ii_positive, xyaxes, 0)
title('final ii positive')
xlabel('spatial frequency (cpd)')
ylabel('orientation (rad)')
subplot(2,2,2)
plot2dstim(final_ii_negative, xyaxes, 0)
title('final ii negative')
xlabel('spatial frequency (cpd)')
ylabel('orientation (rad)')
subplot(2,2,3)
plot2dstim(final_rb_vertical, xyaxes, 0)
title('final rb vertical')
xlabel('spatial frequency (cpd)')
ylabel('orientation (rad)')
subplot(2,2,4)
plot2dstim(final_rb_horizontal, xyaxes, 0)
title('final rb horizontal')
xlabel('spatial frequency (cpd)')
ylabel('orientation (rad)')

exportfig(gcf, [pwd '/output/stim_figs/stim_physical.eps'], 'Color', 'cmyk');

% write the categories to files for later reference
f1 = fopen([pwd '/output/stim_files/final_ii_positive.dat'], 'w');
f2 = fopen([pwd '/output/stim_files/final_ii_negative.dat'], 'w');
f3 = fopen([pwd '/output/stim_files/final_rb_vertical.dat'], 'w');
f4 = fopen([pwd '/output/stim_files/final_rb_horizontal.dat'], 'w');

fprintf(f1,'%i %7.4f %1.4f \n', final_ii_positive');
fprintf(f2,'%i %7.4f %1.4f \n', final_ii_negative');
fprintf(f3,'%i %7.4f %1.4f \n', final_rb_vertical');
fprintf(f4,'%i %7.4f %1.4f \n', final_rb_horizontal');

fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);


% Now we will 
[circlespace,circlespace_frame,meshX,meshY] = make_circlespace(stim_radius);

circlespace_frame = circlespace_frame/2;
x = meshX;
y = meshY;

clear ii_positive
clear final_ii_positive

ii_positive = [1 36 50; 2 64 50];
final_ii_positive(:,1) = ii_positive(:,1);
final_ii_positive(:,2) = 1 + (ii_positive(:,2)/30); %x - axis (cpd)
final_ii_positive(:,3) = (ii_positive(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

for i = 1:2
    
    % frequency in cpd
    f = final_ii_positive(i,2);
    %orientation units: radians
    angle = final_ii_positive(i,3);
    g=cos(angle);
    h=sin(angle);
    alpha = 2*atan(1/70)*(180/pi);
    % We choose the stimulus diameter to match the size of meshX and meshY
    D = sum(circlespace(((size(circlespace,1)+1)/2),:));
    % Wavenumber: cycle/pixel
    k = 2*pi / (D/(alpha * f));
    % Create image to place on screen
    wave=sin((h*x+g*y)*k);
    disc=circlespace.*wave;
    disc = disc + circlespace_frame;
    
    % get cat label
    cat_label = ii_positive(i,1);
    
    % Transform stim dimensions to 0 to 100 space
    x_coord = round(ii_positive(i,2)); % (spatial frequency)
    y_coord = round(ii_positive(i,3)); % (orientation)
    
    % File name has the following format: cat_x_y_fileindex.bmp
%     FILENAME = [pwd '/output/ii_positive/' three_digit_string(cat_label) '_' three_digit_string(x_coord) '_' three_digit_string(y_coord) '_' three_digit_string(i) '.bmp'];
      FILENAME = [pwd '/output//' three_digit_string(cat_label) '_' three_digit_string(x_coord) '_' three_digit_string(y_coord) '_' three_digit_string(i) '.bmp'];
    FMT = 'bmp';
    IMWRITE(disc,FILENAME,FMT)
    
    
    
%     % frequency in cpd
%     f = final_ii_negative(i,2);
%     %orientation units: radians
%     angle = final_ii_negative(i,3);
%     g=cos(angle);
%     h=sin(angle);
%     alpha = 2*atan(1/70)*(180/pi);
%     % We choose the stimulus diameter to match the size of meshX and meshY
%     D = sum(circlespace(((size(circlespace,1)+1)/2),:));
%     % Wavenumber: cycle/pixel
%     k = 2*pi / (D/(alpha * f));
%     % Create image to place on screen
%     wave=sin((h*x+g*y)*k);
%     disc=circlespace.*wave;
%     disc = disc + circlespace_frame;
%     
%     % get cat label
%     cat_label = ii_negative(i,1);
%     
%     % Transform stim dimensions to 0 to 100 space
%     x_coord = round(ii_negative(i,2)); % (spatial frequency)
%     y_coord = round(ii_negative(i,3)); % (orientation)
%     
%     % File name has the following format: cat_x_y_fileindex.bmp
%     FILENAME = [pwd '/output/ii_negative/' three_digit_string(cat_label) '_' three_digit_string(x_coord) '_' three_digit_string(y_coord) '_' three_digit_string(i) '.bmp'];
%     FMT = 'bmp';
%     IMWRITE(disc,FILENAME,FMT)
%     
%     
%     
%     % frequency in cpd
%     f = final_rb_vertical(i,2);
%     %orientation units: radians
%     angle = final_rb_vertical(i,3);
%     g=cos(angle);
%     h=sin(angle);
%     alpha = 2*atan(1/70)*(180/pi);
%     % We choose the stimulus diameter to match the size of meshX and meshY
%     D = sum(circlespace(((size(circlespace,1)+1)/2),:));
%     % Wavenumber: cycle/pixel
%     k = 2*pi / (D/(alpha * f));
%     % Create image to place on screen
%     wave=sin((h*x+g*y)*k);
%     disc=circlespace.*wave;
%     disc = disc + circlespace_frame;
%     
%     % get cat label
%     cat_label = rb_vertical(i,1);
%     
%     % Transform stim dimensions to 0 to 100 space
%     x_coord = round(rb_vertical(i,2)); % (spatial frequency)
%     y_coord = round(rb_vertical(i,3)); % (orientation)
%     
%     % File name has the following format: cat_x_y_fileindex.bmp
%     FILENAME = [pwd '/output/rb_vertical/' three_digit_string(cat_label) '_' three_digit_string(x_coord) '_' three_digit_string(y_coord) '_' three_digit_string(i) '.bmp'];
%     FMT = 'bmp';
%     IMWRITE(disc,FILENAME,FMT)
%     
%     
%     
%     % frequency in cpd
%     f = final_rb_horizontal(i,2);
%     %orientation units: radians
%     angle = final_rb_horizontal(i,3);
%     g=cos(angle);
%     h=sin(angle);
%     alpha = 2*atan(1/70)*(180/pi);
%     % We choose the stimulus diameter to match the size of meshX and meshY
%     D = sum(circlespace(((size(circlespace,1)+1)/2),:));
%     % Wavenumber: cycle/pixel
%     k = 2*pi / (D/(alpha * f));
%     % Create image to place on screen
%     wave=sin((h*x+g*y)*k);
%     disc=circlespace.*wave;
%     disc = disc + circlespace_frame;
%     
%     % get cat label
%     cat_label = rb_horizontal(i,1);
%     
%     % Transform stim dimensions to 0 to 100 space
%     x_coord = round(rb_horizontal(i,2)); % (spatial frequency)
%     y_coord = round(rb_horizontal(i,3)); % (orientation)
%     
%     % File name has the following format: cat_x_y_fileindex.bmp
%     FILENAME = [pwd '/output/rb_horizontal/' three_digit_string(cat_label) '_' three_digit_string(x_coord) '_' three_digit_string(y_coord) '_' three_digit_string(i) '.bmp'];
%     FMT = 'bmp';
%     IMWRITE(disc,FILENAME,FMT)

end