clear all
close all
clc

load_cmd = ['load ' '/Volumes/Media/lab/stimuli_creation/save_disc/pigeon_hybrid/input/hybrid_stim.dat'];
eval(load_cmd);

load_cmd = ['load ' '/Volumes/Media/lab/stimuli_creation/save_disc/pigeon_hybrid/input/hybrid_control_stimuli.dat'];
eval(load_cmd);

stim = hybrid_stim;

plot2dstim(stim)
% title('Hybrid')
% xlabel('spatial frequency')
% ylabel('orientation')

% exportfig(gcf, [pwd '/output/stim_figs/stim_physical.eps'], 'Color', 'cmyk');

% transform these categories into physically meaningful dimensions
% NOTE: I used the same hybrid category stimuli as the  hybrid
% transfer study. And the file I use is already transformed

% fancy figure making now

visual_angle_in_degrees = 5; 
dist_to_screen_cm = 76;

[screen_size_mm(1), screen_size_mm(2)] = Screen('DisplaySize', 0);
screen_width_cm = screen_size_mm(1)/10;
screen_height_cm = screen_size_mm(2)/10;
[screen_width_pixels,screen_height_pixels] = Screen('WindowSize', 0);

screen_pixel_density = screen_width_pixels/screen_width_cm;
stim_width = 2*dist_to_screen_cm*tan((visual_angle_in_degrees/2)*(pi/180));
stim_pixels = stim_width*screen_pixel_density;
stim_radius = ceil(stim_pixels/2);

[circlespace,circlespace_frame,meshX,meshY] = make_circlespace(stim_radius);
 
x = meshX;
y = meshY;

% for i = 1:length(stim)

for i = 1:5
    
    % We use the fact that sin(x) has 1 cycle every 2*pi pixels.  That
    % means that sin(k*x) has k cycles every 2*pi pixels.  Using a conversion
    % factor of stim_pixels/visual_angle_in_degrees and after doing some
    % algebra we see that the following value of k will produce the
    % correct number of cycles per degree.

    cycles_per_degree = stim(i,2);
    k = (2*pi*visual_angle_in_degrees*cycles_per_degree)/stim_pixels;

    % Next we use that the gradient of the curve z = sin(theta)*x+cos(theta)*y
    % has constant length, meaning that the spatial frequency of sin(k*x) will
    % be the same as the spatial frequency of
    % sin(k*(sin(theta)*x+cos(theta)*y)) for all values of theta.

    orientation_angle = stim(i,3);
    g = sin(orientation_angle);
    h = cos(orientation_angle);
    wave=sin((g*x+h*y)*k);

    % Create image to place on screen

    disc = circlespace.*wave;
    disc = disc + circlespace_frame; 
    
    % get cat label
    cat_label = stim(i,1);
    cpd = stim(i,2);
    rad = stim(i,3);
    
    % File name has the following format: cat_fileindex.tif
    FILENAME = [pwd '/output/' num2str(i) '_' num2str(cat_label) '_' num2str(cpd) '_' num2str(rad) '.tif'];
    FMT = 'tiff';
    IMWRITE(disc,FILENAME,FMT, 'Resolution', 2)

end