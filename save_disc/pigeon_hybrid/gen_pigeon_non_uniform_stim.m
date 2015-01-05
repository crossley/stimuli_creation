close all
clear all
clc

region_I_stim = [ones(75,1) 50*rand(75,1) (50*rand(75,1)+50)];
region_II_stim = [2*ones(75,1) (50*rand(75,1)+50) (100*rand(75,1))];
region_III_stim = [50*rand(450,1) 50*rand(450,1)];

A_ind = find(region_III_stim(:,1) < region_III_stim(:,2));
B_ind = find(region_III_stim(:,1) >= region_III_stim(:,2));

region_III_A_stim = region_III_stim(A_ind,:);
region_III_B_stim = region_III_stim(B_ind,:);

region_III_A_stim = [ones(size(region_III_A_stim,1),1) region_III_A_stim];
region_III_B_stim = [2*ones(size(region_III_B_stim,1),1) region_III_B_stim];

new_hybrid_stim = [region_I_stim; region_II_stim; region_III_A_stim; region_III_B_stim];

plot2dstim(new_hybrid_stim)
hold on
plot([0:50 50*ones(1,50)],[0:100],'-k','linewidth', 2)
plot([0:100],[0:100],'-k','linewidth', 2)
hold off
% exportfig(gcf, 'new_hybrid_stim.png', 'format', 'png', 'width', 5, 'height', 5, 'resolution', 720)

num_stim = size(new_hybrid_stim,1);

final_stim = new_hybrid_stim;

final_stim(:,1) = new_hybrid_stim(:,1);
final_stim(:,2) = .25 + (new_hybrid_stim(:,2)/30); %x - axis (cpd)
final_stim(:,3) = (new_hybrid_stim(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

fid = fopen('new_hybrid_stim.dat','w');
fprintf(fid,'%i %f %f \n',final_stim');
fclose(fid);