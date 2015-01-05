% Generate stimuli for Hybrid study
% Part of the stimulus space is explicit 
% Part of the stimulus space is implicit

% The stimuli created below are designed to penalize rule use and have the
% best fitting unidimensional rule in the spatial frequency dimension.

% home

close all;
clear all;

n=5000;
raw_width = 100;
total_stim = 600; % Supposed to be a multiple of 10!
stim_per_type = [total_stim/2 total_stim/10 4*total_stim/10];

rand('seed',sum(100*clock));

% Generate a uniform distrubtion of numbers
x=round((rand(n,1))*raw_width);
y=round((rand(n,1))*raw_width);

stim_count(1:6) = 0;

for i=1:n
    if (x(i) < raw_width/2)
        if (y(i) > x(i))
            stim_count(1) = stim_count(1) + 1;
            stim_1(stim_count(1),:) = [x(i) y(i)];
        end
        if (y(i) < x(i))
            stim_count(2) = stim_count(2) + 1;
            stim_2(stim_count(2),:) = [x(i) y(i)];
        end
    end
    if (x(i) > raw_width/2)
        stim_count(3) = stim_count(3) + 1;
        stim_3(stim_count(3),:) = [x(i) y(i)];
    end
end

stim_1 = stim_1(1:stim_per_type(1),:);
stim_2 = stim_2(1:stim_per_type(2),:);
stim_3 = stim_3(1:stim_per_type(3),:);

figure           
hold on
xyaxes=[0 raw_width 0 raw_width];
plot(stim_1(:,1),stim_1(:,2),'*','markeredgecolor','r')
plot(stim_2(:,1),stim_2(:,2),'*','markeredgecolor','g')
plot(stim_3(:,1),stim_3(:,2),'*','markeredgecolor','b')
legend('A','B','B')
hold off

stim_1 = [ones(stim_per_type(1),1) stim_1]; % Category A member
stim_2 = [2*ones(stim_per_type(2),1) stim_2]; % Category B member
stim_3 = [2*ones(stim_per_type(3),1) stim_3]; % Category B member

temp_exp_stim = [stim_1; stim_2; stim_3];

% Transform stim into useful form

final_stim(:,1) = temp_exp_stim(:,1);
final_stim(:,2) = .25 + (temp_exp_stim(:,2)/30); %x - axis (cpd)
final_stim(:,3) = temp_exp_stim(:,3) * (pi/200); % y - axis (rad)

final_stim=randrows(final_stim);

% Plot sorted stim
figure('Position',[715 25 560 420]);
hold on 
plot2dstim(final_stim,[0.25 3.58 0 pi/2],0)
axis on
legend('A','B')
hold off

%  Write data to file
open_cmd = ['fid = fopen([cd ''/output/hybrid_stim.dat''],''w'');'];
eval(open_cmd);
fprintf(fid,'%i %7.4f %1.4f \n',final_stim');
fclose(fid);

num_a_members = length(find(final_stim(:,1) == 1))
num_b_members = length(find(final_stim(:,1) == 2))