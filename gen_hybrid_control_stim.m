% Create the control stimuli for the hybrid experiment.

load_cmd = ['load ' cd '/output/hybrid_stim.dat'];
eval(load_cmd);

data = hybrid_stim;

% Transform stim into 0-100 space

data(:,2) = (data(:,2)-0.25)*30;
data(:,3) =  data(:,3)*(200/pi);

clear above_data below_data

stim_per_category = 300;
above_count = 0;
below_count = 0;

for i = 1:size(data,1)
    if (data(i,3) > data (i,2))
        above_count = above_count + 1;
        above_data(above_count,:) = data(i,2:3);      
    else
        below_count = below_count + 1;
        below_data(below_count,:) = data(i,2:3);
    end
end

above_indices = randperm(above_count);

above_data = above_data(above_indices(1:stim_per_category),:);

below_indices = randperm(below_count);

extra_below_data = below_data(below_indices(1:(stim_per_category-below_count)),:);

below_data = [below_data(randperm(size(below_data,1)),:); extra_below_data];

figure           
hold on
xyaxes=[0 100 0 100];
plot(above_data(:,1),above_data(:,2),'*','markeredgecolor','r')
plot(below_data(:,1),below_data(:,2),'*','markeredgecolor','g')
legend('A','B')
hold off

above_data = [ones(size(above_data,1),1) above_data];
below_data = [2*ones(size(below_data,1),1) below_data];
final_stim = [above_data; below_data];
final_stim = final_stim(randperm(size(final_stim,1)),:);

final_stim(:,2) = .25 + (final_stim(:,2)/30); %x - axis (cpd)
final_stim(:,3) = final_stim(:,3) * (pi/200); % y - axis (rad)

% Plot sorted stim
figure('Position',[715 25 560 420]);
hold on
plot2dstim(final_stim,[0.25 3.58 0 pi/2],0)
axis on
legend('A','B')
hold off

open_cmd = ['fid = fopen([cd ''/output/hybrid_control_stim.dat''],''w'');'];
eval(open_cmd);
fprintf(fid,'%i %7.4f %1.4f \n',final_stim');
fclose(fid);