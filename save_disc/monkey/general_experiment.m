
% General Setup 

clc 
clear all
close all

% We start by naming the experiment.  *** Important - The program expects
% the stimulus set to be experiment_name_stim.dat ***

experiment_name = 'general';

warning('off')
PsychJavaTrouble 
screen('Preference', 'SkipSyncTests',1);

% Setup for experiment

global window screenRect white black gray xc yc meshX meshY circlespace visual_angle_in_degrees stim_pixels

% Get Subject Number
sub = input(' Subject Number : ');
day = input(' Day Number : ');

% Create data file
outfile = [cd '/data/' experiment_name '_sub_' num2str(sub) '_day_' num2str(day) '.dat'];

if exist(outfile)
 	error(['The file ' outfile ' already exists.']);
end

fid = fopen([outfile],'w');

% We specify the number of practice trials, the number of trials per block and which blocks will be
% feedback blocks.

num_practice_trials = 5;
trials_per_block = 50;

% We specify the desired degrees of visual angle for the stimulus as well
% as how far the subject will be from the screen.

visual_angle_in_degrees = 5; 
dist_to_screen_cm = 76; 

if (day == 1)
    
    training_blocks = 1:12;
    extinction_blocks = [];

else

    training_blocks = 1:12;
    extinction_blocks = [];

end

total_blocks = length(training_blocks) + length(extinction_blocks);
num_trials = total_blocks*trials_per_block;

% We open the screen

[window,screenRect] = screen('OpenWindow',0,0);
hideCursor;
ListenChar(2); 

% Define colors

white = WhiteIndex(window);
black = BlackIndex(window);
gray = [(white+black)/2 (white+black)/2 (white+black)/2];

% Define center

xc = screenRect(3)/2; 
yc = screenRect(4)/2; 

% Obtain screen size and pixel info, then create an appropriate circlespace
% matrix

[screen_size_mm(1), screen_size_mm(2)] = Screen('DisplaySize', 0);
screen_width_cm = screen_size_mm(1)/10;
screen_height_cm = screen_size_mm(2)/10;
[screen_width_pixels,screen_height_pixels] = Screen('WindowSize', 0);

screen_pixel_density = screen_width_pixels/screen_width_cm;
stim_width = 2*dist_to_screen_cm*tan((visual_angle_in_degrees/2)*(pi/180));
stim_pixels = stim_width*screen_pixel_density;
stim_radius = ceil(stim_pixels/2);
      
[circlespace,meshX,meshY] = make_circlespace(stim_radius);

feedback_init(); % prepares sound stimlui for auditory feedback
rand('seed',sum(100*clock))

% Load data

load_cmd = ['load ' cd '/input/' experiment_name '_stim.dat'];
eval(load_cmd)

def_cmd = ['stimuli = ' experiment_name '_stim;'];
eval(def_cmd)

% Begin experiment

begin_experiment;

begin_practice;

practice_stim = stimuli(1:num_practice_trials,:);

experiment_stim = stimuli(randperm(size(stimuli,1)),:);

for trial = 1:min(num_practice_trials, num_trials)
    
        clear_screen;
        show_crosshairs;
        clear_screen;
        show_disc(practice_stim(trial,:));
        [response rt] = get_response_with_feedback(practice_stim(trial,:));
        clear_screen;

end

end_practice;

% Training

for block_num = 1:total_blocks
    
    begin_block(block_num);
    
    if (isempty(find(extinction_blocks == block_num))) 

        for trial_num = 1:trials_per_block

            trial = trials_per_block*(block_num-1) + trial_num;

            % Show crosshairs and stimulus, then get response and give feedback

            clear_screen;
            show_crosshairs;
            clear_screen;
            show_disc(experiment_stim(trial,:));
            [response rt] = get_response_with_feedback(experiment_stim(trial,:));
            clear_screen;

            % Write trial stimuli response1 response2 rt1 rt2 to a data file.
            % data=[trial corrcat orientation spatialfreq response RT]' );

            data(trial,:) = [trial experiment_stim(trial,:) response rt];
            fprintf(fid,'%3i %i %3.2f %1.4f %i %1.4f\n',data(trial,:)' );


        end

    else
        
        for trial_num = 1:trials_per_block
    
        trial = trials_per_block*(block_num-1) + trial_num;
        
        % Show crosshairs and stimulus, then get response but don't give
        % feedback
                
        clear_screen;
        show_crosshairs;
        clear_screen;
        show_disc(experiment_stim(trial,:));
        [response rt] = get_response_no_feedback(experiment_stim(trial,:));
        clear_screen;
        
        % Write trial stimuli response1 response2 rt1 rt2 to a data file.
        % data=[trial corrcat orientation spatialfreq response RT]' );
        
        data(trial,:) = [trial experiment_stim(trial,:) response rt];
        fprintf(fid,'%3i %i %3.2f %1.4f %i %1.4f\n',data(trial,:)' );
        
        end
        
    end
        
    end_block(block_num);

end

end_experiment;

% close things up
ListenChar(0);
ShowCursor;
screen('CloseAll');
fclose(fid);