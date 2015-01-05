clear all
close all
clc

rand('seed',sum(100*clock));

% Category Distributions come from:
% Ashby, F. G., Ell, S. W., & Waldron, E. M. (2003)
% 832X624 resolution on a 15" CRT monitor

% Number of As and Bs

total_stim_num = 600;

if (mod(total_stim_num,2) ~= 0)

    error('Need to have equal number of stim per category')
    
end    
    
stim_per_category = total_stim_num/2;

mu_a = [280 350];
mu_b = [420 350];

cov_a = [330 0
         0 14000];

cov_b = cov_a;
 
sample_a = make_mult_norm(mu_a,cov_a,stim_per_category);
sample_b = make_mult_norm(mu_b,cov_b,stim_per_category);

labeled_sample_a = [ones(length(sample_a),1) sample_a];
labeled_sample_b = [2*ones(length(sample_b),1) sample_b];

data = [labeled_sample_a;labeled_sample_b;];

% figure
% plot2dstim(data,[0 700 0 700],0);

% We now enter a loop in which we first run transform_sample to make sure that the population 
% mean is equal to the sample mean and that the population covariance is equal to the sample
% variance.  We then remove all stimuli whose Mahalanobis distance from the mean is greater than 3.
% If we remove any stimuli then we run the loop again.

run_again = 1;
run_number = 0;

while (run_again)
    
    run_number = run_number + 1;
    run_again = 0;
    
	% Transform the sample so that the sample means and covariance
	% matrices are equal to the population means and covariance matrices.
	
	sample_a = transample(sample_a,mu_a,cov_a);
	sample_b = transample(sample_b,mu_b,cov_b);

    clipped_sample_a = [];
	clipped_sample_b = [];
	    
    % We now clip out the stimuli that are too far from the mean.  If there
    % are any then we will run the loop again
    
	for i = 1:stim_per_category
       
        if (mahal_dist(sample_a(i,:), mu_a, cov_a) <= 3)
	
            clipped_sample_a(size(clipped_sample_a,1)+1,:) = sample_a(i,:);
            
        else
            
            run_again = 1;
            
        end
        
        if (mahal_dist(sample_b(i,:), mu_b, cov_b) <= 3)
            
            clipped_sample_b(size(clipped_sample_b,1)+1,:) = sample_b(i,:);
            
        else
            
            run_again = 1;
            
        end
        
	end
	
	% We add stimuli to make up for those lost in clipping. 
	
	while (size(clipped_sample_a,1) < stim_per_category)
	
        temp_stim = make_mult_norm(mu_a,cov_a,1);
	
        if (mahal_dist(sample_a(i,:), mu_a, cov_a) <= 3)
            
            clipped_sample_a(size(clipped_sample_a,1)+1,:) = temp_stim;
            
        end
        
	end
	
	while (size(clipped_sample_b,1) < stim_per_category)
	
        temp_stim = make_mult_norm(mu_b,cov_b,1);
	
        if (mahal_dist(sample_b(i,:), mu_b, cov_b) <= 3)
            
            clipped_sample_b(size(clipped_sample_b,1)+1,:) = temp_stim;
            
        end
        
    end
    
    sample_a = clipped_sample_a;
    sample_b = clipped_sample_b; 
    
   	labeled_sample_a = [ones(stim_per_category,1) sample_a];
	labeled_sample_b = [2*ones(stim_per_category,1) sample_b];
	
	data = [labeled_sample_a;labeled_sample_b;];
end

% figure
% plot2dstim(data,[0 700 0 700],0);

% pixel_density = get(0,'ScreenPixelsPerInch'); 
pixel_density = 0.017;

transformed_data = data;
transformed_data(:,1) = data(:,1);
transformed_data(:,2) = data(:,2)*(1/255); % covert to inches
transformed_data(:,3) = data(:,3)*(pi/700) + pi/6; % convert to radians

% Inverse Transform
% % % abstract_stim(:,2) = II_line_stimuli(:,2)*255;
% % % abstract_stim(:,3) = (II_line_stimuli(:,3)-pi/6)*700/pi;

% figure
% plot2dstim(transformed_data)

for i = 1:600
    
    horizontal_range = 3.5;
    horizontal_center = horizontal_range/2;
        
    vertical_range = 3.5;
    vertical_center = vertical_range/2;
    
    length = transformed_data(i,2);
    horizontal_length = length*cos(transformed_data(i,3));
    vertical_length = length*sin(transformed_data(i,3));
    
    x0 = horizontal_center - 0.5*horizontal_length;
    x1 = horizontal_center + 0.5*horizontal_length;
    
    y0 = vertical_center - 0.5*vertical_length;
    y1 = vertical_center + 0.5*vertical_length;
    
    % Define line
    x = [x0 x1];
    y = [y0 y1];
    
    figure,hold
    set(gcf,'PaperPositionMode','manual')
    set(gcf,'PaperUnits','inches')
%     set(gcf,'PaperSize',[1 1])
    set(gcf,'PaperPosition',[0 0 3.5 3.5])
%     set(gcf,'Renderer','painters')
    plot(x, y, 'k', 'Linewidth', 10, 'LineSmoothing', 'on')
    set(gca,'Unit','normalized','Position',[0 0 1 1]) 
    axis([0 3.5 0 3.5])
    axis square
    axis off

    cat_label = transformed_data(i,1);
    length = round(data(i,2));
    rad = round(data(i,3));

    FILENAME = [pwd '/output/stim/' three_digit_string(i) '_' three_digit_string(cat_label) '_' three_digit_string(length) '_' three_digit_string(rad) '.bmp'];
    print_cmd = ['print -r100 -dbmp '  FILENAME];
    eval(print_cmd); 

    close
   
end

fid = fopen('RB_line_stimuli.dat','w');
fprintf(fid,'%i %f %f \n',transformed_data');
fclose(fid);

% close all