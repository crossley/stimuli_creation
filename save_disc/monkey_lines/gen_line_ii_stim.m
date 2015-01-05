clear all
close all
clc

rand('seed',sum(100*clock));

% Number of As and Bs

total_stim_num = 600;

if (mod(total_stim_num,2) ~= 0)

    error('Need to have equal number of stim per category')
    
end    
    
stim_per_category = total_stim_num/2;
range = 350;


mu_a = [150 220];
mu_b = [220 150];


cov_a = [30^2 27^2
          27^2 30^2];


cov_b = cov_a;
 
sample_a = make_mult_norm(mu_a,cov_a,stim_per_category);
sample_b = make_mult_norm(mu_b,cov_b,stim_per_category);

labeled_sample_a = [ones(length(sample_a),1) sample_a];
labeled_sample_b = [2*ones(length(sample_b),1) sample_b];

xyaxes = [0 350 0 350];
data = [labeled_sample_a;labeled_sample_b;];
figure('Position',[50 550 560 420],'MenuBar','none','NumberTitle','off','Name','Raw Stim');
plot2dstim(data,xyaxes,0);

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
	figure('Position',[50 550 560 420],'MenuBar','none','NumberTitle','off','Name',['Emerging Stim - Iteration: ' num2str(run_number)] );
	plot2dstim(data,xyaxes,0);
    
end
