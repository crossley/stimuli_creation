clear all
close all
clc

rand('seed',sum(100*clock));

% yes bivariate normals. the variance along x and y is 100 and the covariance is zero
% 
% the means are as follows
%        x       y
% A       72      100
% B       100     128
% C       100     72
% D       128     100

% Number of As and Bs

total_stim_num = 600;

if (mod(total_stim_num,2) ~= 0)

    error('Need to have equal number of stim per category')
    
end    
    
stim_per_category = total_stim_num/4;
range = 100;

x_c = 43;
y_c = range-x_c;

% mu_a = [x_c,y_c];
% mu_b = [y_c,x_c];

mu_a = [72 100];
mu_b = [100 128];
mu_c = [100 72];
mu_d = [128 100];

% mu_a = [32+20,48+20];
% mu_b = [48+20,32+20];

% mu_a = [55 50];
% mu_b = [47 50];
 
cov_a = [100 0
          0 100];
      
cov_b = cov_a;
cov_c = cov_a;
cov_d = cov_a;

d1 = sqrt(2)*x_c/3; 
d2 = sqrt(2)*(y_c-x_c)/7;

% d1 = sqrt(2)*32/3; 
% d2 = sqrt(2)*(48-32)/7;

% cov_a = make_covar_matrix(d1^2,d2^2);

% cov_a(2,1) = -cov_a(2,1);
% cov_a(1,2) = -cov_a(1,2);

% cov_b = cov_a;
 
sample_a = make_mult_norm(mu_a,cov_a,stim_per_category);
sample_b = make_mult_norm(mu_b,cov_b,stim_per_category);
sample_c = make_mult_norm(mu_c,cov_c,stim_per_category);
sample_d = make_mult_norm(mu_d,cov_d,stim_per_category);

labeled_sample_a = [ones(length(sample_a),1) sample_a];
labeled_sample_b = [2*ones(length(sample_b),1) sample_b];
labeled_sample_c = [3*ones(length(sample_c),1) sample_c];
labeled_sample_d = [4*ones(length(sample_d),1) sample_d];

xyaxes = [0 200 0 200];
data = [labeled_sample_a;labeled_sample_b;labeled_sample_c;labeled_sample_d];
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
    sample_c = transample(sample_c,mu_c,cov_c);
	sample_d = transample(sample_d,mu_d,cov_d);

    clipped_sample_a = [];
	clipped_sample_b = [];
    clipped_sample_c = [];
	clipped_sample_d = [];
	    
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
        
        if (mahal_dist(sample_c(i,:), mu_c, cov_c) <= 3)
            
            clipped_sample_c(size(clipped_sample_c,1)+1,:) = sample_c(i,:);
            
        else
            
            run_again = 1;
            
         end
        
        if (mahal_dist(sample_d(i,:), mu_d, cov_d) <= 3)
            
            clipped_sample_d(size(clipped_sample_d,1)+1,:) = sample_d(i,:);
            
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
    
    while (size(clipped_sample_c,1) < stim_per_category)
	
        temp_stim = make_mult_norm(mu_c,cov_c,1);
	
        if (mahal_dist(sample_c(i,:), mu_c, cov_c) <= 3)
            
            clipped_sample_c(size(clipped_sample_c,1)+1,:) = temp_stim;
            
        end
        
    end
    
    while (size(clipped_sample_d,1) < stim_per_category)
	
        temp_stim = make_mult_norm(mu_d,cov_d,1);
	
        if (mahal_dist(sample_d(i,:), mu_d, cov_d) <= 3)
            
            clipped_sample_d(size(clipped_sample_d,1)+1,:) = temp_stim;
            
        end
        
	end
    
    sample_a = clipped_sample_a;
    sample_b = clipped_sample_b; 
    sample_c = clipped_sample_c;
    sample_d = clipped_sample_d; 
    
   	labeled_sample_a = [ones(stim_per_category,1) sample_a];
	labeled_sample_b = [2*ones(stim_per_category,1) sample_b];
    labeled_sample_c = [3*ones(stim_per_category,1) sample_c];
	labeled_sample_d = [4*ones(stim_per_category,1) sample_d];
	
	data = [labeled_sample_a;labeled_sample_b;labeled_sample_c;labeled_sample_d];
	figure('Position',[50 550 560 420],'MenuBar','none','NumberTitle','off','Name',['Emerging Stim - Iteration: ' num2str(run_number)] );
	plot2dstim(data,xyaxes,0);
    
end

% Change to rule-based category labels
% A_ind = find(data(:,2)>50);
% B_ind = find(data(:,2)<=50);
% 
% data(A_ind,1) = 1;
% data(B_ind,1) = 2;



% final_stim(:,1) = data(:,1);
% final_stim(:,2) = .25 + (data(:,2)/30); %x - axis (cpd)
% final_stim(:,3) = (data(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

final_stim(:,1) = data(:,1);
final_stim(:,2) = 1.0 + (data(:,2)/30); %x - axis (cpd)
final_stim(:,3) = (data(:,3) * (pi/200)) + (pi/9); % y - axis (rad)

% xyaxes = [0.25 3.58 pi/9 pi/2+pi/9];
xyaxes = [0.0 10.0 0.0 4.0];
% xyaxes = [];

% Plot Original Data
figure('Position',[50 50 560 420],'MenuBar','none','NumberTitle','off','Name','Final Stim');
plot2dstim(final_stim,xyaxes,0)

% final = [category orientation spatial_frequency]

final_stim = randrows(final_stim);

%  Write data to file
% open_cmd = ['fid = fopen([cd ''/output/' num2str(total_stim_num) '_ii_stim.dat''],''w'');'];
% open_cmd = ['fid = fopen([cd ''/output/RB_II_Train_stim_unix.dat''],''w'');'];
% eval(open_cmd);
% fprintf(fid,'%i %7.4f %1.4f \n',final_stim');
% fclose(fid);

