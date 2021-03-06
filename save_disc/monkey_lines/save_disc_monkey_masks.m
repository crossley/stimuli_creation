clear all
close all
clc

rand('seed',sum(100*clock));

% Category Distributions come from:
% Ashby, F. G., Ell, S. W., & Waldron, E. M. (2003)
% 832X624 resolution on a 15" CRT monitor

load II_line_stimuli.dat

% max_length = max(II_line_stimuli(:,2));
% min_length = min(II_line_stimuli(:,2));
% 
% max_angle = max(II_line_stimuli(:,3));
% min_angle = min(II_line_stimuli(:,3));

max_length = 3.5;
min_length = 1;

max_angle = 3.5;
min_angle = 2.0;

for i = 1:600
    
    % Generate five random line coordinates
    x = max_length*rand(15,2) + min_length;
    y = max_angle*rand(15,2) + min_angle;
    
%     figure, hold
%     plot(x(1,:),y(1,:),'k','Linewidth', 10, 'LineSmoothing', 'on')
%     plot(x(2,:),y(2,:),'k','Linewidth', 10, 'LineSmoothing', 'on')
%     plot(x(3,:),y(3,:),'k','Linewidth', 10, 'LineSmoothing', 'on')
%     plot(x(4,:),y(4,:),'k','Linewidth', 10, 'LineSmoothing', 'on')
%     plot(x(5,:),y(5,:),'k','Linewidth', 10, 'LineSmoothing', 'on')
    
    figure,hold
    set(gcf,'PaperPositionMode','manual')
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperPosition',[0 0 3.5 3.5])
%     set(gcf,'Renderer','painters')
    plot(x(1,:),y(1,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(2,:),y(2,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(3,:),y(3,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(4,:),y(4,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(5,:),y(5,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(6,:),y(6,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(7,:),y(7,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(8,:),y(8,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(9,:),y(9,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(10,:),y(10,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(11,:),y(11,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(12,:),y(12,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(13,:),y(13,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(14,:),y(14,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    plot(x(15,:),y(15,:),'k','Linewidth', 5, 'LineSmoothing', 'on')
    set(gca,'Unit','normalized','Position',[0 0 1 1]) 
%     axis([0 3.5 0 3.5])
    axis square
    axis off

    FILENAME = [pwd '/output/masks/' three_digit_string(i) '.bmp'];
    print_cmd = ['print -r100 -dbmp  ' FILENAME];
    eval(print_cmd); 

    close
   
end

% close all