function [circlespace,circlespace_frame,meshX,meshY] = make_circlespace(r)
% gives a circle with radius r 

h = 2*r+1;

circlespace = zeros(h,h);
circlespace_frame = zeros(h,h);

center_point = r+1;

for i = 1:h

    for j = 1:h
    
        if (((i-center_point)^2 + (j-center_point)^2) <= (r/2)^2)

            circlespace(i,j) = 1;

        else
            
            circlespace_frame(i,j) = 1;
        
        end
    end

end

[meshX,meshY] = meshgrid(1:h,1:h);