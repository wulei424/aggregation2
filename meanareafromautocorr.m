function [A, radial_average, radii_list, Cnorm]=meanareafromautocorr(I,maxr)

  
if nargin<2
maxr = 50; % round(0.4 * size(I,1));
end

% compute autocorrelation
R=xcorr2(I);

% Create the meshgrid to be used in resampling
[X,Y] = meshgrid(1:size(R,1),1:size(R,2));

center = (size(R)+3)/2;

dr=1;
radii_list = 1:dr:maxr;

for r = 1:length(radii_list);
    
    radius=radii_list(r);
    
    % To avoid redundancy, sample at roughly 1 px distances
    num_pxls = 2*pi*radius;
    theta = 0:1/num_pxls:2*pi;
    
    x = center(1) + radius*cos(theta);
    y = center(2) + radius*sin(theta);


    sampled_radial_slice = interp2(X,Y,R,x,y);
    radial_average(r) = mean(sampled_radial_slice);
    
end

Cnorm = max(R(:));
r=find(radial_average./Cnorm < 1/exp(1), 1, 'First');

if isempty(r)
    warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
    A= pi* radii_list(end)^2;
else
    A= pi* radii_list(r)^2;
end
    

