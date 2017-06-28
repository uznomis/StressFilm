% this program processes stress film image
%% import image and plot
[name,path] = uigetfile('*.tif');
if name == 0
    return
end
I = imread([path,name]);
grayI = rgb2gray(I);
figure;
imshow(I);

%% locate center of circle and radii
hold on
[x,y] = ginput(3);
[r,c] = fit_circle_through_3_points([x,y]);
[X,Y] = ginput(3);
[R,C] = fit_circle_through_3_points([X,Y]);
xc = (c(1)+C(1))/2;
yc = (c(2)+C(2))/2;
draw_circle(gca,xc,yc,r,'b');
draw_circle(gca,xc,yc,R,'b');
hold off

%% calculate intensity pixel by pixel and group by sector
sectorStep = pi/36;
pixelBySector = cell(1,round(2*pi/sectorStep));
currentProgress = 0;
for i = 1:size(I,1)
    progressPercent = floor(i/size(I,1)*100);
    if mod(progressPercent,10) == 0 && progressPercent ~= currentProgress
        disp([num2str(progressPercent),'%']);
        currentProgress = progressPercent;
    end
    for j = 1:size(I,2)
        dist = sqrt((j-xc)^2+(i-yc)^2); % i is y, j is x
        if r<dist && dist<R
            [theta,rho] = cart2pol(j-xc,i-yc);
            if theta < 0
                theta = theta + 2*pi;
            end
            pixelBySector{ceil(theta/sectorStep)} = ...
                [pixelBySector{ceil(theta/sectorStep)} grayI(i,j)];
        end
    end
end

%% shows mean intensity for each sector
load('stressFilmIntensities.mat');
stresses = zeros(length(pixelBySector),1);
imshow(I);
hold on
for i = 1:length(pixelBySector)
    [x1,y1] = pol2cart(sectorStep*i, r);
    [x2,y2] = pol2cart(sectorStep*i, R);
    [x2t,y2t] = pol2cart(sectorStep*(i-0.5), R);
    line([x1,x2]+xc,[y1,y2]+yc);
    intensity = mean(pixelBySector{i})*p(1) + p(2);
    stresses(i,1) = intensity;
    text(x2t+xc,y2t+yc,num2str(intensity));
end
hold off
disp('stresses are');
disp(stresses);
return

%% pick rectangle and calculate mean intensity
[x,y] = ginput(2);
rect = grayI(round(y(1)):round(y(2)),round(x(1)):round(x(2)));
hold on
line(x,y);
hold off
disp(mean2(rect));

%% regression to calculate y = a*x+b for conversion to stress
p = polyfit(intensities(:,1),intensities(:,2),1);