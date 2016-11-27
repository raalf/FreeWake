%=========================================================================%
%simple tool that creates a wakeplot in matlab.
%
%reads: timestep files
%uses: Main_FindDVEs.exe
%creates: a figure
%all files must be located in the output folder
%axis, view, size, position and framerate are set at the end of the file
%%========================================================================%
function [] = wakeplot_new(step,sym)
if nargin == 0
    
    
    clear
    clc
    step = input('Which timestep to plot?\n');
    
    sym = input('Do you want symmetry?(y/n)\n','s');
    %run wakeplot and get axis settings (to get full axis)
end

h = figure(2);
hold on
% clf(1); %this clf might cause an error if we are making a gif!!!!!!
%  set(gca,'Color',[0 0 0])
% set(gca,'CameraViewAngleMode','manual')
%  set(gcf,'renderer','zbuffer')
%opengl('software')



%Check to see if timestep file exists
check = strcat('timestep',num2str(step),'.txt');
if exist(check,'file') ==0
    error('Timestep file not found.')
end

%create input file for DVE generator
fid = fopen( 'temp.txt', 'wt' );
fprintf (fid,'%d',step);

%run Main_FindDVEs to get coordinates of all DVEs for this timestep
[out,cmdout] =system('Main_FindDVEs0.exe < temp.txt');

%open DVEs.dat
fp = fopen('DVEs.dat','r');

if fp == -1
    error('Cannot open the output of Main_FindDVEs.exe')
end
%look for timestep
loc=0;
while strcmp(loc,'timestep') ==0
    loc =fscanf(fp,'%s',1);
end
timestep = fscanf(fp,'%d',1);

%look for number of wings
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
numwings = fscanf(fp,'%d',1);

%look for n
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
n = fscanf(fp,'%d',1);

%look for m
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
m = fscanf(fp,'%d',1);

numelements = n*m;

%look for wake
loc=0;
while strcmp(loc,'%WAKE') ==0
    loc =fscanf(fp,'%s',1);
end

data= zeros(numwings*timestep*n,13);

%find ZZZZ
loc=0;
while strcmp(loc,'ZZZZ') ==0
    loc =fscanf(fp,'%s',1);
end

%read wake data
for x = 1:(numwings*(timestep+1)*n)
    %read wing
    data(x,1) = fscanf(fp,'%d' ,1);
    %skip ,
    fseek(fp,1,0);
    %read x coords of each corner
    data(x,2:5) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read y coords of each corner
    data(x,6:9) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read z coords of each corner
    data(x,10:13) = fscanf(fp,'%lf %lf %lf %lf' ,4);
end
data = flipud(data);

%find WING
loc=0;
while strcmp(loc,'%WING') ==0
    loc =fscanf(fp,'%s',1);
end

%find ZZZZ
loc=0;
while strcmp(loc,'ZZZZ') ==0
    loc =fscanf(fp,'%s',1);
end
wing = zeros(numelements*numwings,12);

%read wing data
for x = 1:(numelements*numwings)
    %read x coords of each corner
    wing(x,1:4) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read y coords of each corner
    wing(x,5:8) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read z coords of each corner
    wing(x,9:12) = fscanf(fp,'%lf %lf %lf %lf' ,4);
end

% fclose('all');
delete('./temp.txt');
delete('./DVEs.dat');

%colors
% cmap = jet(numwings);
cmap = parula(100);
% cmap = flipud(cmap);
colormap(cmap);
%plot wing

hold on

%plot wing
% x = [1:(numelements*numwings)];
patch(wing(:,1:4)',wing(:,5:8)',wing(:,9:12)',[0.4 0.4 0.4],'FaceAlpha',0.8,'EdgeColor','k');
if strcmp(sym,'y')==1
    patch(wing(:,1:4)',-wing(:,5:8)',wing(:,9:12)',[0.4 0.4 0.4],'FaceAlpha',0.8,'EdgeColor','k');
end
hold off


hold on
%plot wakes
count=1;

% surf(reshape(data(:,2),n,timestep),reshape(data(:,6),n,timestep),reshape(data(:,10),n,timestep));
%color by wake number
% id=(data(:,1)==3 );
patch(data(:,2:5)',data(:,6:9)',data(:,10:13)',(data(:,1))','FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.7)
    if strcmp(sym,'y')==1
       patch(data(:,2:5)',-data(:,6:9)',data(:,10:13)',(data(:,1))','FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.7)
    end

%color by parameter (in progress)
% patch(data(:,2:5)',data(:,6:9)',data(:,10:13)',data(:,6:9)','FaceAlpha',0.5,'EdgeColor','k')
%     if strcmp(sym,'y')==1
%        patch(data(:,2:5)',-data(:,6:9)',data(:,10:13)',data(:,6:9)','FaceAlpha',0.5,'EdgeColor','k')
%     end


hold off
%modify axis
axis equal
axis tight
%AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX])
% axis([xlimit ylimit zlimit]);
%  axis([-300 1600 0 1600 -1900 400]);

%VIEW(AZ,EL)

view([-23 16]);

%set frame position and size

set(h,'renderer','painters');
%ZOOM(factor)
%     zoom(1)

%record frame for movie
%     mov(step) = getframe(gcf);


% j=figure(2)
% clf(2)
% hold on
% for x = 1:(numelements*numwings)
%     fill3(wing(x,1:4),wing(x,5:8),wing(x,9:12),'-k','FaceAlpha',0.6,'EdgeColor','k');
%     if strcmp(sym,'y')==1
%         fill3(wing(x,1:4),-wing(x,5:8),wing(x,9:12),'-k','FaceAlpha',0.6,'EdgeColor','k');
%     end
% end    
% 
% 
% for instep = 1:timestep
%     for x = count:(count+(numwings*n))-1
%         fill3(data(x,2:5),data(x,6:9),data(x,10:13),cmap(data(x,1),:),'FaceAlpha',0.5,'EdgeColor','k')
%         if strcmp(sym,'y')==1
%             fill3(data(x,2:5),-data(x,6:9),data(x,10:13),cmap(data(x,1),:),'FaceAlpha',0.5,'EdgeColor','k')
%         end
%     end
%     count = count+(numwings*n);
% end
% 
% hold off
% %modify axis
% axis('equal')
% 
% %AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX])
% % axis([xlimit ylimit zlimit]);
% %  axis([-300 1600 0 1600 -1900 400]);
% 
% %VIEW(AZ,EL)
% 
% view([-23 16]);
% 
% %set frame position and size
% 
% set(j,'renderer','painters');
% %ZOOM(factor)
% %     zoom(1)
% 
% %record frame for movie
% %     mov(step) = getframe(gcf);
% 
ax = gca;
z = zoom;
setAxes3DPanAndZoomStyle(z,ax,'camera')

fclose all

