%=========================================================================%
%simple tool that creates a wake animation
%
%reads: timestep files
%uses: Main_FindDVEs.exe
%creates: wakeframes.avi
%all files must be located in the output folder
%axis, view, size, position and framerate are set at the end of the file
%%========================================================================%
clear
laststep = input('Which timestep to plot up to?\n');

%run wakeplot and get axis settings (to get full axis)
figure(1);
run wakeplot_color.m

xlimit = xlim;
ylimit = ylim;
ylimit(2) = ylimit(2) + ((ylimit(2)-ylimit(1))/10);
zlimit = zlim;
zlimit(2) = zlimit(2) + abs((zlimit(2)-zlimit(1))/5);
clf(1);
%create and set up figure
writerObj = VideoWriter('wakeframes.avi');
writerObj.FrameRate = 5;
open(writerObj);


h = figure(1);
%  set(gca,'Color',[0 0 0])
set(gca,'CameraViewAngleMode','manual')
%  set(gcf,'renderer','zbuffer')
%opengl('software')

for step = 1:laststep
    
    %Check to see if timestep file exists
    check = strcat('timestep',num2str(step),'.txt');
    if exist(check,'file') ==0
        error('Timestep file not found. All timestep files up to the request need to be availble.')
    end
    
    %create input file for DVE generator
    fid = fopen( 'temp.txt', 'wt' );
    fprintf (fid,'%d',step);
    
    %run Main_FindDVEs to get coordinates of all DVEs for this timestep
    system('Main_FindDVEs.exe < temp.txt');
    
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
    for x = 1:(numwings*timestep*n)
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
    
    fclose('all');
    delete('./temp.txt');
    delete('./DVEs.dat');
    
    %colors
    cmap = jet(numwings);
    %     cmap = jet(laststep);
    cmap = flipud(cmap);
    
    %plot wing
    hold on
    for x = 1:(numelements*numwings)
        fill3(wing(x,1:4),wing(x,5:8),wing(x,9:12),'-k','FaceAlpha',0.6,'EdgeColor','k')
    end
    
    %plot wakes
    count=1;
    for instep = 1:timestep
        for x = count:(count+(numwings*n))-1
            fill3(data(x,2:5),data(x,6:9),data(x,10:13),cmap(data(x,1),:),'FaceAlpha',0.5,'EdgeColor','k')
        end
        count = count+(numwings*n);
    end
    
    %modify axis
    axis('equal')
    
    %AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX])
    axis([xlimit ylimit zlimit]);
    %  axis([-300 1600 0 1600 -1900 400]);
    
    %VIEW(AZ,EL)
    view([-23 16])
    
    %ZOOM(factor)
%     zoom(1)
    
    %set frame position and size
    set(h,'Position',[300 100 1200 600])
    
    %record frame for movie
    %     mov(step) = getframe(gcf);
    frame = getframe;
    writeVideo(writerObj,frame);
    
    %set up figure for next iteration
    clf
end
close(writerObj);

