%=========================================================================%
%simple tool that creates a wakeplot in matlab.
%
%reads: timestep files
%uses: Main_FindDVEs2020.exe
%creates: a figure
%all files must be located in the output folder
%axis, view, size, position and framerate are set at the end of the file
%%========================================================================%
function [] = wakeplot_2020(step,sym,outpath)
if nargin == 0
    
    
    clear
    clc
    step = input('Which timestep to plot?\n');
    
    sym = input('Do you want symmetry?(y/n)\n','s');
    outpath = 'input';
    h = figure(33);
    clf(33);
    %run wakeplot and get axis settings (to get full axis)
end


hold on
h = figure(33);
% clf(33) might cause an error if we are making a gif!!!!!!

%  set(gca,'Color',[0 0 0])
% set(gca,'CameraViewAngleMode','manual')
%  set(gcf,'renderer','zbuffer')
%opengl('software')



%Check to see if timestep file exists
check = strcat(outpath,'/timestep',num2str(step),'.txt');
if exist(check,'file') ==0
    error('Timestep file not found.')
end

%create input file for DVE generator
fid = fopen( 'temp.txt', 'wt' );
fprintf (fid,'%d',step);
fprintf (fid,'\n%s',outpath);
fclose(fid);

%run Main_FindDVEs2020 to get coordinates of all DVEs for this timestep
[out,cmdout] =system('Main_FindDVEs2020.exe < temp.txt');

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

%look for nosurface
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
nosurface = fscanf(fp,'%d',1);

numelements = nosurface;

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
wing = zeros(numelements,12);

%read wing data
for x = 1:(numelements)
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
fclose(fp);
% fclose('all');
delete('./temp.txt');
delete('./DVEs.dat');

%colors
% cmap = jet(numwings);
cmap = parula(100);
% cmap = flipud(cmap);
colormap(cmap);
%plot wing

subplot(1,2,1)
hold on

%plot wing
% x = [1:(numelements*numwings)];
patch(wing(:,1:4)',wing(:,5:8)',wing(:,9:12)',[0.4 0.4 0.4],'FaceAlpha',0.8,'EdgeColor','k');
if strcmp(sym,'y')==1
    patch(wing(:,1:4)',-wing(:,5:8)',wing(:,9:12)',[0.4 0.4 0.4],'FaceAlpha',0.8,'EdgeColor','k');
end
hold on
elecg= ([mean(wing(:,1:4),2),mean(wing(:,5:8),2),mean(wing(:,9:12),2)]);
plot3(elecg(:,1),elecg(:,2),elecg(:,3),'m*');

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

%% add quiver if available
try
    tempdir = cd(outpath);
%quiver.txt
dataLines = [1, Inf]
opts = delimitedTextImportOptions("NumVariables", 6);
opts.DataLines = dataLines;
opts.Delimiter = "\t";
opts.VariableNames = ["x1", "y1", "z1", "u1", "v1", "w1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";



%read direction unit vectors from quiver1
direct = readtable('quiver1.txt', opts);
quiver3(direct.x1,direct.y1,direct.z1,direct.u1,direct.v1,direct.w1);

%read local freestream vel from quiver2
uvel = readtable('quiver2.txt', opts);
quiver3(uvel.x1,uvel.y1,uvel.z1,uvel.u1,uvel.v1,uvel.w1);

%read span force and moment from quiver3
forcemom = readtable('quiver3.txt', opts);
quiver3(forcemom.x1,forcemom.y1,forcemom.z1,forcemom.u1,forcemom.v1,forcemom.w1,0.1);

rotforcemom = readtable('quiver4.txt', opts);


totalforcemom = readtable('quiver5.txt', opts);

cd(tempdir);
clear tempdir
catch
    cd(tempdir);
clear tempdir
end
%%
hold off
%modify axis
axis equal
axis tight
%AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX])
% axis([xlimit ylimit zlimit]);
%  axis([-300 1600 0 1600 -1900 400]);

grid on
%VIEW(AZ,EL)

view([-23 16]);
xlabel('X');
ylabel('Y');
%set frame position and size

% set(h,'renderer','painters');
%ZOOM(factor)
%     zoom(1)

ax = gca;
z = zoom;
setAxes3DPanAndZoomStyle(z,ax,'camera')

%%
%read omega, alpha and bank from input file:
ls=0;
fid = fopen( strcat('./',outpath,'/',num2str(outpath),'.txt'), 'r' );
maxtime = step;

ls=0;
while strcmp(ls,'deltime')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
deltime = fscanf(fid,'%s',1);
deltime = str2num(deltime);

ls=0;
while strcmp(ls,'Uinf')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
uinf = fscanf(fid,'%s',1);
uinf = str2num(uinf);

while strcmp(ls,'alpha=')==0
    ls = fscanf(fid,'%s',1);
end
alphad = fscanf(fid,'%d',1);

ls=0;
while strcmp(ls,'circling')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
circling = fscanf(fid,'%d',1);

ls=0;
while strcmp(ls,'horizontal')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
horiz = fscanf(fid,'%d',1);

ls=0;
while strcmp(ls,'phi')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
bankd = fscanf(fid,'%d',1);

ls=0;
while strcmp(ls,'gradient')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
grad = fscanf(fid,'%s',1);
grad = str2num(grad);

ls=0;
while strcmp(ls,'cg')==0
    ls = fscanf(fid,'%s',1);
end
ls=0;
while strcmp(ls,'=')==0
    ls = fscanf(fid,'%s',1);
end
accg(1) = str2num(fscanf(fid,'%s',1));
accg(2) = str2num(fscanf(fid,'%s',1));
accg(3) = str2num(fscanf(fid,'%s',1));
fclose(fid);


if grad == 0
    grad = (9.81*tand(bankd)) / (uinf);
end


if circling == 0
    grad =0;
    horiz = 0;
    bankd = 0;
end


if horiz == 1
    alphad= 0;
end


% omega = 1.5456000000000001;
omega = grad*deltime*(maxtime+1);

subplot(1,2,2)

wingcgavg= [totalforcemom.x1(1) totalforcemom.y1(1) totalforcemom.z1(1)];
% r1= -(repmat(wingcgavg(:,1),[1,4])-repmat(wingcg(:,1),[1,4]));
% r2=-(repmat(wingcgavg(:,2),[1,4])-repmat(wingcg(:,2),[1,4]))
% r3=-(repmat(wingcgavg(:,3),[1,4])-repmat(wingcg(:,3),[1,4]))

%% rotate wing
for count = 1:4
    if count ==1
        tempx= wing(:,1:4)-repmat(wingcgavg(:,1),[1,4]);
        tempy= wing(:,5:8)-repmat(wingcgavg(:,2),[1,4]);
        tempz= wing(:,9:12)-repmat(wingcgavg(:,3),[1,4]);
        
    elseif count ==2
        tempx= elecg(:,1)-(wingcgavg(:,1));
        tempy= elecg(:,2)-(wingcgavg(:,2));
        tempz= elecg(:,3)-(wingcgavg(:,3));
                
    elseif count ==3
        tempx= rotforcemom.x1-(wingcgavg(:,1));
        tempy= rotforcemom.y1-(wingcgavg(:,2));
        tempz= rotforcemom.z1-(wingcgavg(:,3));
    elseif count ==4  
        tempx= uvel.u1;  
        tempy= uvel.v1; 
        tempz= uvel.w1;
    end
    rotx(:,:) = tempx .* cos(omega) + tempy .* sin(omega);
    roty(:,:) = -tempx .* sin(omega) + tempy .* cos(omega);
    rotz(:,:) = tempz;
    %
    
    tempx= rotx;
    tempy= roty;
    tempz= rotz;
    
    rotx(:,:) = tempx;
    roty(:,:) = tempy .* cosd(bankd) + tempz .* sind(bankd);
    rotz(:,:) = -tempy .* sind(bankd) + tempz .* cosd(bankd);
    %
    
    tempx= rotx;
    tempy= roty;
    tempz= rotz;
    
    rotx(:,:) = tempx .* cosd(alphad) + tempz .* sind(alphad);
    roty(:,:) = tempy;
    rotz(:,:) = -tempx .* sind(alphad) + tempz .* cosd(alphad);
    
    tempx= rotx;
    tempy= roty;
    tempz= rotz;
    
    clear rotx roty rotz
    if count ==1
        wingx(:,:) = tempx;
        wingy(:,:) = tempy;
        wingz(:,:) = tempz;
    elseif count ==2
        newwingcg(:,1)= tempx;
        newwingcg(:,2) = tempy;
        newwingcg(:,3) = tempz;
    elseif count ==3
        forcecent(:,1) = tempx;
        forcecent(:,2) = tempy;
        forcecent(:,3) = tempz;
    elseif count ==4
        uvelrot(:,1) = tempx;
        uvelrot(:,2) = tempy;
        uvelrot(:,3) = tempz;
    end
end
%% rotate all vectors:

tempx= direct.u1(:);
tempy= direct.v1(:);
tempz= direct.w1(:);

forcex(:,:) = tempx .* cos(omega) + tempy .* sin(omega);
forcey(:,:) = -tempx .* sin(omega) + tempy .* cos(omega);
forcez(:,:) = tempz;

tempx= forcex;
tempy= forcey;
tempz= forcez;

forcex(:,:) = tempx;
forcey(:,:) = tempy .* cosd(bankd) + tempz .* sind(bankd);
forcez(:,:) = -tempy .* sind(bankd) + tempz .* cosd(bankd);

tempx= forcex;
tempy= forcey;
tempz= forcez;
forcex(:,:) = tempx .* cosd(alphad) + tempz .* sind(alphad);
forcey(:,:) = tempy;
forcez(:,:) = -tempx .* sind(alphad) + tempz .* cosd(alphad);



hold on

% newwingcg =  ([mean(wingx,2),mean(wingy,2),mean(wingz,2)]);
patch(wingx(:,:)',wingy(:,:)',wingz(:,:)',[0.4 0.4 0.4],'FaceAlpha',0.8,'EdgeColor','k');


quiver3(forcecent(:,1),forcecent(:,2),forcecent(:,3),rotforcemom.u1(:),rotforcemom.v1(:),rotforcemom.w1(:));


temp= ([repelem(newwingcg,3,1);unique(forcecent, 'rows', 'stable')]);

try
    quiver3(temp(:,1),temp(:,2),temp(:,3),forcex,forcey,forcez);
catch
end
quiver3(newwingcg(:,1),newwingcg(:,2),newwingcg(:,3),uvelrot(:,1),uvelrot(:,2),uvelrot(:,3));
xlabel('X');
ylabel('Y');
span = max(max(newwingcg(:,:))) - min(min(newwingcg(:,:)));
% accg = repmat([mean(newwingcg(:,1)),mean(newwingcg(:,2)),mean(newwingcg(:,3))],[2,1])

V = [totalforcemom.u1,totalforcemom.v1,totalforcemom.w1];
V = V/norm(V) * span;
quiver3(0,0,0,V(1,1),V(1,2),V(1,3));
quiver3(0,0,0,V(2,1),V(2,2),V(2,3));
% zlim([-1 1]);
grid on
box on
axis tight
axis equal
ax = gca;
z = zoom;
setAxes3DPanAndZoomStyle(z,ax,'camera')
fclose all

