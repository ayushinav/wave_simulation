% 1-D wave propagation using Finite Difference Method

%% Given
dx= 1;   % Each spatial step in mm
dt= 0.001;   % Each time step in sec
T= 0.4;     % Total time of simulation (sec)
X= 200;      % The length of the medium (mm)
nT= T/dt;   % Number of time steps
nX= X/dx;   % Number of spatial steps
U= zeros(nX, nT); % U(x,t)= wave equation parameter

f= 20;    % Frequency (Hz)
v= 1000;     % Velocity (mm/s)
c2= (v*dt/dx)^2;
%% Initial condition
% U(x,0)= sin(2*pi*f*t)
U(1,1:floor(v/f/dx))= 1*sin(2*pi*f*dt*(1:floor(v/f/dx)));
% du/dt |(x=0)= 0
U(2,1:floor(v/f/dx))= U(1,1:floor(v/f/dx));

%% Wave Propagation (FDM)

for j= 3:nT-1
    for i= 3:nX-1
        U1= 2*U(i,j-1)- U(i,j-2);
        U2= U(i+1,j-1)- 2*U(i,j-1)+ U(i-1,j-1);
        U(i,j)= U1+ c2*U2;
    end
    U(1,j)= 0;  % Boundary condition
    U(nX-1,j)= 0; % Boundary condition
end

%% Snaps
% times = [1 20 50 100]; 
% for i = 1:length(times)
%   figure(i);
%   k = times(i);
%   plot(1:nX,U(k,:),'linewidth',2);
%   grid on;
%   axis([1 nX -2 2]);
%   xlabel('X axis','fontSize',14);
%   ylabel('Wave Amplitude','fontSize',14);              
%   titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(k*dt), 'second'];
%   title(titlestring ,'fontsize',14);                              
%   
% end
%% Simulation
for j = 1:nT              
  plot((1:nX)*dx,U(:,j));
  grid on;
  %axis([1 nX -2 2]);
  %xlim([1 50])
  ylim([-2.5 2.5])
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(j), ' TIME = ',num2str(j*dt), 'second'];
  title(titlestring ,'fontsize',14);                            
  fh = figure(5);
  %set(fh, 'color', 'white'); 
  frame(j)=getframe(gcf);
           
end
%movie(F,T,1);

%% saving animation

vidObj = VideoWriter('teset.avi', 'Uncompressed AVI');
open(vidObj)
writeVideo(vidObj, frame)
close(vidObj)

