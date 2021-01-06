%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main script method of caracteristics
% project work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all
%% set global variables to be used in other functions
global Ap a dt dx Tcl g f mode n D V H

%% Choose between normalVersion/surgeTank/airChamber
%mode = "airChamber"; %% normalVersion/surgeTank/airChamber
%mode = "surgeTank";
mode = "normalVersion";
K_ut = 0.004; % 0.004 to 0.0054
K_ux = 0.033; % 0.033 to 0.05
theta = 0.5; 
%% Choose friction 
if mode == "normalVersion" 
    %friction = "normalFriction";
    friction = "transientFriction";
               
    comparisonFriction = "no";          %comparison of the transient with normal friction
    %comparisonFriction = "yes";
K_ut = 0.004; % 0.004 to 0.0054
K_ux = 0.033; % 0.033 to 0.05
theta = 0.5; 
fprintf('Mode is set to: %s, The friction is set to: %s.\n',mode,friction);
else
fprintf('Mode is set to: %s.\n',mode);   
end

%% input data

K=2.05*10^9; % water bulk modulus [N/m2]
g=9.81; % gravity acceleration [m/s2]
rho=1000; % water density [kg/m3]
visc=1.004*10^(-6); % kinematic viscosity [m2/s]

% Pipe data
H0 =  120 ;   % reservoir water level [m]
D  =    1.6 ;   % Pipe diameter [m]
Ap = pi*D^2/4 ; % Pipe area
a  = 1000.0 ;   % wave velocity [m/s]
n  =   25  ;   % control volumes
f  =  0.05 ; % friction coefficient
L  = 1100.0 ;   % Pipe lenght [m]
L_1 = 886.0 ;
L_2 = 86.59 ;
L_3 = 127.41 ;
mu = 2*L/a  ;   % time for complete sequence of events

% Pipe properties
rough=0.00015; % pipe roughness [m] (steel)
rough=0.0003; % pipe roughness [m] (reinforced concrete)
rough=0.000115; % pipe roughness [m] (ductile iron)
rough=0.00003; % pipe roughness [m] (GRP)
e=0.012; % pipe thickness [m]
E=210*10^9; % Young?s modulus of the pipe [N/m2]

 
% Position of the surge tank
    %%for surge tank in the middle 
    nsT = round(n/2)+1; % if n=5, nsT=4
    AsT = 4 ;% AsT Area of surge tank 20 m^2
    maxit = 1000; %maximum iteration steps
    tol = 1e-5; %tolerance

% Valve data
Tcl= 3 ;        %time of closure of the valve in seconds [s]
zeta = 1;
Dv = 0.4;       %Valve diameter [m]
Av = pi*Dv^2/4 ; % Valve area   [m^2]
Cd = 0.8;




% simulation data
dx = L/n ;      % lenght of control volumes
dt = dx/a;      %length/wave velocity, CFl=dt*a/dx
tend = 200 ;  % simulation duration
tsteps = round(tend/dt) ; % simulated time steps

V_initial = sqrt(2*g*H0/(f*L/D+zeta*(Ap/Av)^2+1)); % remember that this velocity comes from Bernouli's energy conservation equation.

%% Calculations 

[H,V] = intial_set_up_moc(tsteps,H0,L,zeta,Av,V_initial);

fprintf('initial set up finished\n')

[H,V] = running_moc(friction,tsteps,theta,K_ut,K_ux,H0,V_initial,nsT,maxit,AsT,tol,Cd,Av);

fprintf('simulation finished\n')

%CFL = max(V)*dt/dx;%%ask again?
%printf(CFL)

[maxH,minH] =  maximum_data_calc(mu,L,H0);



%% %%% plot data %%%
%%%%%%%%%%%%%%%%%%%%
if mode == "normalVersion"
    
% set x-axis data for plots over the lenght of the pipe
xaxis = 0:dx:L;
% set x-axis data for plots over time
t = 0 : dt : (tsteps-1)*dt ;

for i=1:n+1
    H_max_env(i) = max(H(i,:));
    H_min_env(i) = min(H(i,:));
end

for i=1:n+1
    V_max_env(i) = max(V(i,:));
    V_min_env(i) = min(V(i,:));
end
    if comparisonFriction == "yes"
        % plot comparison of friction
    figure(22);     
  
    subplot(1,2,1) 
    
    plot(t,H_normal(1,:),'r'), xlabel('t [s]'), ylabel('H [m]')
    hold on
    plot(t,H_normal(n+1,:),'g')
    hold on
    plot(t,H_transient(1,:),'b--')
    hold on
    plot(t,H_transient(n+1,:),'m--')
    hold off
    legend('inlet, normal', 'valve, normal','inlet, transient', 'valve, transient', 'Location', 'best')
  
    subplot(1,2,2)
    plot(t,V_normal(1,:),'r'), xlabel('t [s]'), ylabel('V [m/s]')
    hold on
    plot(t,V_normal(n+1,:),'g')
    hold on
    plot(t,V_transient(1,:),'b--')
    hold on
    plot(t,V_transient(n+1,:),'m--')
    hold off
   legend('inlet, normal', 'valve, normal','inlet, transient', 'valve, transient', 'Location', 'best')
    sgtitle(['dx = ',num2str(dx),'  ','n = ',num2str(n),'  ','dt =',num2str(dt),'  ','tend = ',num2str(tend),' K_u_t= ',num2str(K_ut),' K_u_x= ',num2str(K_ux)])
    saveas(gcf,['dx = ',num2str(dx),'  ','n = ',num2str(n),'  ','dt =',num2str(dt),'  ','tend = ',num2str(tend),' K_ut= ',num2str(K_ut),' K_ux= ',num2str(K_ux),'.bmp'])
    else

% plot H at valve over time
  figure();     % More information on ploting in: https://de.mathworks.com/help/matlab/ref/plot.html
  % suptitle('Valve closure applying MOC') The command suptitle is not supporte on the basic package
  subplot(2,2,1) 
  plot(t,H(1,:)), xlabel('t [s]'), ylabel('H [m]')
  hold on
  plot(t,H(n+1,:)), xlabel('t [s]'), ylabel('H [m]')
  legend('at the inlet', 'at the valve', 'Location', 'best')
  
subplot(2,2,2)
    hold off
    plot(t,V(1,:)), xlabel('t [s]'), ylabel('V [m/s]')
    hold all
    plot(t,V(n+1,:)), xlabel('t [s]'), ylabel('V [m/s]')
    legend('at the inlet', 'at the valve', 'Location', 'best')
    
    % plot H at different time steps over the lenght of the pipe
    % reset hold for new plots
    for i =1:(1) % changing the number in parenthesis allows to animate the simulation from the first to the desired timestep
        subplot(2,2,3)       
            plot(xaxis,H(:,i )), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)]);
        subplot(2,2,4)
            plot(xaxis,V(:,i )), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
        pause(.0002) 
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
        subplot(2,2,3)
            plot(xaxis,H(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)])
            hold on
        subplot(2,2,4)
            plot(xaxis,V(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
            hold on
        subplot(2,2,3)    
            plot(xaxis,H_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            plot(xaxis,H_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxis,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxis,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end   
    end  
elseif mode == "surgeTank"

% set x-axis data for plots over the lenght of the pipe
xaxisold = 0:dx:L;
xaxis = zeros(1,n+2) ;
for i=1:nsT
    xaxis(i) = xaxisold(i);
end
for i=nsT+1:n+2
    xaxis(i) = xaxisold(i-1);
end
% set x-axis data for plots over time
t = 0 : dt : (tsteps-1)*dt ;

for i=1:n+2
    H_max_env(i) = max(H(i,:));
    H_min_env(i) = min(H(i,:));
end

for i=1:n+2
    V_max_env(i) = max(V(i,:));
    V_min_env(i) = min(V(i,:));
end
% plot H at valve over time
  figure();     % More information on ploting in: https://de.mathworks.com/help/matlab/ref/plot.html
  % suptitle('Valve closure applying MOC') The command suptitle is not supporte on the basic package
  subplot(2,2,1) 
  plot(t,H(1,:)), xlabel('t [s]'), ylabel('H [m]')
  hold on
  plot(t,H(n+2,:)), xlabel('t [s]'), ylabel('H [m]')
  legend('at the inlet', 'at the valve', 'Location', 'best')
  
subplot(2,2,2)
    hold off
    plot(t,V(1,:)), xlabel('t [s]'), ylabel('V [m/s]')
    hold all
    plot(t,V(n+2,:)), xlabel('t [s]'), ylabel('V [m/s]')
    legend('at the inlet', 'at the valve', 'Location', 'best')
    
    % plot H at different time steps over the lenght of the pipe
    % reset hold for new plots
    for i =1:(1) % changing the number in parenthesis allows to animate the simulation from the first to the desired timestep
        subplot(2,2,3)       
            plot(xaxis,H(:,i )), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)]);
        subplot(2,2,4)
            plot(xaxis,V(:,i )), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
        pause(.0002) 
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
        subplot(2,2,3)
            plot(xaxis,H(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)])
            hold on
        subplot(2,2,4)
            plot(xaxis,V(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
            hold on
        subplot(2,2,3)    
            plot(xaxis,H_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            plot(xaxis,H_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxis,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxis,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end   
  
elseif mode == "airChamber"
    disp(' not implemented yet')
    return
else 
    disp('no mode')
    return   
end

 if comparisonFriction == "yes"
     
 else
 % plot H at different time steps over the lenght of the pipe
    % reset hold for new plots
    for i =1:(1) % changing the number in parenthesis allows to animate the simulation from the first to the desired timestep
        
        subplot(2,2,3)       
            plot(xaxis,H(:,i )), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)]);
        subplot(2,2,4)
            plot(xaxis,V(:,i )), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
        pause(.0002) 
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
        subplot(2,2,3)
            plot(xaxis,H(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)])
            hold on
        subplot(2,2,4)
            plot(xaxis,V(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
            hold on
        subplot(2,2,3)    
            plot(xaxis,H_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            plot(xaxis,H_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxis,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxis,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end    
 end


