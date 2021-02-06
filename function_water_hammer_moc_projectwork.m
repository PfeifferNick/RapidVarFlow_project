clear;clc;close all
%todo
% friction is dependent on velocity and velocity on friciton?? 

%% set global variables to be used in other functions
global Ap a dt Tcl g f mode L D dx tend n

%% Choose between normalVersion/surgeTank/airChamber
%mode = "airChamber"; %% normalVersion/surgeTank/airChamber
mode = "surgeTank";
mode = "normalVersion";
%mode = "airChamber";

%% Choose between turbine and pumping mode 
% system = pump/turbine
system = "pump"
%system = "turbine"


%% Choose friction 
    friction = "normalFriction";
    friction = "transientFriction";

K_ut = 0.004; % 0.004 to 0.0054
K_ux = 0.033; % 0.033 to 0.05
theta = 0.5; 
fprintf('Mode is set to: %s, The friction is set to: %s.\n',mode,friction);
% else
% fprintf('Mode is set to: %s.\n',mode);   
% end

%% control volumes
%   choose number larger than 5. 
n  =  50 ;  %  > 5
%% input data
% Pipe data
H0 =  120 ;   % reservoir water level [m]
D  =    1.6 ;   % Pipe diameter [m]
Ap = pi*D^2/4 ; % Pipe area
L  = 1100.0 ;   % Pipe lenght [m]
L_1 = 886.0 ;
L_2 = 86.59 ;
L_3 = 127.41 ;

a  = 1000.0 ;   % wave velocity [m/s]
f  =  0.2; % friction coefficient
mu = 2*L/a  ;   % time for complete sequence of events
K=2.05*10^9; % water bulk modulus [N/m2]
rho=1000; % water density [kg/m3]
visc=1.004*10^(-6); % kinematic viscosity [m2/s]


if mode == "surgeTank"
    %%for surge tank in the middle 
    % Position of the surge tank
    nsT = 15%round(n/2)+1; % if n=5, nsT=4
    AsT = 20 ;% AsT Area of surge tank 20 m^2
    maxit = 1000; %maximum iteration steps
    tol = 1e-5; %tolerance
end


if mode == "airChamber"
    %%for airChamber in the middle 
    % Position of the air chamber
    nsT = 6%round(n/2)+1; % if n=5, nsT=4 %needs to be close to reservoir 
    Vol_air =2300; % volumen of air chamber [m^3]
    A_air = 4; % Area of air chamber [m^2]
    m_pol = 1.2; %polytrophic constant , assumption: air
    H_atm = 10;%10 mWS 1bar
    maxit = 1000; %maximum iteration steps
    tol = 1e-5; %tolerance
end

% Valve data
Tcl= 1e9 ;        %time of closure of the valve
Tcl= 18 ;        %time of closure of the valve
Dv  = 0.4 ;   % Valve diameter [m]
Av = pi*Dv^2/4 ;       % Valve opening area [m?]
zeta = 1;
g = 9.81 ;      % graviational constant
Cd = 0.8 ;

% simulation data
dx = L/n ;      % lenght of control volumes
dt = dx/a;      %length/wave velocity, CFl=dt*a/dx
tend = 100 ;  % simulation duration
tsteps = round(tend/dt) ; % simulated time steps

V_initial = sqrt(2*g*H0/(f*L/D+zeta*(Ap/Av)^2+1)); % remember that this velocity comes from Bernouli's energy conservation equation.

%% Separation of the pipe into control volumes
n1 = round(L_1/dx);
n2 = round(L_2/dx);
n3 = n-n1-n2;
sum_n = n1+n2+n3;

if sum_n== n 
    disp('separation of control volumes worked')
else
    return
end

%% Pipe properties

e=0.012; % pipe thickness [m]
E=210*10^9; % Youngs modulus of the pipe [N/m2]
%u??

% calculation of pipe friction factor
Re=V_initial*D/visc; % reynolds number 

% rough=0.00015; % pipe roughness [m] (steel)
rough=0.0003; % pipe roughness [m] (reinforced concrete)
% rough=0.000115; % pipe roughness [m] (ductile iron)
% rough=0.00003; % pipe roughness [m] (GRP)

% Pipe 1 
rough1 = 0.00015; % pipe roughness [m] (steel)
f1=0.25/(log(rough1/(D*3.7)+5.74/Re^0.9))^2; % friction factor
% Pipe 2 
rough2 = 0.0003; % pipe roughness [m] (reinforced concrete)
f2=0.25/(log(rough2/(D*3.7)+5.74/Re^0.9))^2; % friction factor
% Pipe 3 
rough3 = 0.0003; % pipe roughness [m] (reinforced concrete)
f3=0.25/(log(rough3/(D*3.7)+5.74/Re^0.9))^2; % friction factor

f = zeros(n+1,1);


for i = 1:n+1
    if i <= n1 
    f(i,1) = f1;
    elseif i <= n1+n2 
    f(i,1) = f2;
    else  
    f(i,1) = f3;
    end
end

%f=0.25/(log(rough/(D*3.7)+5.74/Re^0.9))^2; % friction factor
%f =0
% phi according to the pipe supporting condition (see table 1 in paper:https://scielo.conicyt.cl/pdf/oyp/n20/art07.pdf
% Pipe Case 2 Pipe anchored against any axial movement
% u is the Poissons ratio see Table 3 from paper line 43 (Larock et al., 2000)

u = 0.30; % Steel accordign to roughness coefficient
phi = (1/(1+e/D))*(1-u^2+2*(e/D)*(1+u)*(1+e/D));
% wave speed calculation
alpha = sqrt((K/rho)/(1+phi*(D*K)/(E*e)));
fprintf('wave speed alpha = %g, choosen wave spped = %g\n',alpha,a);

%% Numerical Caluclations - MoC

% either pump or turbine system
if system == "pump"
[H,V] = function_PumpSimulation(visc,H0,rough,e,K,phi,E,friction,rho,K_ut,K_ux,theta,V_initial,tsteps,Ap);                                
elseif system == "turbine"
    
%initialization 
if mode == "normalVersion"
H = zeros(n+1,tsteps);
V = zeros(n+1,tsteps);
elseif mode == "surgeTank"
H = zeros(n+1,tsteps);
V = zeros(n+2,tsteps);
elseif mode == "airChamber"
H = zeros(n+1,tsteps);
V = zeros(n+2,tsteps);
else 
    disp('no mode')
    return   
end


%% set initial conditions
% loop over lenght of pipe
if mode == "normalVersion"
for i = 1:n+1
    if i==1
        H(i,1) = H0;
    else
   H(i,1) =  H0-(i-1)*dx*f(i)/D*V_initial^2/(2*g);    %H(i,1) =  % remember that the initial conditions come from steady
           %state and from fluid mechanics you know how to
           %calculate the piezometric head 
    end
   V(i,1) = V_initial;
end
 
elseif mode == "surgeTank"
  for i = 1:n+1
    if i==1
        H(i,1) = H0;
    else
   H(i,1) =  H0-(i-1)*dx*f(i)/D*V_initial^2/(2*g);    %H(i,1) =  % remember that the initial conditions come from steady
                    %state and from fluid mechanics you know how to calculate the
        %piezometric head 
    end
  end
  for i = 1:n+2
   V(i,1) = V_initial;
  end
  
   
elseif mode == "airChamber"
  for i = 1:n+1
    if i==1
        H(i,1) = H0;
    else
   H(i,1) =  H0-(i-1)*dx*f(i)/D*V_initial^2/(2*g);    %H(i,1) =  % remember that the initial conditions come from steady
                    %state and from fluid mechanics you know how to calculate the
        %piezometric head 
    end
  end
  for i = 1:n+2
   V(i,1) = V_initial;
  end
end

 
%% calculate H and V
% loop over time steps
if mode == "normalVersion"      
for j = 1:tsteps-1 
    
    % left boundary condition
    H(1,j+1) = H0; 
    V(1,j+1) = V(2,j)-g/a*(H(2,j)-H0)-f(1)*dt/(2*D)*V(2,j)*abs(V(2,j));
    % calculation of pressure wave
    % loop over lenght of pipe
    if friction == "normalFriction"
        for i = 2:n
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f(i)*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
            V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));
        end
    elseif friction == "transientFriction"
        for i = 2:n
            if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-((f(i)*dt*a)/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j)-V(i+1,j)+V(i+1,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
            
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j)+V(i+1,j)-V(i+1,j)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            
            else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f(i)*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j-1)-V(i+1,j)+V(i+1,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
        
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j-1)+V(i+1,j)-V(i+1,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            end
        end
    else
       disp("choose friction");
       return
    end
    
    % right boundary condition
    
    % the valve closing scheme has been already implemented in the function:
    % [V]=water_hammer_moc_boundary(type,timestep,initial fluid velocity)
    % This function only allows two types of closure: 10 for sinusoidal and 20 for linear. 
    % More on writing matlab functions in: https://de.mathworks.com/help/matlab/ref/function.html
    % More on calling functions in matlab: https://de.mathworks.com/help/matlab/learn_matlab/calling-functions.html
    
%    help=(j-1)*dt/Tcl*pi();
%     if help <= pi()
%         V(n+1,j+1) = V_initial*0.5*(cos(help)+1);
%     else
%         V(n+1,j+1) = 0;
%     end 
%     H(n+1,j+1) = H(n,j)+a/g*(V(n,j)-V(n+1,j+1))-(a/g)*((f(n+1)*dt)/(2*D))*V(n+1,j)*abs(V(n+1,j));
 
    
    %right b.c. walve closure: based on Torricelli’s Law
    t=j*dt;
    phi = (Cd*Av)/Ap; 
    if (1-t/Tcl) >= 0 
    T_closure = (1-t/Tcl)^0.5;
    else 
    T_closure = 0;
    end
    % is j correct?
    %1-t/Tcl
    if (1-t/Tcl) >= 0
    Fric=((f(n)/(2*D))*V(n,j)*abs(V(n,j))*dt);
    Coef_a=  (g/a)^2;
    Coef_b = 2*(g/a)*Fric-2*V(n,j)*(g/a)-2*(g/a)^2*H(n,j)-(phi*T_closure)^2*2*g;
    Coef_c = -2*Fric*V(n,j)+V(n,j)^2+Fric^2+2*(g/a)*H(n,j)*V(n,j)-2*Fric*(g/a)*H(n,j)+(g*H(n,j)/a)^2;
   mat(j,1)=Coef_a;
   mat(j,2)=Coef_b;
   mat(j,3)=Coef_c;
    H(n+1,j+1) = (-Coef_b-sqrt(Coef_b^2-4*Coef_a*Coef_c))/(2*Coef_a);
    
    
    V(n+1,j+1) = phi*T_closure*sqrt(2*g*H(n+1,j+1)); %H0 ?H(n,j+1) this is
    
    else 
    V(n+1,j+1) = 0;
    H(n+1,j+1) = H(n,j)+(a/g)*(V(n,j)-V(n+1,j+1))-(a/g)*((f(n)*dt)/(2*D))*V(n,j)*abs(V(n,j));
    end
end

 
elseif mode == "surgeTank"
    count = 1;
for j = 1:tsteps-1 
    
    % left boundary condition
    H(1,j+1) = H0; 
    V(1,j+1) = V(2,j)-g/a*(H(2,j)-H0)-f(1)*dt/(2*D)*V(2,j)*abs(V(2,j));
    
    if friction == "normalFriction"
    % calculation of pressure wave % loop over lenght of pipe
    for i = 2:nsT
        
        if i == nsT
        H(nsT,j+1) = H(nsT,j);%H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
        
        Hcheck = zeros (2,maxit);
        
            while true
        %Hcheck(1,count) = Hcheck(2,count-1);
        Hcheck(1,count) = H(nsT,j+1);
        V(nsT,j+1) = V(nsT-1,j)+g/a*(H(nsT-1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT-1,j)*abs(V(nsT-1,j))*dt);
        V(nsT+1,j+1)= V(nsT+2,j)-g/a*(H(nsT+1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT+2,j)*abs(V(nsT+2,j))*dt);
        H(nsT,j+1) = H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
                
        Hcheck(2,count) = H(nsT,j+1);
        
        diff = abs(Hcheck(1,count)-Hcheck(2,count));
        
          if (diff<=tol) || (count== maxit) 
            break;
          else
        count=count+1;
          end
       end   
              
        else
            
        H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f(i)*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
        V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));
        
        end
    end
    for i = nsT+1:n
        
        %Note V is changed because it has n+2 values and H only n+1
        H(i,j+1) = (H (i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-(f(i)*dt*a/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)));
        V(i+1,j+1) = (V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)));
        
    end
    elseif friction == "transientFriction"      
        
        % calculation of pressure wave % loop over lenght of pipe
    for i = 2:nsT
        
        if i == nsT
        H(nsT,j+1) = H(nsT,j);%H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
        
        Hcheck = zeros (2,maxit);
        
            while true
        %Hcheck(1,count) = Hcheck(2,count-1);
        Hcheck(1,count) = H(nsT,j+1);
        V(nsT,j+1) = V(nsT-1,j)+g/a*(H(nsT-1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT-1,j)*abs(V(nsT-1,j))*dt);
        V(nsT+1,j+1)= V(nsT+2,j)-g/a*(H(nsT+1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT+2,j)*abs(V(nsT+2,j))*dt);
        H(nsT,j+1) = H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
                
        Hcheck(2,count) = H(nsT,j+1);
        
        diff = abs(Hcheck(1,count)-Hcheck(2,count));
        
          if (diff<=tol) || (count== maxit) 
            break;
          else
        count=count+1;
          end
            end   
              
        else
            
                if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-((f(i)*dt*a)/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j)-V(i+1,j)+V(i+1,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
            
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j)+V(i+1,j)-V(i+1,j)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            
                else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f(i)*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j-1)-V(i+1,j)+V(i+1,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
        
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j-1)+V(i+1,j)-V(i+1,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
                end  
        end
    end
    for i = nsT+1:n
        %Note V is changed because it has n+2 values and H only n+1
                        if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-((f(i)*dt*a)/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i,j)-V(i,j)-V(i+2,j)+V(i+2,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i,j))*abs((V(i+1,j)-V(i,j))/dx)-sign(V(i+2,j))*abs((V(i+1,j)-V(i+2,j))/dx)));
            
            V(i+1,j+1) = ((V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i+1,j)+(1-theta)*(V(i,j)-V(i,j)+V(i+2,j)-V(i+2,j)))+((K_ux*a/dx))*...
                        (sign(V(i,j))*abs(V(i+1,j)-V(i,j))+sign(V(i+2,j))*abs(V(i+1,j)-V(i+2,j)))))*(1/(1+theta*K_ut));
            
                        else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-(f(i)*dt*a/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i,j)-V(i,j-1)-V(i+2,j)+V(i+2,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i,j))*abs((V(i+1,j)-V(i,j))/dx)-sign(V(i+2,j))*abs((V(i+1,j)-V(i+2,j))/dx)));
        
            V(i+1,j+1) = ((V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i+1,j)+(1-theta)*(V(i,j)-V(i,j-1)+V(i+2,j)-V(i+2,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i,j))*abs(V(i+1,j)-V(i,j))+sign(V(i+2,j))*abs(V(i+1,j)-V(i+2,j)))))*(1/(1+theta*K_ut));
                        end
    end
        
    else
       disp("choose friction");
       return
    end
        % right boundary condition
    
    % the valve closing scheme has been already implemented in the function:
    % [V]=water_hammer_moc_boundary(type,timestep,initial fluid velocity)
    % This function only allows two types of closure: 10 for sinusoidal and
    % 20 for linear. 
    % More on writing matlab functions in: https://de.mathworks.com/help/matlab/ref/function.html
    % More on calling functions in matlab: https://de.mathworks.com/help/matlab/learn_matlab/calling-functions.html
%    help=(j-1)*dt/Tcl*pi();
%     if help <= pi()
%         V(n+2,j+1) = V_initial*0.5*(cos(help)+1);
%     else
%         V(n+2,j+1) = 0;
%     end 

    %right b.c. walve closure: based on Torricelli’s Law
    t=j*dt;
    phi = (Cd*Av)/Ap; 
    if (1-t/Tcl) >= 0 
    T_closure = (1-t/Tcl)^0.5;
    else 
    T_closure = 0;
    end
    % is j correct?
    %1-t/Tcl
    if (1-t/Tcl) >= 0
    Fric=((f(n+1)/(2*D))*V(n,j)*abs(V(n,j))*dt);
    Coef_a=  (g/a)^2;
    Coef_b = 2*(g/a)*Fric-2*V(n+1,j)*(g/a)-2*(g/a)^2*H(n,j)-(phi*T_closure)^2*2*g;
    Coef_c = -2*Fric*V(n+1,j)+V(n+1,j)^2+Fric^2+2*(g/a)*H(n,j)*V(n+1,j)-2*Fric*(g/a)*H(n,j)+(g*H(n,j)/a)^2;
   mat(j,1)=Coef_a;
   mat(j,2)=Coef_b;
   mat(j,3)=Coef_c;
    H(n+1,j+1) = (-Coef_b-sqrt(Coef_b^2-4*Coef_a*Coef_c))/(2*Coef_a);
    1-t/Tcl
    
    V(n+2,j+1) = phi*T_closure*sqrt(2*g*H(n+1,j+1)); %H0 ?H(n,j+1) this is
    
    else 
    V(n+2,j+1) = 0;
    H(n+1,j+1) = H(n,j)+a/g*(V(n+1,j)-V(n+2,j+1))-(a/g)*((f(n+1)*dt)/(2*D))*V(n+1,j)*abs(V(n+1,j));
    end
end
  
elseif mode == "airChamber"
  count = 1;  
  H_air =  H_atm+H(nsT,1) ; %initial absolute pressure head at air chamber
for j = 1:tsteps-1 
    
    % left boundary condition
    H(1,j+1) = H0; 
    V(1,j+1) = V(2,j)-g/a*(H(2,j)-H0)-f(1)*dt/(2*D)*V(2,j)*abs(V(2,j));
   
  if friction == "normalFriction"
    % calculation of pressure wave % loop over lenght of pipe
    for i = 2:nsT
        if i == nsT
        H(nsT,j+1) = H(nsT,j);%H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
        
        Hcheck = zeros (2,maxit);
        
            while true
        %Hcheck(1,count) = Hcheck(2,count-1);
        Hcheck(1,count) = H(nsT,j+1);
        V(nsT,j+1) = V(nsT-1,j)+g/a*(H(nsT-1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT-1,j)*abs(V(nsT-1,j))*dt);
        V(nsT+1,j+1)= V(nsT+2,j)-g/a*(H(nsT+1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT+2,j)*abs(V(nsT+2,j))*dt);
        H(nsT,j+1) = H(nsT,j)+((A_air*dt*m_pol)/2)*((H(i,j)+H_atm)^(1+1/m_pol)/(H_air^(1/m_pol)*Vol_air))*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
                
        Hcheck(2,count) = H(nsT,j+1);
        
        diff = abs(Hcheck(1,count)-Hcheck(2,count));
        
          if (diff<=tol) || (count== maxit) 
            break;
          else
        count=count+1;
             end
            end   
              
        else         
        H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f(i)*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
        V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));       
        end
    end
    for i = nsT+1:n
        %Note V is changed because it has n+2 values and H only n+1
        H(i,j+1) = (H (i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-(f(i)*dt*a/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)));
        V(i+1,j+1) = (V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)));      
    end
    elseif friction == "transientFriction"           
              % calculation of pressure wave % loop over lenght of pipe
    for i = 2:nsT
        
        if i == nsT
        H(nsT,j+1) = H(nsT,j);%H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
        
        Hcheck = zeros (2,maxit);
        
            while true
        %Hcheck(1,count) = Hcheck(2,count-1);
        Hcheck(1,count) = H(nsT,j+1);
        V(nsT,j+1) = V(nsT-1,j)+g/a*(H(nsT-1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT-1,j)*abs(V(nsT-1,j))*dt);
        V(nsT+1,j+1)= V(nsT+2,j)-g/a*(H(nsT+1,j)-H(nsT,j+1))-(f(i)/(2*D))*(V(nsT+2,j)*abs(V(nsT+2,j))*dt);
        H(nsT,j+1) = H(nsT,j)+((A_air*dt*m_pol)/2)*((H(i,j)+H_atm)^(1+1/m_pol)/(H_air^(1/m_pol)*Vol_air))*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
                    
        Hcheck(2,count) = H(nsT,j+1);
        
        diff = abs(Hcheck(1,count)-Hcheck(2,count));
        
          if (diff<=tol) || (count== maxit) 
            break;
          else
        count=count+1;
          end
            end   
              
        else
            
                if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-((f(i)*dt*a)/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j)-V(i+1,j)+V(i+1,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
            
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j)+V(i+1,j)-V(i+1,j)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            
                else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f(i)*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j-1)-V(i+1,j)+V(i+1,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
        
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j-1)+V(i+1,j)-V(i+1,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
                end  
        end
    end
    for i = nsT+1:n
        %Note V is changed because it has n+2 values and H only n+1
                        if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-((f(i)*dt*a)/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i,j)-V(i,j)-V(i+2,j)+V(i+2,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i,j))*abs((V(i+1,j)-V(i,j))/dx)-sign(V(i+2,j))*abs((V(i+1,j)-V(i+2,j))/dx)));
            
            V(i+1,j+1) = ((V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i+1,j)+(1-theta)*(V(i,j)-V(i,j)+V(i+2,j)-V(i+2,j)))+((K_ux*a/dx))*...
                        (sign(V(i,j))*abs(V(i+1,j)-V(i,j))+sign(V(i+2,j))*abs(V(i+1,j)-V(i+2,j)))))*(1/(1+theta*K_ut));
            
                        else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-(f(i)*dt*a/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i,j)-V(i,j-1)-V(i+2,j)+V(i+2,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i,j))*abs((V(i+1,j)-V(i,j))/dx)-sign(V(i+2,j))*abs((V(i+1,j)-V(i+2,j))/dx)));
        
            V(i+1,j+1) = ((V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f(i)*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i+1,j)+(1-theta)*(V(i,j)-V(i,j-1)+V(i+2,j)-V(i+2,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i,j))*abs(V(i+1,j)-V(i,j))+sign(V(i+2,j))*abs(V(i+1,j)-V(i+2,j)))))*(1/(1+theta*K_ut));
                        end
    end
    else
       disp("choose friction");
       return
    end
    % right boundary condition
    % the valve closing scheme has been already implemented in the function:
    % [V]=water_hammer_moc_boundary(type,timestep,initial fluid velocity)
    % This function only allows two types of closure: 10 for sinusoidal and
    % 20 for linear. 
    % More on writing matlab functions in: https://de.mathworks.com/help/matlab/ref/function.html
    % More on calling functions in matlab: https://de.mathworks.com/help/matlab/learn_matlab/calling-functions.html
    
   help=(j-1)*dt/Tcl*pi();
    if help <= pi()
        V(n+2,j+1) = V_initial*0.5*(cos(help)+1);
    else
        V(n+2,j+1) = 0;
    end 
     H(n+1,j+1) = H(n,j)+a/g*(V(n+1,j)-V(n+2,j+1))-(a/g)*((f(n+1)*dt)/(2*D))*V(n+1,j)*abs(V(n+1,j));   
end   

else 
    disp('no mode')
    return  
end


%CFL = max(V)*dt/dx;%%ask him again?
%printf(CFL)

%% calculate maximum data

% maximum pressure difference
maxH = max(max(H())) ;
output = sprintf('maximum pressure height at valve: %0.2f m', maxH) ;
disp(output)
minH = min(min(H())) ;
output = sprintf('minimum pressure height at valve: %0.2f m', minH) ;
disp(output)
excessH = maxH-H0 ;
output = sprintf('excess pressure height at valve: %0.2f m', excessH) ;
disp(output)

% comparison with maximum pressure equation
if mu/Tcl >= 1
    deltaH = a*V(n+1,1)/g ;
else
    deltaH = L/Tcl*V(n+1,1)/g ;
end
output = sprintf('excess pressure height at valve from equations in 4.2.1: %0.2f m'...
    , deltaH) ;
disp(output)


%% plot data



if mode == "normalVersion"
    
% set x-axis data for plots over the lenght of the pipe
xaxis = 0:dx:L;
xaxisV = xaxis;
xaxisH = xaxis;
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
            plot([0,100],[886,100],[20,972.59],[20,1100],'DisplayName','Pipe Elevation'), 
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxis,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxis,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end   
  
elseif mode == "surgeTank"

% set x-axis data for plots over the lenght of the pipe
xaxisold = 0:dx:L;
xaxisV = zeros(1,n+2) ;
xaxisH = zeros(1,n+1) ;
for i=1:nsT
    xaxisV(i) = xaxisold(i);
    xaxisH(i) = xaxisold(i);
end
for i=nsT+1:n+2
    xaxisV(i) = xaxisold(i-1);
end
for i=nsT+1:n+1
    xaxisH(i) = xaxisold(i-1);
end
% set x-axis data for plots over time
t = 0 : dt : (tsteps-1)*dt ;

for i=1:n+1
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
  plot(t,H(1,:),'b'), xlabel('t [s]'), ylabel('H [m]')
  hold on
  plot(t,H(n+1,:),'r'), 
  hold on
  plot(t,H(nsT-1,:)), 
  hold on
  plot(t,H(nsT+1,:)), 
  hold on
  plot(t,H(nsT,:)), 
  
  legend('at the inlet', 'at the valve','before the surge tank','after the surge tank','at the surge tank', 'Location', 'best')
  
subplot(2,2,2)
    hold off
    plot(t,V(1,:),'b'), xlabel('t [s]'), ylabel('V [m/s]')
    hold all
    plot(t,V(n+2,:),'r'), 
    hold on
    plot(t,V(nsT-1,:)), 
    hold on
    plot(t,V(nsT+1,:)), 
    hold on
    plot(t,V(nsT,:)), 
    legend('at the inlet', 'at the valve','before the surge tank','after the surge tank','at the surge tank', 'Location', 'best')
    
    % plot H at different time steps over the lenght of the pipe
    % reset hold for new plots
  for i =1:(1) % changing the number in parenthesis allows to animate the simulation from the first to the desired timestep
        subplot(2,2,3)       
            plot(xaxisH,H(:,i )), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)]);
        subplot(2,2,4)
            plot(xaxisV,V(:,i )), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
        pause(.0002) 
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
        subplot(2,2,3)
            plot(xaxisH,H(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)])
            hold on
        subplot(2,2,4)
            plot(xaxisV,V(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
            hold on
        subplot(2,2,3)    
            plot(xaxisH,H_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            plot(xaxisH,H_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxisV,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxisV,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end 
  
elseif mode == "airChamber"
    % set x-axis data for plots over the lenght of the pipe
xaxisold = 0:dx:L;
xaxisV = zeros(1,n+2) ;
xaxisH = zeros(1,n+1) ;
for i=1:nsT
    xaxisV(i) = xaxisold(i);
    xaxisH(i) = xaxisold(i);
end
for i=nsT+1:n+2
    xaxisV(i) = xaxisold(i-1);
end
for i=nsT+1:n+1
    xaxisH(i) = xaxisold(i-1);
end
% set x-axis data for plots over time
t = 0 : dt : (tsteps-1)*dt ;

for i=1:n+1
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
  plot(t,H(n+1,:)), 
  hold on
  plot(t,H(nsT-1,:)), 
  hold on
  plot(t,H(nsT+1,:)), 
  hold on
  plot(t,H(nsT,:)), 
  
  legend('at the inlet', 'at the valve','before the surge tank','after the surge tank','at the surge tank', 'Location', 'best')
  
subplot(2,2,2)
    hold off
    plot(t,V(1,:)), xlabel('t [s]'), ylabel('V [m/s]')
    hold all
    plot(t,V(n+2,:)), 
    hold on
    plot(t,V(nsT-1,:)), 
    hold on
    plot(t,V(nsT+2,:)), 
    hold on
    plot(t,V(nsT,:)), 
    
    legend('at the inlet', 'at the valve','before the surge tank','after the surge tank','at the surge tank', 'Location', 'best')
    
    % plot H at different time steps over the lenght of the pipe
    % reset hold for new plots
    for i =1:(1) % changing the number in parenthesis allows to animate the simulation from the first to the desired timestep
        subplot(2,2,3)       
            plot(xaxisH,H(:,i )), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)]);
        subplot(2,2,4)
            plot(xaxisV,V(:,i )), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
        pause(.0002) 
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
        subplot(2,2,3)
            plot(xaxisH,H(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)])
            hold on
        subplot(2,2,4)
            plot(xaxisV,V(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
            hold on
        subplot(2,2,3)    
            plot(xaxisH,H_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            plot(xaxisH,H_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxisV,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxisV,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end
else 
    disp('no mode')
    return   
end

 % plot H at different time steps over the lenght of the pipe
    % reset hold for new plots
    for i =1:(1) % changing the number in parenthesis allows to animate the simulation from the first to the desired timestep
        
        subplot(2,2,3)       
            plot(xaxisH,H(:,i )), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)]);
        subplot(2,2,4)
            plot(xaxisV,V(:,i )), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
        pause(.0002) 
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
        subplot(2,2,3)
            plot(xaxisH,H(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            axis([0,inf,(min(H_min_env)*1.1),(max(H_max_env)*1.1)])
            hold on
        subplot(2,2,4)
            plot(xaxisV,V(:,1),'DisplayName','steady state'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            axis([0,inf,(min(V_min_env)*1.1),(max(V_max_env)*1.1)]);
            hold on
        subplot(2,2,3)    
            plot(xaxisH,H_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            plot(xaxisH,H_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('H [m]')
            legend('Location','southoutside','orientation','horizontal')
        subplot(2,2,4)    
            plot(xaxisV,V_max_env,'DisplayName','max envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            plot(xaxisV,V_min_env,'DisplayName','min envelope'), xlabel ('Pipe Length [m]'),ylabel('V [m/s]')
            legend('Location','southoutside','orientation','horizontal')
    end
else 
   disp("choose system");
       return 
end % if system