function[H,V] = function_PumpSimulation(visc,H0,rough,e,K,phi,E,friction,rho,K_ut,K_ux,theta,V_initial,tsteps)

global Ap a dt Tcl g f mode L D dx tend n



% Pumping set up

Qo=.50; % required pumping discharge [m3/s]
USWL=20; % upstream water level [m]
DSWL=H0; % downstream water level [m]
PCL_US=100; % Pipe center line downstream [m]
PCL_DS=20; % Pipe center line upstream [m]

pipe_slope=(PCL_DS-PCL_US)/L; % used to plot the pipe (visual purposes)

% calculation of pipe friction factor

Vo=Qo/Ap;  % steady state flow velocity [m/s]
Re=Vo*D/visc; % reynolds number 
f=0.25/(log(rough/(D*3.7)+5.74/Re^0.9))^2; % friction factor

% head loss calculation

Hf=f*L/D*Vo^2/(2*g);
H_static=DSWL-USWL;  
TDH=Hf+H_static;

% phi according to the pipe supporting condition (see table 1 in paper:https://scielo.conicyt.cl/pdf/oyp/n20/art07.pdf
% Pipe Case 2 Pipe anchored against any axial movement
% u is the Poissons ratio see Table 3 from paper line 43 (Larock et al., 2000)

u = 0.30; % Steel accordign to roughness coefficient
phi = (1/(1+e/D))*(1-u^2+2*(e/D)*(1+u)*(1+e/D));
% wave speed calculation
alpha = sqrt((K/rho)/(1+phi*(D*K)/(E*e)));
fprintf('wave speed alpha = %g\n',alpha);

%% End of Student Input %%

%% Pump curve %%

% Pump curve characteristics
H1= TDH*1.20; % Shut off head [m] --> arbitrarly chosen
Q1= 0; % Shut off discharge [m3/s] --> zero by definition
H2= TDH; % operating point head [m] --> defined by operation point
Q2= Qo; % operating point discharge [m3/s] --> defined by operation point
H3= TDH*0.65; % run off head [m] --> arbitrarly chosen
Q3= Qo*1.75; % run off discharge [m3/s] --> arbitrarly chosen

Q_pump=[Q1 Q2 Q3];
H_pump=[H1 H2 H3];

% Calculation of pump curve equation coefficients

Coeff=polyfit(Q_pump,H_pump,2);

%% Student Input%%

% calculate coefficients usable in the characteristic equations (Equation
% 6.3.2 in the lecture notes). Note: to work with a coefficient of an
% array just use "name_of_the_array(n)" where n is the number of the
% element. In the case of "Coeff" elements are 1,2 and 3.

Ap= Coeff(1)*Ap^2; % not Ap'
Bp= Coeff(2)*Ap;
Cp= Coeff(3)+USWL;
%% End of Student Input %%

% Plot pump curve and system curve

Q=0:Q3/20:Q3;
for i=1:length(Q)
    H_fit(i)=Coeff(1)*Q(i)^2+Coeff(2)*Q(i)^1+Coeff(3);
    H_sys(i)=f*L/D*(Q(i)/Ap)^2/(2*g)+H_static;
end

figure(1)
title('Pump & System Curves')
plot(Q_pump,H_pump,'o')
hold on
plot(Q,H_fit)
plot(Q,H_sys)
hold off
xlabel('flow [m3/s]')
ylabel('Head [m]')
legend('pump curve points','pump curve','system curve')

%% Begin of Student Input %%
% domain discretization

tmstps=  round(tend/dt);% total number of timesteps

CFL = alpha*dt/dx;
fprintf('Courant number = %g\n',CFL);
% Initial conditions

for i=1:n+1
    H(i,1)= DSWL;  % the pressure head along the complete system equals the reservior level
    V(i,1)= 0; % system is shut off 
  
    pipe(i)=PCL_US+(i-1)*dx*pipe_slope; % Pipe elevation profile (this is just to show the pipe profile in the plot)
end

for j=1:tmstps-1
    % left boundary condition (constant speed pump)
    % coefficients of the quadratic formula. It is easier to have short
    % expressions to be able to chack any typing mistake.
    a=  (g/alpha)*Ap;
    b = (g/alpha)*Bp-1;
    c = (V(2,j)-(g/alpha)*H(2,j)+(g/alpha)*Cp-(f/(2*D))*V(2,j)*abs(V(2,j))*dt);
    
    
    V(1,j+1)=(-b-sqrt(b^2-4*a*c))/(2*a);
    H(1,j+1)=Ap*V(1,j+1)^2+Bp*V(1,j+1)+Cp;
  
    % Inner domain
    if friction == "normalFriction"
    for i=2:n
      H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(alpha/(2*g))*(V(i+1,j)-V(i-1,j))-(f*dt*alpha/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
      V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*alpha)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));
    end
    elseif friction == "transientFriction"
        for i = 2:n
            if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-((f*dt*a)/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j)-V(i+1,j)+V(i+1,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
            
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j)+V(i+1,j)-V(i+1,j)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            
            else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j-1)-V(i+1,j)+V(i+1,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
        
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j-1)+V(i+1,j)-V(i+1,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            end
        end
    else
       disp("choose friction");
       return
    end
    % Right boundary condition (constant level reservoir)
          H(n+1,j+1)= DSWL;
          V(n+1,j+1)= V(n,j)+(g/alpha)*(H(n,j)-H(n+1,j+1))-f*dt/(2*D)*V(n,j)*abs(V(n,j));
          
end

%% End of Student Input %%

% calculate maximum and minimum envelopes

for i=1:n+1
    H_max(i)=max(H(i,:));
    H_min(i)=min(H(i,:));
end

% define temporal axis

t=0:dt:dt*(tmstps-1);
x=0:dx:L;

pipex=[0,127.41,214,1100];
pipey=[20,20,120,120];
% define pipe profile

figure(2)
for i=1:20:tsteps
    plot(pipex,pipey)
    hold all
    plot(x,H_max)
    plot(x,H_min)
    plot(x,H(:,i))
    label_1=['Head at ',num2str((i-1)*dt), ' s'];
    legend('pipe profile','max head', 'min head',label_1, 'location', 'southeast')
    hold off
    %ylim([60 400])
    xlabel('pipe length [m]')
    ylabel('elevation/head [m]')
    pause(0.1)
end

figure(3)
title('Pump & System Curves')
plot(Q_pump,H_pump,'o')
hold on
plot(Q,H_fit)
plot(Q,H_sys)
hold off
xlabel('flow [m3/s]')
ylabel('Head [m]')
legend('pump curve points','pump curve','system curve')

  % plot H at valve over time
  figure(4);     % More information on ploting in: https://de.mathworks.com/help/matlab/ref/plot.html
  %set(gcf, 'Position', get(0, 'Screensize')); % if you dont like having a plot window that uses the whole screen just comment this line out.
  % suptitle('Valve closure applying MOC') The command suptitle is not supporte on the basic package
  subplot(2,2,1) 
  plot(t,H(1,:)), xlabel('t [s]'), ylabel('H [m]')
  hold on
  plot(t,H(n+1,:)), xlabel('t [s]'), ylabel('H [m]')
  legend('at the pump', 'at the reservoir', 'Location', 'best')
  
subplot(2,2,2)
    hold off
    plot(t,V(1,:)), xlabel('t [s]'), ylabel('V [m/s]')
    hold all
    plot(t,V(n+1,:)), xlabel('t [s]'), ylabel('V [m/s]')
    
end
