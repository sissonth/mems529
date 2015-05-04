%% true parameters
% DYNAMCIS
%omega= 60*(pi/180); % angular velocity % UNITS?
%omega=-5.2:0.01:5.2; % roughly -50 RPM to 50 RPM

omega=10*(pi/180)
V_dc= 5; %DC voltage
V_ac= 7; %AC voltage

% MATERIAL
E= 160e9; % youngs modulus of poly-silicon
density = 2329; %kg/m3

% GEOMETRY
% spring geometry 
t_b1= 30e-6; % beam 1 thickness
t_b2= 30e-6; % beam 2 thickness
w_b1= 2e-6; % beam 1 width
w_b2= 2e-6; % beam 2 width
l_b1= 125e-6; % beam 1 length 
l_b2= 125e-6; % beam 2 length

% outer frame geometry
outer_Lx=300e-6
outer_Ly=300e-6
outer_thick=t_b1;
outer_width=100e-6;

% innter frame geometry
inner_Lx=outer_Lx-2*outer_width-2*4*2*w_b1;   %100e-6
inner_Ly=0.9*(outer_Ly-2*outer_width-2*4*2*w_b2);   %100e-6 % the 0.9 is a safty factor for space
inner_thick=t_b1;
inner_width=30e-6;
inner_centerbeam_width=20e-6;

% comb geometry
% drive combs
N_drive = 200; % number of drive combs
length_drive=5e-6; % length of drive combs
width_drive=2e-6; % width of drive combs
thickness_drive = t_b1; % finger thickness
gap_drive = 2e-6;  % gap between comb fingers


% sense combs
N_sense = 15; % number of drive combs
length_sense=30e-6; %length of sense combs
width_sense=3e-6; %length of sense combs
thickness_sense = t_b1; % finger thickness
overlap_sense = 50e-6; % comb finger overlap
gap_sense = 2e-6;  % gap between comb fingers

h=20e-6; % gap between proof mass and subsrate 

% characteristic length
%Lc=1.9870e-06; %characterstic length
Lc=2e-6;

% pressure & temp.
%P=101325; % 1 atm ** LIKELY TO CHANGE ***
P=1000;
T=298; % temerature

% permitivity of vacuum/air
perm= 8.85e-12;

%% parameters
   M_outerframe=(2*outer_width*outer_Ly+2*outer_width*(outer_Lx-2*outer_width))*outer_thick*density;
   M_innerframe=(2*inner_width*inner_Ly+2*inner_width*(inner_Lx-2*inner_width)+(inner_Ly-2*inner_width)*inner_centerbeam_width)*inner_thick*density;
   M_m = M_outerframe+M_innerframe; % mass of inner and outer frames
   
   M_b = ((4*4*t_b1*w_b1*l_b1)+(4*4*t_b2*w_b2*l_b2))*density; % suspension beam mass
   
   M_drive=N_drive*width_drive*length_drive*thickness_drive*density;
   M_sense=N_sense*width_sense*length_sense*thickness_sense*density;
   M_c = M_drive+M_sense; % drive and sense combs mass
   
   M_e = 0; % etch holes mass
   
   M = M_m+M_b+M_c-M_e; % total mass
   
%M=4e-9; % mass in kg


%% stiffness coefficients
Kx=((4*E*t_b1*w_b1^3)/l_b1^3)+((4*E*t_b2*w_b2^3)/l_b2^3)
Ky=((4*E*t_b1*w_b1^3)/l_b1^3)+((4*E*t_b2*w_b2^3)/l_b2^3)
%Kx=30.95; % N/m
%Ky=30.95; % N/m

%% mode frequencies
% drive mode frequency
fx=(1/(2*pi))*sqrt(Kx/M);
% sense mode frequency
fy=(1/(2*pi))*sqrt(Ky/M);

%% knudsen number (Kn)

%%%%%%%% calculate lambda
%constants
R=8.3145; % gas constant
D=3.711e-10; % air molecule diameter
Na=6.022e23; %avagadros number, constant 

% mean free path calculation
lambda=(R*T)/(sqrt(2)*pi*D^2*Na*P);

% resulting kundsen number
Kn=lambda/Lc;

%% DAMPING coefficients
b=1.458e-5; %constant for air
S=110.4; %Kelvin
mu=(b*T^(3/2))/(T+S); % air viscosity

% areas
A_proof_drive=(2*inner_width*inner_Ly+2*inner_width*(inner_Lx-2*inner_width)+(inner_Ly-2*inner_width)*inner_centerbeam_width)+N_drive*width_drive*length_drive;
A_proof_sense=(2*outer_width*outer_Ly+2*outer_width*(outer_Lx-2*outer_width))+N_sense*width_sense*length_sense;

A_comb_drive=N_drive*thickness_drive*length_drive; % INCLUDES ALL COMBS
A_comb_sense=N_sense*thickness_sense*length_sense; % INCLUDES ALL COMBS


C_proof_sense=(mu*A_proof_sense)/((1+2*Kn)*h);
C_comb_sense=(mu*A_comb_sense)/(gap_sense+lambda);

C_proof_drive=(mu*A_proof_drive)/((1+2*Kn)*h);
C_comb_drive=(mu*A_comb_drive)/(gap_drive+lambda);

C_drive=C_proof_drive+C_comb_drive; % damping 
C_sense=C_proof_sense+C_comb_sense; % damping

%% electrostatic force from comb drive
%note that the 2*N is for the number of capcitors formed.
Fd=(2.28*2*N_drive*perm*thickness_drive*V_dc*V_ac)/gap_drive;

%% some maths via first order approximation methods

% drive static displacement
Xstatic=Fd/Kx;

% drive quality factor
Qdrive=sqrt(Kx*M)/C_drive;

% drive displacement
Xdrive=Xstatic*Qdrive;

%% comb drive velocity
Vdrive=Xdrive*2*pi*fx;

%% coriolis force
Fc=2*M*omega*Vdrive;

%% some other maths
Ystatic=Fc/Ky;
Qsense=sqrt(Ky*M)/C_sense; 
Ysense=Ystatic*Qsense; % sense displacement

%% change in capacitance
% change in capacitance
dCap=(perm*overlap_sense*thickness_sense*Ysense)/(gap_sense^2)

%% ac voltage reading
%Vdc_sense=5;
%current=Vdc_sense*dCap/fy;

% ASIC Voltage Output
Cref=0.15e-12; % reference capacitance
Vout=dCap/Cref


%% plotting

% figure
% plot(omega*180/pi,Vout*1000)
% grid on
% ylabel('Output Voltage (mV)');
% xlabel('Angular Velocity (deg/sec)')

omega*180/pi
Vout*1000

%%%%% NOTES
% CAN ADD TEMPERATURE PLOT WITH YOUNGS MODULUS EQUATION CHANGE
% CAN DO PLOT OF QUALITY FACTOR




