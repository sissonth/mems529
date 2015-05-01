%% true parameters
% DYNAMCIS
omega= 60*(pi/180); % angular velocity % UNITS?

% MATERIAL
E= 160e9; % youngs modulus of poly-silicon


% GEOMETRY
% spring geometry 
% % t_b1= ; % beam 1 thickness
% % t_b2= ; % beam 2 thickness
% % w_b1= ; % beam 1 width
% % w_b2= ; % beam 2 width
% % l_b1= ; % beam 1 length 
% % l_b2= ; % beam 2 length

% comb geometry
l_pc=50e-6; % comb finger overlap


%% parameters

% M_m= ; % mass of inner and outer frames
% M_b= ; % suspension beam mass
  M_c= ; % drive and sense combs mass
  M_e= 0; % etch holes mass
% M=M_m+M_b+M_c-M_e; % total mass
M=4e-9; % mass in kg


%% stiffness coefficients
% Kx=((4*E*t_b1*w_b1^3)/l_b1^3)+((4*E*t_b2*w_b2^3)/l_b2^3);
% Ky=((4*E*t_b1*w_b1^3)/l_b1^3)+((4*E*t_b2*w_b2^3)/l_b2^3);
Kx=30.95; % N/m
Ky=30.95; % N/m

%% mode frequencies
% drive mode frequency
fx=(1/(2*pi))*sqrt(Kx/M)
% sense mode frequency
fy=(1/(2*pi))*sqrt(Ky/M)

%% knudsen number (Kn)
Lc=1.9870e-06; %characterstic length

%%%%%%%% calculate lambda
% variables
P=101325; % 1 atm ** LIKELY TO CHANGE ***
%constants
R=8.3145; % gas constant
D=3.711e-10; % air molecule diameter
Na=6.022e23; %avagadros number, constant
T=298; % temerature 
% mean free path calculation
lambda=(R*T)/(sqrt(2)*pi*D^2*Na*P);

% resulting kundsen number
Kn=lambda/Lc;

%% DAMPING coefficients
b=1.458e-5; %constant for air
S=110.4; %Kelvin
T=298;
mu=(b*T^(3/2))/(T+S); % air viscosity

N= 500; % number of combs
%A_combs= ; % comb area
A_proof = 500e-6*500e-6; %proof mass area
g_c = 5e-6;% gap between comb finger
h=5e-6; % gap between proof mass and subsrate 

C_proof=(mu*A_proof)/((1+2*Kn)*h);
C_comb=(mu*N*A_combs)/(g_c+lambda);

C_drive=C_proof+C_comb; % damping 
C_sense=C_proof+C_comb; % damping

%% electrostatic force from comb drive
perm= 8.85e-12; %permitivity of air (vacuum, really). 
V_dc= 15; %DC voltage
V_ac= 10; %AC voltage
t_c= 10e-6; % finger thickness

Fd=(2.28*N*perm*t_c*V_dc*V_ac)/g_c;
%%
Fd/2.28/t_c*g_c/perm/500/15
%% some maths via first order approximation methods

% drive static displacement
Xstatic=Fd/Kx;

% drive quality factor
Qdrive=sqrt(Kx*M/C_drive);

% drive displacement
Xdrive=Xstatic*Qdrive;

%% comb drive velocity
Vdrive=Xdrive*2*pi*fx;

%% coriolis force
Fc=2*M*omega*Vdrive;

%% some other maths
Ystatic=Fc/Ky;
Qsense=sqrt(Ky*M/C_sense); 
Ysense=Ystatic*Qsense; % sense displacement

%% change in capacitance
% constants
perm= 8.85e-12; % permittivity of vacuum (air is close enough). 
% change in capacitance
dCap=(perm*l_pc*t_c*Ysense)/(g_c^2)




%% guessing some numbers and calculating dCap
% Ysense=6.75e-9; % meters
% t_c= 5e-6; % comb finger thickness
% perm= 8.85e-12; % permittivity of air
% l_pc= 50e-6; % comb finger overlap
% g_c= 2e-6; % gap between comb finger
% 
% dCap=(perm*l_pc*t_c*Ysense)/(g_c^2)
















