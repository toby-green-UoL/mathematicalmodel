%%%%%%%%%%%%%%%%
% light source %
%%%%%%%%%%%%%%%%

% select the light source
light_source_user_defined = true;
light_source_sun_moving = false;
light_source_dataset = false;
light_source_solar_simulator = false;

% parameters to set when light source is user_defined
if light_source_user_defined == true
    direct_normal_irradiance = 800; % [W/m^2]
    diffuse_irradiance = 300; % [W/m^2]
    zenith_degrees = 0; % the angle between the sun and the vertical
    zenith_radians = zenith_degrees * pi/180; %convert degrees to radians
    azimuth_degrees = 270; % clockwise from north
    azimuth_radians = azimuth_degrees * pi/180; %convert degrees to radians
    time_hours_start = 10; % this is for running the model
    time_hours_stop = 16; % this is for running the model
end

% parameters to set when light source is sun_moving
if light_source_sun_moving == true
    height_above_sea_level = 1200; % [m]
    latitude_of_site_deg = 0.345; % latitude of Uganda
    day_number = 1; % 1 Jan is day 1
    time_hours_start = 11;
    time_hours_stop = 13;
end

% parameters to set when light source is dataset
if light_source_dataset == true
    FID = fopen('LOG06.txt'); % filename
    marker = 9; % this is the position of the && marker. It is usually 9.
    clear_sky_threshold = 0.8; % factor for clear sky threshold
    
    % the following parameters are used to calculate the clear sky
    % global irradiance
    height_above_sea_level = 1200; % [m]
    latitude_of_site_deg = 0.345; % latitude of Uganda
    day_number = 1; % 1 Jan is day 1
    time_hours_start = 11;
    time_hours_stop = 13;
    
    irradiance_calibration = 1.5; % [W m^-2] of each point 
    time_calibration = 0.1; % timestep of each sample
end

% parameters to set when light source is solar_simulator
if light_source_solar_simulator == true
    
    time_hours_start = 11; % this is for running the model
    time_hours_stop = 13; % this is for running the model
    
    direct_normal_irradiance = 226; % W/m^2
    diffuse_irradiance = 127; % W/m^2
    zenith_radians = 0;
    azimuth_radians = 0;
    
    % set lamp radius and centres for the solar simulator
    % all units are m
    lamps_radius = 0.15;
    lamp_1_centre = [0.35,0.5, 1.2]; 
    lamp_2_centre = [0.65,0.5, 1.2]; 
    lamp_3_centre = [0.2,0.2, 1.2]; 
    lamp_4_centre = [0.5,0.2, 1.2]; 
    lamp_5_centre = [0.8,0.2, 1.2]; 
    lamp_6_centre = [0.2,0.8, 1.2]; 
    lamp_7_centre = [0.5,0.8, 1.2]; 
    lamp_8_centre = [0.8,0.8, 1.2];
end

zenith_max_degrees = 60; % horizon for diffuse irradiance

% intensity data with CSR = 10
% taken from Neumann et al, Transactions of the ASME, p198, v124, 2002.
intensity_data = [1,0.99,0.96,0.90,0.79,0.049,0.025,0.016,0.009,0.007,0.005,0,0,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry and properties of the collector %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

collector_linear_fresnel = true;
collector_compound_parabolic = false;
collector_parabolic_trough = false;

% parameters to set for Linear Fresnel
% note that mirrors are all the same size
% currently this is set for the lab demonstrator in Leeds
if collector_linear_fresnel == true
    collector_design_linear_fresnel = true;
    linear_fresnel_tracking = true;
    number_of_mirrors = 20;
    mirror_size_x = 0.1;
    mirror_size_y = 1.9;
    mirror_centre_y = 0;
    mirror_centre_z = 0.2;
    mirror_first_x = -1;
    mirror_spacing_x = 0.1;
end

mirror_reflectivity = 1;
mirror_soiling = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry and properties of the absorber and reactor %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Select the Absorber Type %%%
absorber_design_flat_panel = false;
absorber_design_pipe = true;

% parameters to set when the absorber is a flat panel
if absorber_design_flat_panel == true
    % all units are m
    absorber_panel_centre_point = [1.5,0.5,1.2];
    absorber_panel_size_x = 0.2;
    absorber_panel_size_y = 1.0;
    % tilt of the absorber normal, clockwise from vertical
    % these can be positive or negative values
    absorber_panel_tilt_x_deg = 40;
    absorber_panel_tilt_y_deg = 0; % keep at zero for linear designs
    absorber_panel_x_res = 0.02;
    absorber_panel_y_res = 0.02;
    element_thickness_panel = 0.0005;

    % Calculate the emissivity values
    % Thermal wavelengths refer to ambient and absorber temperatures
    % Optical wavelengths refer to solar temperatures
    % Emissivity values taken from: https://www.thermoworks.com/emissivity-table
    emissivity_thermal_panel = 0.65*ones(floor(absorber_panel_size_x/absorber_panel_x_res)+1,floor(absorber_panel_size_y/absorber_panel_y_res)+1); %oxidised copper
    emissivity_thermal_panel(2:9,2:49) = 0.95; %black paint
    
    emissivity_optical_panel = 0.65*ones(floor(absorber_panel_size_x/absorber_panel_x_res)+1,floor(absorber_panel_size_y/absorber_panel_y_res)+1); %oxidised copper
    emissivity_optical_panel(2:9,2:49) = 0.95; %black paint
    
    % Set Reactor Properties
    reactor_depth = 0.01;
    reactor_volume = absorber_panel_size_x*absorber_panel_size_y*reactor_depth;
    reactor_max_fill = 0.75; % fraction full of the reactor chambre
end

% parameters to set when the absorber is a pipe
if absorber_design_pipe == true

    % absorber pipe in the y-axis direction
    % this is set for the lab demonstrator
    absorber_pipe_x = 0; % centre point, not end point
    absorber_pipe_y = 0; % centre point, not end point
    absorber_pipe_z = 0.3; % centre point, not end point
    absorber_pipe_radius = 0.05; % radius
    absorber_pipe_length = 1.85; % length
    % for heat transfer code
    absorber_pipe_y_res = 0.05; % longitudinal resolution (m)
    absorber_pipe_theta_res = 10; % angular resolutin (degrees)
    % absorber pipe thickness
    element_thickness_pipe = 0.01;

    % Calculate the emissivity values
    % Thermal wavelengths refer to ambient and absorber temperatures
    % Optical wavelengths refer to solar temperatures
    % Emissivity values taken from: https://www.thermoworks.com/emissivity-table
    emissivity_thermal_pipe = 0.65*ones(absorber_pipe_length/absorber_pipe_y_res+1,360/absorber_pipe_theta_res); %oxidised copper
    emissivity_thermal_pipe(3:17,9:27) = 0.95; %black paint
    emissivity_optical_pipe = 0.65*ones(absorber_pipe_length/absorber_pipe_y_res+1,360/absorber_pipe_theta_res); %oxidised copper
    emissivity_optical_pipe(3:17,9:27) = 0.95; %black paint
end

%%% Water and Ambient Temperatures %%%
water_ambient_temperature = 292; % K
air_ambient_temperature = 295; % K

%%% convergence of lumped parameter iterations %%%
thermal_iterations = 100; % number of iterations
thermal_factor = 0.01; % decrease this for convergence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties of the biomass %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Continuous or Batch Process %%%

% the continuous process is for water heating only
% the interior of the absorber is divided into a series of elements
% in the y direction. Water shuttles along from one element to the next.
process_continuous = false;

% the batch process is for biomass processing only
% the interior of the absorber is treated as a single element
% so the assumption is effective mixing either by convection or mechanical
% mixing
process_batch = true;

if process_continuous == true
    water_flow_rate = 3.33E-6; % m^3 s^-1 flow rate of solar sim is 1L/5min
end

%%% Set Parameters for First Order Rate Equation(s) %%%

% data from Energies 2019, 12, 516
% Olive trimmings
%reactor_E = 42930; % activation energy J/mol
%reactor_A = 2510/60;

% data from: Luo Ge et al Kinetics of the pyrolytic and hydrothermal decomposition of water hyacinth. Bioresource Technology. 2011;102(13):6990-4.
% Water hyacinth
% reactor_E = 147110; % activation energy J/mol
% reactor_A = 1.74E16 / 3600; % pre-exponential factor /s
% n = 1;

% data from Energies 2019, 12, 516
% Olive trimmings
%reactor_E = 42930; % activation energy J/mol
%reactor_A = 2510/60;


% Data from (Lucian M, Volpe M & Fiori 2019) kinetic model consisting of 6 non-linbear differential equations  representing
%the carbon molar balance for different components 

% % Initial carbon molar concentrations, in mol / L
C_B = 38; % biomass
C_L = 0; % liquid intermediate
C_G1 = 0; % gas 1
C_HC1 = 0; % hydrochar
C_G2 = 0; % gas 2
C_HC2 = 0; % solid byproduct
% 
% %Pre-exponential factors 
% % olive trimmings
% k0_1 = 82.32/3600; %(s^-1)
% k0_2 = (2.51*10.^3)/3600;
% k0_3 = 7.99/3600;
% k0_4 = (5.38*10.^4)/3600;
% k0_5 = 1.41/3600;
% n = 1.5; % reaction order associated with k_5
% grape marc
%k0_1 = 41.55/3600;  %(s^-1)
%k0_2 = (1.69*10.^3)/3600;
%k0_3 = 11.28/3600;
%k0_4 = 0.0086/3600;
%k0_5 = 52.78/3600;
% Opuntia Ficus Indica %(s^-1)
%k0_1= 4.31/3600;
%k0_2 = 14.7/3600;
%k0_3= 9.37/3600;
%k0_4 = (2.4*10^4)/3600;
%k0_5 = (1.14*10^4)/3600;
% Water Hyacinth
k0_1 = 0; %(s^-1)
k0_2 = 0;
k0_3 = 1.74E16 / 3600;
k0_4 = 0;
k0_5 = 0;
n = 1; % reaction order associated with k_5
%activation energies (J/mol) 
%olive trimmings 
% Ea_1 = 22.03*(10.^3);
% Ea_2 = 42.93*(10.^3);
% Ea_3 = 7.75*(10.^3);
% Ea_4 = 67.35*(10.^3);
% Ea_5 = 10.37*(10.^3);
%grape marc (J/mol)
%Ea_1 = 20.23*10.^3;
%Ea_2 = 43.14*10.^3;
%Ea_3 = 9.18*10.^3;
%Ea_4= 4.11*10.^3;
%Ea_5 = 24.41*10.^3;
% Opuntia Ficus Indica (J/mol)
%Ea_1 = 9.82*10.^3;
%Ea_2 = 21.33*10.^3;
%Ea_3 = 8.39*10.^3;
%Ea_4 = 58.92*10.^3;
%Ea_5 = 7.44*10.^3;
%Water Hyacinth 
Ea_1 = 0;
Ea_2 = 0;
Ea_3 = 147.11*(10.^3);
Ea_4 = 0;
Ea_5 = 0;

%%%%%%%%%%%%%%%%%%%
% Time Parameters %
%%%%%%%%%%%%%%%%%%%

%%% Set Time Parameters %%%

% all time units are seconds
time_step_universal = 10; % this must be an integer fraction of the next 2 time steps
time_step_ray_tracing_calculation = 60; % 20 minutes
time_step_thermal_calculation = 30; % 2 minutes


%%%%%%%%%%%%%%%
% Ray Tracing %
%%%%%%%%%%%%%%%

%%% Set Number of Rays %%% 
 
% 200 is good for a figure
% 10000 is good for smooth absorber radiation patterns
number_of_rays_direct = 200;
number_of_rays_diffuse = 0 ;

% Do you wish to display the rays? 
display_rays = true;

%%% Maximum Zentith angle for diffuse radiation %%%

zenith_max_radians = 60 * (pi/180);

%%% Aperture and Viewing Area %%%

% Define an aperture through which all relevant rays will pass

rays_aperture_min_x = -1;
rays_aperture_max_x = 1.2;
rays_aperture_min_y = -0.95;
rays_aperture_max_y = 0.95;
rays_aperture_z = 0;

if light_source_solar_simulator == true
    rays_aperture_min_x = 0;
    rays_aperture_max_x = 1;
    rays_aperture_min_y = 0;
    rays_aperture_max_y = 1;
    rays_aperture_z = 0.2; % the plane of the mirrors
end


% Set the viewing area
floor_point = [0,0,0];
floor_normal = [0,0,1];
incident_ray_length = 1; % distance from aperture outwards


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material Properties and Constants %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set the Physical Constants %%%

% Stefan-Boltzmann Constant
SB = 5.67E-8;
% Natural Gas Constant 
R = 8.314; % J / mol K

%%% Set the Material Properties %%%
%cp_copper = 385; % J per kg per K
cp_reactor = 503;
cp_water = 4181; % liquid, J per kg per K
mix_water = 0.9; %for batch
% density_copper = 8960; % kg per m^3
density_reactor = 8050;
density_water = 997; % kg per m^3
% conductivity_copper = 38.4; % W / (m K) reduced by a factor of 10
conductivity_reactor = 8;
cp_biomass = 1455; %wh
density_biomass = 150; %  WH kg per m^3 
mix_biomass = 1 - mix_water; %for batch

% h_c_air_copper = 13.5; % convection to air, W /(m^2 K)
h_c_air_reactor = 11.3;                       %https://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
% h_c_water_copper = 400; % convection to water, W /(m^2 K)% Comparison of Heat Transfer Coefficients in Free and ForcedConvection using Circular Annular Finned Tubes, 
h_c_water_reactor = 380;                            %Dr. Abdul Jabbar, IJAIEM, 5,4, 2016 

                            
                            

