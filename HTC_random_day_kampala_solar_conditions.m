clear all;

%%%%%%%%%%%%%%%%
% Light source %
%%%%%%%%%%%%%%%%

%%% Select the Light Source %%%
light_source_user_defined = false;
light_source_sun_moving = false;
light_source_dataset = true;
light_source_solar_simulator = false;

%%%%%%%%%%%%%%%%
% Irradiance %
%%%%%%%%%%%%%%%%
    
%%% Global Horizontal Irradiance %%%

% GHI for fixed sun & moving sun
if (light_source_user_defined == true) || (light_source_sun_moving == true)
global_horizontal_irradiance = 1000; %user input
end

% GHI for data set
if light_source_dataset == true
     nchannels = 9;
     FID = fopen('july19.txt'); % filename
     marker = 9; % this is the position of the && marker. It is usually 9.
     
     raw_data = fread(FID,[2*(nchannels), inf],'uint8'); 
     ST = fclose(FID);
     [~,nsamples] = size(raw_data);
     measured_irradiance = zeros(8,nsamples);
     disp([num2str(nsamples) ' samples found, equivalent to ' num2str(nsamples/(10*60*60)) ' hours.']);

     day_number = 190;

     if(marker == 9) % channel 9
         channel = 1;
         measured_irradiance(1,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:); 
         channel = 2;
         measured_irradiance(2,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:); 
         channel = 3;
         measured_irradiance(3,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:);
         channel = 4;
         measured_irradiance(4,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:); 
         channel = 5;
         measured_irradiance(5,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:); 
         channel = 6;
         measured_irradiance(6,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:);
         channel = 7;
         measured_irradiance(7,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:); 
         channel = 8;
         measured_irradiance(8,:) = 256*raw_data(2*channel-1,:) + raw_data(2*channel,:); 
        disp('Data taken from logger channels 1,2,3,4,5,6,7,8.');
     end

     %%%%%%%%%%%%%%%%%%%%%%%%
     % plot irradiance data %
     %%%%%%%%%%%%%%%%%%%%%%%%

    time_x_axis = linspace(0,24,300001); % x axis
    single_day_irradiance_data = 1.075*measured_irradiance(2:4,600000:900000); % y axis
    figure(9);
    % Create the plot
    plot(time_x_axis,single_day_irradiance_data) 
    % Set axis limits
    xlim([0 24])
    ylim([0 1000])
    % Label the plot
    title('Solar irradiance in Kampala during a 24hr period in January 2019')
    xlabel('time (hours)')
    ylabel('Irradiance (W/m2)');

    real_day_GHI = max(single_day_irradiance_data,[],2);
    global_horizontal_irradiance = mean(real_day_GHI);    
    
     %%%%%%%%%%%%%%%%%%%%
     % Clear sky or not %
     %%%%%%%%%%%%%%%%%%%%

     % Convert into the approrite time step
     % Live data records 10 times a second
     % Find the mean irradiance across 1 minute
     % Find the mean value across 600 data points
     % Set the mean frequency to 600 if averaging over 1 minute.
    mean_frequency = 600; 
    mean_data_point = single_day_irradiance_data >= mean_frequency; 
     % Multiply by 1000 to find true value
    ture_mean_data_point = 1000*(mean_data_point);
     % Set the threshold irradiance value
     % The value that determines whether clear sky or not
    threshold_irradience = 500;
    % Create the new array
    size_new_array = size(ture_mean_data_point) ;
    New_data_array = arrayfun(@(k) mean(ture_mean_data_point(:,k:min(size_new_array(2),k+threshold_irradience-1)),2), 1:threshold_irradience:size_new_array(2), 'un', 0) ;
    New_data_array = [New_data_array{:}] ;
    
    new_time_x_axis = linspace(0,24,601);
    figure (10);
    plot(new_time_x_axis,New_data_array); 
    % Set axis limits
    xlim([0 24]);

    
    
end

%%% Set Locational Parameters %%%

% Set the elevation above sea level
height_above_sea_level = 1200; %user input for h
Fh1 = exp(-height_above_sea_level / 8000);
Fh2 = (-height_above_sea_level / 1250);

% Set the day of year
day_of_year = 190; %example day

% Set the latitude of the operation location
latitude_of_site_deg = 0.345; %latitude of Uganda
latitude_of_site_rad = pi*(latitude_of_site_deg/180);

% Calculate the zenith max
zenith_max = 60 * 2 * pi / 360; % horizon, measured from vertical

% State the intensity data
% intensity data with CSR = 10
% taken from Neumann et al, Transactions of the ASME, p198, v124, 2002.
intensity_data = [1,0.99,0.96,0.90,0.79,0.049,0.025,0.016,0.009,0.007,0.005,0,0,0,0];

%%% Calculate DNI and DI %%%

% Calculate direct and diffuse irradiance for a fixed sun solar source
if light_source_user_defined == true
    time_hours = 13; %user input
    hour_angle = 360 * (time_hours - 12) / 24;
    declination = -23.45 * cosd((day_of_year+10)*360/365);
    zenith = acosd(sind(declination)*sind(latitude_of_site_deg) + cosd(declination)*cosd(latitude_of_site_deg)*cosd(hour_angle));
    azimuth = acosd((sind(declination)*cosd(latitude_of_site_deg) - (cosd(declination)*sind(latitude_of_site_deg)*cosd(hour_angle)))/sind(zenith));
        if hour_angle > 0
           azimuth = 360 - azimuth;
        end
  % convert to radians
    zenith = zenith * 2 * pi / 360; % the angle between the sun and the vertical
    azimuth = azimuth * 2 * pi / 360; % clockwise from north
   

    %Air mass
    air_mass = 1 / cos(zenith);
                %M.J. Reno et al, Global Horizontal Irradiance Clear Sky Models: Implementation and Analysis, Sandia Report, SAND2012-2389, Sandia National Laboratories, USA, 2012 
    
    %extraterrestrial normal inicdent irradiance
    extraterrestrial_normal_incident_irradiace = 1367.7 * (1 + 0.033*cos(((2*pi)/365)* day_of_year));
                %M.J. Reno et al, Global Horizontal Irradiance Clear Sky Models: Implementation and Analysis, Sandia Report, SAND2012-2389, Sandia National Laboratories, USA, 2012

    %calcluate direct and diffuse irradiance
    direct_normal_irradiance = extraterrestrial_normal_incident_irradiace * 0.7^(air_mass^0.678);
    diffuse_irradiance = global_horizontal_irradiance - (direct_normal_irradiance * cos(zenith));
                %M.J. Reno et al, Global Horizontal Irradiance Clear Sky Models: Implementation and Analysis, Sandia Report, SAND2012-2389, Sandia National Laboratories, USA, 2012
    
end

% Calculate direct and diffuse irradiance for a moving sun solar source
% if light_source_sun_moving == true 
% see time stepping code within main code

% Calculate direct and diffuse irradiance for the data set
% to be added

% Calculate direct and diffuse irradiance for the solar simulator
if light_source_solar_simulator == true

    lamp_average_direct_mW = 40; %measured value, user input
    lamp_average_diffuse_mW = 22.5; %measured value, user input
    
    lamp_average_direct_irradiance = (lamp_average_direct_mW/1000)*(1/0.000177); %0.000177 is the area of measuring device. Converting mW into W/m2
    lamp_average_diffuse_irradiance = (lamp_average_diffuse_mW/1000)*(1/0.000177);
    
    direct_normal_irradiance = lamp_average_direct_irradiance; %W/m^2
    diffuse_irradiance = lamp_average_diffuse_irradiance; %W/m^2
    zenith = 0; %sun overhead
    azimuth = 0; %doesn't matter
    
    %Set lamp radius and centres for the solar simulator
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry and Properties of the Collector %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear Fresnel Collector Design
% note that mirrors are all the same size
% currently this is set for the lab demonstrator in Leeds
collector_design_linear_fresnel = true;
linear_fresnel_tracking = true;
number_of_mirrors = 45;
mirror_size_x = 0.14;
mirror_size_y = 10;
mirror_centre_y = 0;
mirror_centre_z = 0.2;
mirror_first_x = -3;
mirror_spacing_x = 0.14;

mirror_reflectivity = 1;
mirror_soiling = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry and Properties of the Absorber %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Select the Absorber Type %%%
absorber_design_flat_panel = false;
absorber_design_pipe = true;

%%% Flat Panel Absorber Properties %%%

% Set size and location
absorber_centre_point = [1.5,0.5,1.2];
absorber_size_x = 0.2;
absorber_size_y = 1.0;

% Set the tilt of the absorber normal, clockwise from vertical. These can be positive or negative values
absorber_tilt_x_deg = 0;
absorber_tilt_y_deg = 0;
absorber_panel_x_res = 0.02;
absorber_panel_y_res = 0.02;
element_thickness_panel = 0.0005;

% Calculate the emissivity values
% Thermal wavelengths refer to ambient and absorber temperatures
% Optical wavelengths refer to solar temperatures
% Emissivity values taken from: https://www.thermoworks.com/emissivity-table
emissivity_thermal_panel = 0.65*ones(floor(absorber_size_x/absorber_panel_x_res)+1,floor(absorber_size_y/absorber_panel_y_res)+1); %oxidised copper
emissivity_optical_panel = 0.95*ones(floor(absorber_size_x/absorber_panel_x_res)+1,floor(absorber_size_y/absorber_panel_y_res)+1); %black paint

% Set Reactor Properties
reactor_depth = 0.01;
reactor_volume = absorber_size_x*absorber_size_y*reactor_depth;
reactor_max_fill = 0.75; %you wouldn't fill a reactor chambre to the very top

%%% Pipe Absorber Properties %%%

% Absorber pipe is in the y-axis direction
% Currently set for the lab demonstrator

% Set size and location 
if absorber_design_pipe == true

    % absorber pipe in the y-axis direction
    % this is set for the lab demonstrator
    absorber_pipe_x = 0; % centre point, not end point
    absorber_pipe_y = 0; % centre point, not end point
    absorber_pipe_z = 0.3; % centre point, not end point
    absorber_pipe_radius = 0.07; % radius
    absorber_pipe_length = 10; % length
    % for heat transfer code
    absorber_pipe_y_res = 0.05; % longitudinal resolution (m)
    absorber_pipe_theta_res = 10; % angular resolutin (degrees)
    % absorber pipe thickness
    element_thickness_pipe = 0.0005;

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
%reactor_E = 147110; % activation energy J/mol
%reactor_A = 1.74E16 / 3600; % pre-exponential factor /s
%n = 1;

% data from Energies 2019, 12, 516
% Olive trimmings
%reactor_E = 42930; % activation energy J/mol
%reactor_A = 2510/60;


% Data from (Lucian M, Volpe M & Fiori 2019) kinetic model consisting of 6 non-linbear differential equations  representing
%the carbon molar balance for different components 

% Initial carbon molar concentrations, in mol / L
C_B = 40.174; % biomass
C_L = 0; % liquid intermediate
C_G1 = 0; % gas 1
C_HC1 = 0; % hydrochar
C_G2 = 0; % gas 2
C_HC2 = 0; % solid byproduct

%Pre-exponential factors 
% olive trimmings
k0_1 = 82.32/3600; %(s^-1)
k0_2 = (2.51*10.^3)/3600;
k0_3 = 7.99/3600;
k0_4 = (5.38*10.^4)/3600;
k0_5 = 1.41/3600;
n = 1.5; % reaction order associated with k_5
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
%activation energies (J/mol) 
%olive trimmings 
Ea_1 = 22.03*(10.^3);
Ea_2 = 42.93*(10.^3);
Ea_3 = 7.75*(10.^3);
Ea_4 = 67.35*(10.^3);
Ea_5 = 10.37*(10.^3);
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

%%%%%%%%%%%%%%%%%%%
% Time Parameters %
%%%%%%%%%%%%%%%%%%%

%%% Set Time Parameters %%%
time_hours_start = 10;
time_hours_stop = 17;

% all time units are seconds
time_step_universal = 20; % this must be an integer fraction of the next 2 time steps
time_step_ray_tracing_calculation = 120; % 20 minutes
time_step_thermal_calculation = 60; % 2 minutes

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

rays_aperture_min_x = -3;
rays_aperture_max_x = 3;
rays_aperture_min_y = -5;
rays_aperture_max_y = 5;
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
