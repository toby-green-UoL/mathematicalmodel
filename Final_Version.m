% This is an integrated model for solar powered biomass processing
% It is targeted at low-temperature processes
% Only linear focussing designs are accommodated
% Some code has been commented out for focussing to a point

% Code written by Toby Green and Rolf Crook
% We thank Valerie Dupont, Andy Ross, Peter Heggs and Darron Dixon Hardy for contibutions and advice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters from another file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
final_parameters;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initiate arrays and variables before time stepping loop %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The ray tracing part of the model associates a power (not energy) with
% each ray. The power is valid throughout the timestep.

% arrays used by mirrors
mirror_centre_point = zeros(number_of_mirrors,3);
mirror_point = zeros(number_of_mirrors,3);
mirror_corner_1 = zeros(number_of_mirrors,3);
mirror_corner_2 = zeros(number_of_mirrors,3);
mirror_corner_3 = zeros(number_of_mirrors,3);
mirror_corner_4 = zeros(number_of_mirrors,3);
%rot_mirror_corner_1 = zeros(number_of_mirrors,3);
%rot_mirror_corner_2 = zeros(number_of_mirrors,3);
%rot_mirror_corner_3 = zeros(number_of_mirrors,3);
%rot_mirror_corner_4 = zeros(number_of_mirrors,3);
%rot_mirror_max_x = zeros(number_of_mirrors,1);
%rot_mirror_min_x = zeros(number_of_mirrors,1);
%rot_mirror_max_y = zeros(number_of_mirrors,1);
%rot_mirror_min_y = zeros(number_of_mirrors,1);
mirror_max_x = zeros(number_of_mirrors,1);
mirror_min_x = zeros(number_of_mirrors,1);
mirror_max_y = zeros(number_of_mirrors,1);
mirror_min_y = zeros(number_of_mirrors,1);
mirror_tilt_x_deg = zeros(number_of_mirrors,1);
%mirror_tilt_y_deg = zeros(number_of_mirrors,1);
mirror_tilt_x = zeros(number_of_mirrors,1); % radians
%mirror_tilt_y = zeros(number_of_mirrors,1); % radians
mirror_normal = zeros(number_of_mirrors,3);
mirror_normal_unit = zeros(number_of_mirrors,3);
mirror_normal_direction = zeros(number_of_mirrors,1);

if absorber_design_flat_panel == true
    % power received by the absorber
    absorber_panel_array = zeros(floor(absorber_panel_size_x/absorber_panel_x_res)+1,floor(absorber_panel_size_y/absorber_panel_y_res)+1);
    
    element_area_panel = absorber_panel_x_res*absorber_panel_y_res;
    
    % heat capacities [J k^-1]
    c_water_reactor = absorber_panel_size_x*absorber_panel_size_y*reactor_depth*cp_water*density_water;
    c_element_panel = element_area_panel*element_thickness_panel*cp_copper*density_copper;
    c_water_element = absorber_panel_y_res*absorber_panel_size_x*reactor_depth*cp_water*density_water;

    element_conductivity_panel_x = conductivity_copper*element_thickness_panel*absorber_panel_y_res/absorber_panel_x_res; % in the x direction
    element_conductivity_panel_y = conductivity_copper*element_thickness_panel*absorber_panel_x_res/absorber_panel_y_res; % in the y direction
    element_conductivity_panel_d = conductivity_copper*absorber_panel_x_res*absorber_panel_y_res/element_thickness_panel; % in the depth direction

    absorber_temperature_panel = air_ambient_temperature*ones(floor(absorber_panel_size_x/absorber_panel_x_res)+1,floor(absorber_panel_size_y/absorber_panel_y_res)+1);
    absorber_temperature_panel_2 = zeros(floor(absorber_panel_size_x/absorber_panel_x_res)+1,floor(absorber_panel_size_y/absorber_panel_y_res)+1); % new temperatures after calculation
    
    water_temperature_panel = water_ambient_temperature*ones(absorber_panel_size_y/absorber_panel_y_res+1,1);

    
end

if absorber_design_pipe == true
    % power received by the absorber
    absorber_pipe_array = zeros(absorber_pipe_length/absorber_pipe_y_res+1,360/absorber_pipe_theta_res);
    
    element_area_pipe = absorber_pipe_y_res*2*3.141*absorber_pipe_radius*(absorber_pipe_theta_res/360);

    % heat capacities [J k^-1]
    c_water_reactor = absorber_pipe_length*3.141*(absorber_pipe_radius^2)*cp_water*density_water;
    c_element_pipe = element_area_pipe*element_thickness_pipe*cp_copper*density_copper;
    c_water_element = absorber_pipe_y_res*3.141*absorber_pipe_radius.^2*cp_water*density_water;

    element_conductivity_pipe_y = conductivity_copper*element_thickness_pipe*(2*3.141*absorber_pipe_radius*(absorber_pipe_theta_res/360))/absorber_pipe_y_res; % in the y direction
    element_conductivity_pipe_tan = conductivity_copper*element_thickness_pipe*absorber_pipe_y_res/(2*3.141*absorber_pipe_radius*(absorber_pipe_theta_res/360)); % in the tangential direction
    element_conductivity_pipe_rad = conductivity_copper*absorber_pipe_y_res*(2*3.141*absorber_pipe_radius*(absorber_pipe_theta_res/360))/element_thickness_pipe; % in the radial direction

    absorber_temperature_pipe = air_ambient_temperature*ones(absorber_pipe_length/absorber_pipe_y_res+1,360/absorber_pipe_theta_res);
    absorber_temperature_pipe_2 = zeros(absorber_pipe_length/absorber_pipe_y_res+1,360/absorber_pipe_theta_res); % new temperatures after calculation
    
    water_temperature_pipe = water_ambient_temperature*ones(absorber_pipe_length/absorber_pipe_y_res+1,1);
end

water_temperature_reactor = water_ambient_temperature;

%reaction_a = 0; % reaction progress (from 0 to 1)

if process_continuous == true
    time_water_shuttle_pipe = absorber_pipe_y_res*3.141*absorber_pipe_radius.^2 / water_flow_rate;
    time_left_water_shuttle = time_water_shuttle_pipe;
end

if process_batch == true
    C_total = C_B + C_L + C_G1 + C_HC2 + C_G2 + C_HC2;
end

time_seconds_start = 3600*time_hours_start;
time_seconds_stop = 3600*time_hours_stop;

log_index = 0; % index to the following logs
log_max_index = (time_seconds_stop-time_seconds_start)/time_step_thermal_calculation;
reactor_temperature_log = zeros(log_max_index,1);
time_seconds_log = zeros(log_max_index,1);
time_minutes_log = zeros(log_max_index,1);
time_hours_log = zeros(log_max_index,1);
reactor_power_in_log = zeros(log_max_index,1);
global_irradiance_log = zeros(log_max_index,1);
C_B_log = zeros(log_max_index,1);
C_B_dt_log = zeros(log_max_index,1);

C_L_log = zeros(log_max_index,1);
C_HC1_log = zeros(log_max_index,1);
C_HC2_log = zeros(log_max_index,1);
C_G1_log = zeros(log_max_index,1);
C_G2_log = zeros(log_max_index,1);
X_B_log = zeros(log_max_index,1);
X_L_log = zeros(log_max_index,1);
X_HC1_log = zeros(log_max_index,1);
X_HC2_log = zeros(log_max_index,1);
X_G1_log = zeros(log_max_index,1);
X_G2_log = zeros(log_max_index,1);
% reactor_a_log = zeros(log_max_index,1);
zenith_log = zeros(log_max_index,1);
azimuth_log = zeros(log_max_index,1);
aperture_incident_power_log = zeros(log_max_index,1);
absorber_incident_power_log = zeros(log_max_index,1);
%aperture_power_log = zeros(log_max_index,1);

%%%%%%%%%%%%%%%%%%%%
% RAY TRACING CODE %
%%%%%%%%%%%%%%%%%%%%

% calculate the geometry of the absorber, for a flat panel only
% tilt in both x and y directions is accommodated
% for linear designs, the tilt in the y direction should be zero
if absorber_design_flat_panel == true
 
    % calculate corners of the absorber
    % these will be used to determine if the intersect is on the absorber
    % convert tilts from degrees to radians
    absorber_panel_tilt_x_deg = absorber_panel_tilt_x_deg*(pi/180);
    absorber_panel_tilt_y_deg = absorber_panel_tilt_y_deg*(pi/180);
        
    % calculate normal vector of the plane of the absorber
    absorber_normal = [tan(absorber_panel_tilt_x_deg),tan(absorber_panel_tilt_y_deg),1];
    absorber_normal_unit = absorber_normal./norm(absorber_normal);
        
    % calculate the absorber normal angle as viewed from above
    [absorber_angle,~] = cart2pol(absorber_normal(1),absorber_normal(2));
        
    % calculate corners of the absorber from the normal and absorber_size
    U = cross([0,0,1],absorber_normal_unit);
    U_unit = U./norm(U);
    W = cross(absorber_normal_unit,U_unit);
    W_unit = W./norm(W);
    if(absorber_panel_tilt_x_deg == 0 && absorber_panel_tilt_y_deg == 0)
        U_unit = [0,1,0];
        W_unit = [1,0,0];
    end
    absorber_corner_1 = absorber_panel_centre_point+(absorber_panel_size_y/2)*U_unit+(absorber_panel_size_x/2)*W_unit;
    absorber_corner_2 = absorber_panel_centre_point+(absorber_panel_size_y/2)*U_unit-(absorber_panel_size_x/2)*W_unit;
    absorber_corner_3 = absorber_panel_centre_point-(absorber_panel_size_y/2)*U_unit-(absorber_panel_size_x/2)*W_unit;
    absorber_corner_4 = absorber_panel_centre_point-(absorber_panel_size_y/2)*U_unit+(absorber_panel_size_x/2)*W_unit;
        
    % rotate the absorber, as viewed from above, so the sides are parallel with x
    % and y axis
    absorber_rotation_matrix = [cos(absorber_angle) -sin(absorber_angle) 0; sin(absorber_angle) cos(absorber_angle) 0; 0 0 1];
        
    rot_absorber_corner_1 = (absorber_corner_1-absorber_panel_centre_point)*absorber_rotation_matrix+absorber_panel_centre_point;
    rot_absorber_corner_2 = (absorber_corner_2-absorber_panel_centre_point)*absorber_rotation_matrix+absorber_panel_centre_point;
    rot_absorber_corner_3 = (absorber_corner_3-absorber_panel_centre_point)*absorber_rotation_matrix+absorber_panel_centre_point;
    rot_absorber_corner_4 = (absorber_corner_4-absorber_panel_centre_point)*absorber_rotation_matrix+absorber_panel_centre_point;
        
    % calculate the mirror extent, to be used to determine if the interect lies
    % on the mirror
    rot_absorber_max_x = max([rot_absorber_corner_1(1),rot_absorber_corner_2(1),rot_absorber_corner_3(1),rot_absorber_corner_4(1)]);
    rot_absorber_min_x = min([rot_absorber_corner_1(1),rot_absorber_corner_2(1),rot_absorber_corner_3(1),rot_absorber_corner_4(1)]);
    rot_absorber_max_y = max([rot_absorber_corner_1(2),rot_absorber_corner_2(2),rot_absorber_corner_3(2),rot_absorber_corner_4(2)]);
    rot_absorber_min_y = min([rot_absorber_corner_1(2),rot_absorber_corner_2(2),rot_absorber_corner_3(2),rot_absorber_corner_4(2)]);
        
    % point on the plane of the absorber(x,y,z)
    absorber_point = absorber_panel_centre_point;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculations prior to ray tracing %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start of main timestepping loop
for time_seconds = time_seconds_start:time_step_universal:time_seconds_stop
    
    time_minutes = time_seconds /60;
    time_hours = time_minutes / 60;
    
    % see if ray tracing is required
    if(rem(time_seconds,time_step_ray_tracing_calculation) == 0)
        
        disp(['Ray tracing at time ',num2str(time_seconds),' s']);
        
        
        % recalculate azimuth and zenith angle
        % azimuth is degrees clockwise from north
        if (light_source_sun_moving == true || light_source_dataset == true)
            hour_angle = 360 * (time_hours - 12) / 24;
            declination = -23.45 * cosd((day_number+10)*360/365);
            zenith_degrees = acosd(sind(declination)*sind(latitude_of_site_deg) + cosd(declination)*cosd(latitude_of_site_deg)*cosd(hour_angle));
            azimuth_degrees = acosd((sind(declination)*cosd(latitude_of_site_deg) - (cosd(declination)*sind(latitude_of_site_deg)*cosd(hour_angle)))/sind(zenith_degrees));
            if hour_angle > 0
                azimuth_degrees = 360 - azimuth_degrees;
            end
            % convert to radians
            zenith_radians = zenith_degrees * 2 * pi / 360; % the angle between the sun and the vertical
            azimuth_radians = azimuth_degrees * 2 * pi / 360; % clockwise from north
            
            % the following equations are included in
            % M.J. Reno et al, Global Horizontal Irradiance Clear Sky Models: Implementation and Analysis, Sandia Report, SAND2012-2389, Sandia National Laboratories, USA, 2012 
            % this is all assuming clear sky conditions
            
            air_mass = 1 / cos(zenith_radians);
    
            extraterrestrial_normal_incident_irradiance = 1367.7 * (1 + 0.033*cos(((2*pi)/365)* day_number));

            % Meinel model
            direct_normal_irradiance = extraterrestrial_normal_incident_irradiance * 0.7^(air_mass^0.678);

            % an alternative calculation for DNI, the Daneshyar-Paltridge-Proctor (DPP) model
            % direct_normal_irradiance = 950.2*(1 - exp(0-0.075*(90 - zenith_degrees)));
            
            % the Daneshyar-Paltridge-Proctor (DPP) model
            diffuse_irradiance = 14.29 + 21.04 * ((pi/2) - zenith_radians); % diffuse horizontal irradiance
             
            global_horizontal_irradiance = diffuse_irradiance + (direct_normal_irradiance * cos(zenith_radians));
    
        end
        
        % calculate the centre points and tilts of the mirrors
        % this is for an array of flat mirrors
        for mirror = 1:number_of_mirrors
            % mirror centre point
            mirror_centre_point(mirror,:) = [mirror_first_x+mirror_spacing_x*mirror,mirror_centre_y,mirror_centre_z];
            
            % tilt of the mirror normal, clockwise from vertical
            % these can be positive or negative values
            if absorber_design_flat_panel == true
                delta_x = absorber_panel_centre_point(1) - mirror_centre_point(mirror,1);
                delta_z = absorber_panel_centre_point(3) - mirror_centre_point(mirror,3);
            end
            
            if absorber_design_pipe == true
                delta_x = absorber_pipe_x - mirror_centre_point(mirror,1);
                delta_z = absorber_pipe_z - mirror_centre_point(mirror,3);
            end
            
            mirror_tilt_x_deg(mirror) = 0.5*atand(delta_x/delta_z);
            % mirror_tilt_y_deg(mirror) = 0; % zero for linear designs
            
            if linear_fresnel_tracking == true
                solar_delta_r = tan(zenith_radians); % zenith is in radians
                [solar_delta_x,solar_delta_y] = pol2cart(0.5*pi - azimuth_radians,tan(zenith_radians));
                solar_delta_z = 1;
                
                theta_solar_x = atand(solar_delta_x/solar_delta_z);
                theta_absorber_x = atand(delta_x/delta_z);
                
                mirror_tilt_x_deg(mirror) = 0.5*(theta_absorber_x+theta_solar_x);
                % mirror_tilt_y_deg(mirror) = 0; % zero for linear designs
            end
            
        end
        
        % calculate the corner points of the mirror
        for mirror = 1:number_of_mirrors
            % calculate corners of the mirror
            % these will be used to determine if the intersect is on the mirror
            % convert tilts from degrees to radians
            mirror_tilt_x(mirror) = mirror_tilt_x_deg(mirror)*(pi/180);
            % mirror_tilt_y(mirror) = mirror_tilt_y_deg(mirror)*(pi/180);
            
            % calculate normal vector of the plane of the mirror
            % mirror_normal(mirror,:) = [tan(mirror_tilt_x(mirror)),tan(mirror_tilt_y(mirror)),1];
            mirror_normal(mirror,:) = [tan(mirror_tilt_x(mirror)),0,1];
            mirror_normal_unit(mirror,:) = mirror_normal(mirror,:)./norm(mirror_normal(mirror,:));
            
            % calculate the mirror normal angle as viewed from above
           % [mirror_normal_direction(mirror),~] = cart2pol(mirror_normal(mirror,1),mirror_normal(mirror,2));
            
            % calculate corners of the mirror from the normal and mirror_size
            U = cross([0,0,1],mirror_normal_unit(mirror,:));
            U_unit = U./norm(U);
            W = cross(mirror_normal_unit(mirror,:),U_unit);
            W_unit = W./norm(W);
            if(mirror_tilt_x_deg(mirror) == 0)
                U_unit = [0,1,0];
                W_unit = [1,0,0];
            end
            mirror_corner_1(mirror,:) = mirror_centre_point(mirror,:)+(mirror_size_y/2)*U_unit+(mirror_size_x/2)*W_unit;
            mirror_corner_2(mirror,:) = mirror_centre_point(mirror,:)+(mirror_size_y/2)*U_unit-(mirror_size_x/2)*W_unit;
            mirror_corner_3(mirror,:) = mirror_centre_point(mirror,:)-(mirror_size_y/2)*U_unit-(mirror_size_x/2)*W_unit;
            mirror_corner_4(mirror,:) = mirror_centre_point(mirror,:)-(mirror_size_y/2)*U_unit+(mirror_size_x/2)*W_unit;
            
            mirror_max_x(mirror) = max([mirror_corner_1(mirror,1),mirror_corner_2(mirror,1),mirror_corner_3(mirror,1),mirror_corner_4(mirror,1)]);
            mirror_min_x(mirror) = min([mirror_corner_1(mirror,1),mirror_corner_2(mirror,1),mirror_corner_3(mirror,1),mirror_corner_4(mirror,1)]);
            mirror_max_y(mirror) = max([mirror_corner_1(mirror,2),mirror_corner_2(mirror,2),mirror_corner_3(mirror,2),mirror_corner_4(mirror,2)]);
            mirror_min_y(mirror) = min([mirror_corner_1(mirror,2),mirror_corner_2(mirror,2),mirror_corner_3(mirror,2),mirror_corner_4(mirror,2)]);
            
            
            % the following commented code is for designs that focus to a point
            
            % rotate the mirror, as viewed from above, so the sides are parallel with x
            % and y axis
            %mirror_rotation_matrix = [cos(mirror_normal_direction(mirror)) -sin(mirror_normal_direction(mirror)) 0; sin(mirror_normal_direction(mirror)) cos(mirror_normal_direction(mirror)) 0; 0 0 1];
            
            %rot_mirror_corner_1(mirror,:) = (mirror_corner_1(mirror,:)-mirror_centre_point(mirror,:))*mirror_rotation_matrix+mirror_centre_point(mirror,:);
            %rot_mirror_corner_2(mirror,:) = (mirror_corner_2(mirror,:)-mirror_centre_point(mirror,:))*mirror_rotation_matrix+mirror_centre_point(mirror,:);
            %rot_mirror_corner_3(mirror,:) = (mirror_corner_3(mirror,:)-mirror_centre_point(mirror,:))*mirror_rotation_matrix+mirror_centre_point(mirror,:);
            %rot_mirror_corner_4(mirror,:) = (mirror_corner_4(mirror,:)-mirror_centre_point(mirror,:))*mirror_rotation_matrix+mirror_centre_point(mirror,:);
            
            % calculate the mirror extent, to be used to determine if the interect lies
            % on the mirror
            %rot_mirror_max_x(mirror) = max([rot_mirror_corner_1(mirror,1),rot_mirror_corner_2(mirror,1),rot_mirror_corner_3(mirror,1),rot_mirror_corner_4(mirror,1)]);
            %rot_mirror_min_x(mirror) = min([rot_mirror_corner_1(mirror,1),rot_mirror_corner_2(mirror,1),rot_mirror_corner_3(mirror,1),rot_mirror_corner_4(mirror,1)]);
            %rot_mirror_max_y(mirror) = max([rot_mirror_corner_1(mirror,2),rot_mirror_corner_2(mirror,2),rot_mirror_corner_3(mirror,2),rot_mirror_corner_4(mirror,2)]);
            %rot_mirror_min_y(mirror) = min([rot_mirror_corner_1(mirror,2),rot_mirror_corner_2(mirror,2),rot_mirror_corner_3(mirror,2),rot_mirror_corner_4(mirror,2)]);
            
            % point on the plane of the mirror(x,y,z)
            mirror_point(mirror,:) = mirror_centre_point(mirror,:);
        end
        
        
        
        % aperture and power through aperture calculations
        aperture_area = (rays_aperture_max_x - rays_aperture_min_x) * (rays_aperture_max_y - rays_aperture_min_y);
        incident_power_direct = direct_normal_irradiance * aperture_area * cos(zenith_radians) * mirror_reflectivity * mirror_soiling;
        incident_power_diffuse = diffuse_irradiance * aperture_area;
        
        
        %%% initiate vectors and matrices for the rays %%%
        number_of_rays_total = number_of_rays_direct + number_of_rays_diffuse;
        
        % point on each incident ray located on the aperture
        % there are many rays, so this is a 2 dimensional matrix
        % the first index indentifies the ray
        % the second index provides the spatial dimension (x,y,z)
        % using the monte-carlo method, so pick random points
        ray_aperture_point = zeros(number_of_rays_total,3);
        ray_incident_point = zeros(number_of_rays_total,3);
        
        ray_aperture_point(:,1) = rays_aperture_min_x + rand(number_of_rays_total,1)*(rays_aperture_max_x - rays_aperture_min_x); % x
        ray_aperture_point(:,2) = rays_aperture_min_y + rand(number_of_rays_total,1)*(rays_aperture_max_y - rays_aperture_min_y); % y
        ray_aperture_point(:,3) = rays_aperture_z; % z
        
        
        % direction of the incident ray (x,y,z)
        % azimuth (compass direction)of each ray will be random, and evenly distributed, in rad
        % the zenith is the angle from the vertical in mrad
        
        ray_incident_direction = zeros(number_of_rays_total,3);
        ray_intensity = zeros(number_of_rays_total,1);
        
        
        % calculate the direction vectors for the dni rays
        % these will all be unit vectors, away from the sun
        unscaled_power_direct = 0;
        for ray = 1:number_of_rays_direct
            % first assume the sun is directly overhead
            % work with cartesian coordinates, to get an even coverage
            ray_incident_direction(ray,1) = 2*(rand - 0.5)*sin(0.01); % +-10 mrad range
            ray_incident_direction(ray,2) = 2*(rand - 0.5)*sin(0.01); % +-10 mrad range
            ray_incident_direction(ray,3) = -1; % away from the sun
            
            if light_source_solar_simulator == true
                 ray_incident_direction(ray,1) = 2*(rand - 0.5)*sind(3.5); % +-10 mrad range
                 ray_incident_direction(ray,2) = 2*(rand - 0.5)*sind(3.5); % +-10 mrad range
                 ray_incident_direction(ray,3) = -1; % away from the sun
            end
            
            ray_incident_radius = sqrt(ray_incident_direction(ray,1)^2+ray_incident_direction(ray,2)^2);
            ray_incident_zenith = 1000*atan(ray_incident_radius); % in mrad
            
            % calculate the ray intensity
             if light_source_solar_simulator == false 
                ray_intensity(ray,1) = intensity_data(ceil(ray_incident_zenith));
             end
             
             
             if light_source_solar_simulator == true
                %find the distance between ray and solar simulator
                %lamp1
                ray_ss_distance_1 = sqrt(((ray_aperture_point(ray,1)-lamp_1_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_1_centre(2))^2));
                if ray_ss_distance_1 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
               
                %lamp2
                ray_ss_distance_2 = sqrt(((ray_aperture_point(ray,1)-lamp_2_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_2_centre(2))^2));
                if ray_ss_distance_2 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
               
                %lamp3
                ray_ss_distance_3 = sqrt(((ray_aperture_point(ray,1)-lamp_3_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_3_centre(2))^2));
                if ray_ss_distance_3 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
                
                %lamp4
                ray_ss_distance_4 = sqrt(((ray_aperture_point(ray,1)-lamp_4_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_4_centre(2))^2));
                if ray_ss_distance_4 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
                
                %lamp5
                ray_ss_distance_5 = sqrt(((ray_aperture_point(ray,1)-lamp_5_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_5_centre(2))^2));
                if ray_ss_distance_5 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
                
                %lamp6
                ray_ss_distance_6 = sqrt(((ray_aperture_point(ray,1)-lamp_6_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_6_centre(2))^2));
                if ray_ss_distance_6 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
                
                %lamp7
                ray_ss_distance_7 = sqrt(((ray_aperture_point(ray,1)-lamp_7_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_7_centre(2))^2));
                if ray_ss_distance_7 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
               
                %lamp8
                ray_ss_distance_8 = sqrt(((ray_aperture_point(ray,1)-lamp_8_centre(1))^2)+((ray_aperture_point(ray,2)-lamp_8_centre(2))^2));
                if ray_ss_distance_8 < lamps_radius
                    ray_intensity(ray,1) = 100;
                end
                        
             end
            
                if (light_source_sun_moving == true )|| (light_source_user_defined == true)||(light_source_dataset == true)
                % second rotate about x axis, equivalent to east-west
                rotate_x = [1 0 0; 0 cos(-zenith_radians) -sin(-zenith_radians); 0 sin(-zenith_radians) cos(-zenith_radians)];
                % rotate about z axis, clockwise
                rotate_z = [cos(-azimuth_radians) -sin(-azimuth_radians) 0; sin(-azimuth_radians) cos(-azimuth_radians) 0; 0 0 1];
                ray_incident_direction(ray,:) = (rotate_z * (rotate_x * ray_incident_direction(ray,:)'))';
                end
            
            %ray_incident_azimuth = 2*pi*rand; % in rad
            %ray_incident_zenith = 10*rand; % in mrad
            
            %ray_incident_azimuth = azimuth_deg * 2 * 3.141 / 360; % in rad
            %ray_incident_zenith = zenith_deg * 2 * 3.141 / 360; % in rad
            
            %ray_incident_unit_radius = tan(ray_incident_zenith);
            %[ray_incident_direction(ray,1),ray_incident_direction(ray,2)] = pol2cart(ray_incident_azimuth,ray_incident_unit_radius);
            %ray_incident_direction(ray,3) = -sqrt(1 - ray_incident_unit_radius^2); % z
            
            % associate a relative intensity with each ray
            
            %ray_intensity(ray,1) = 1;
            %calculate the unscaled power passing through the aperture
            %unscaled_power_in = unscaled_power_in + ray_intensity(ray) * cos(ray_incident_zenith/1000);
            unscaled_power_direct = unscaled_power_direct + ray_intensity(ray);
        end
        power_scale_direct = incident_power_direct / unscaled_power_direct;
        
        % calculate the direction vectors for for diffuse irradiance
        % these will all be unit vectors, away from the sun
        unscaled_power_diffuse = 0;
        for ray = number_of_rays_direct+1:number_of_rays_total
            % work in polar coordinates to calculate the random direction
            ray_incident_azimuth = rand * 2 * pi; % in rad
            ray_incident_zenith = asin(rand * sin(zenith_max_radians)); % in rad
            ray_incident_radius = sin(ray_incident_zenith); % from top down perspective
            
            % polar to cartesian
            [ray_incident_direction(ray,1),ray_incident_direction(ray,2)] = pol2cart(ray_incident_azimuth,ray_incident_radius);
            ray_incident_direction(ray,3) = -sqrt(1 - ray_incident_radius^2); % z
            
            % assume intensity is the same from all directions
            ray_intensity(ray,1) = 1;
            %calculate the unscaled power passing through the aperture
            %unscaled_power_in = unscaled_power_in + ray_intensity(ray) * cos(ray_incident_zenith/1000);
            unscaled_power_diffuse = unscaled_power_diffuse + ray_intensity(ray,1);
        end
        power_scale_diffuse = incident_power_diffuse / unscaled_power_diffuse;
        
        % calculate the ray incident points for both direct and diffuse
        ray_incident_point(:,1) = ray_aperture_point(:,1)-incident_ray_length*ray_incident_direction(:,1);
        ray_incident_point(:,2) = ray_aperture_point(:,2)-incident_ray_length*ray_incident_direction(:,2);
        ray_incident_point(:,3) = ray_aperture_point(:,3)-incident_ray_length*ray_incident_direction(:,3);
        
        
        %%% calculate mirror & ray intersect points %%%
        
        % does the ray intersect the mirror? true(1) or false(0)
        % initiate to false(0)
        mirror_intersect = zeros(number_of_rays_total,1);
        absorber_intersect = zeros(number_of_rays_total,1);
        
        % point of intersection between ray and mirror plane
        ray_intersect_point = zeros(number_of_rays_total,3);
        
        % point on reflected ray (not the intersect)
        ray_reflected_point = zeros(number_of_rays_total,3);
        ray_reflected_direction = zeros(number_of_rays_total,3);
        
        % point of intersection between reflected ray and absorber plane
        ray_absorber_point = zeros(number_of_rays_total,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform the ray tracing %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % also calculate power incident on absorber
        unscaled_absorber_power_direct = 0;
        unscaled_absorber_power_diffuse = 0;
        
        % Set ray length scalar
        ray_length_scalar = 2;
        
        % collector part
        if collector_design_linear_fresnel == true
            for ray = 1:number_of_rays_total
                
                % dispay ray number every 10,000 so that the user can see progress
                if(rem(ray,10000) == 0)
                    disp(['Calculating ray number ',num2str(ray)]);
                end
                
                for mirror = 1:number_of_mirrors
                    % if an intersect has already been detected, then quit loop
                    if(mirror_intersect(ray) == true)
                        break;
                    end
                    % find intersection between ray and mirror plane
                    m = dot((mirror_point(mirror,:) - ray_incident_point(ray,:)),mirror_normal(mirror,:)) / dot(ray_incident_direction(ray,:),mirror_normal(mirror,:));
                    ray_intersect_point(ray,:) = ray_incident_point(ray,:) + m*ray_incident_direction(ray,:);
                    
                    % calculate if the intersection is on the mirror
                    % there are several ways of doing this
                    % rotate the intersect so it's parallel to x and y axis
                    % this could be faster if precalculated
                    mirror_rotation_matrix = [cos(mirror_normal_direction(mirror)) -sin(mirror_normal_direction(mirror)) 0; sin(mirror_normal_direction(mirror)) cos(mirror_normal_direction(mirror)) 0; 0 0 1];
                    
                    rot_ray_intersect_point = (ray_intersect_point(ray,:)-mirror_centre_point(mirror,:))*mirror_rotation_matrix+mirror_centre_point(mirror,:);
                    
                    % and set mirror_intersect to true or false
                    %if (rot_ray_intersect_point(1) > rot_mirror_min_x(mirror)) && (rot_ray_intersect_point(1) < rot_mirror_max_x(mirror)) && (rot_ray_intersect_point(2) > rot_mirror_min_y(mirror)) && (rot_ray_intersect_point(2) < rot_mirror_max_y(mirror))
                    if (ray_intersect_point(ray,1) > mirror_min_x(mirror)) && (ray_intersect_point(ray,1) < mirror_max_x(mirror)) && (ray_intersect_point(ray,2) > mirror_min_y(mirror)) && (ray_intersect_point(ray,2) < mirror_max_y(mirror))
                        mirror_intersect(ray) = true;
                    else
                        mirror_intersect(ray) = false;
                    end
                    
                    if(mirror_intersect(ray) == true)
                        
                        % find the reflected point (on the other side of the mirror)
                        q = dot((mirror_point(mirror,:) - ray_incident_point(ray,:)),mirror_normal(mirror,:)) / dot(mirror_normal(mirror,:),mirror_normal(mirror,:));
                        ray_reflected_point_b = ray_incident_point(ray,:) + 2*q*mirror_normal(mirror,:);
                        
                        % find the reflected point (on the original side of the mirror)
                        ray_reflected_point(ray,:) = ray_length_scalar*(ray_intersect_point(ray,:) - ray_reflected_point_b)+ray_intersect_point(ray,:);
                        
                        ray_reflected_direction(ray,:) = ray_intersect_point(ray,:) - ray_reflected_point_b;
                    end
                end
                % deal with cases where there is no intersect with the mirror
                if (mirror_intersect(ray) == false)
                    % calculate the intersect with the floor
                    % find the reflected point (on the other side of the mirror)
                    m = dot((floor_point - ray_incident_point(ray,:)),floor_normal) / dot(ray_incident_direction(ray,:),floor_normal);
                    ray_intersect_point(ray,:) = ray_incident_point(ray,:) + m*ray_incident_direction(ray,:);
                end
            end
        end
        
        % absorber part
        if absorber_design_pipe == true
            absorber_pipe_array = zeros(absorber_pipe_length/absorber_pipe_y_res+1,360/absorber_pipe_theta_res);
            for ray = 1:number_of_rays_total
                
                % only continue if the ray intersects a mirror
                if(mirror_intersect(ray) == true)
                    
                    % solving a quadratic equation
                    % need to accommodate different rays in some of the matrices
                    % below
                    a = ray_reflected_direction(ray,1)^2 + ray_reflected_direction(ray,3)^2;
                    b = 2*(ray_reflected_direction(ray,1)*(ray_intersect_point(ray,1)-absorber_pipe_x)+ray_reflected_direction(ray,3)*(ray_intersect_point(ray,3)-absorber_pipe_z));
                    c = (ray_intersect_point(ray,1)-absorber_pipe_x)^2 + (ray_intersect_point(ray,3)-absorber_pipe_z)^2 - absorber_pipe_radius^2;
                    det = b^2 - 4 * a * c;
                    
                    if(det > 0)
                        
                        m = (-b - sqrt(det)) / (2 * a);
                        ray_absorber_point(ray,:) = ray_intersect_point(ray,:) + m*ray_reflected_direction(ray,:);
                        
                        if(imag(ray_absorber_point(ray,3)) ~= 0)
                            disp(['error - imaginary component in ray_absorber_point 3 detected']);
                        end
                        
                        if(imag(ray_absorber_point(ray,1)) ~= 0)
                            disp(['error - imaginary component in ray_absorber_point 1 detected']);

                        end
                        
                        absorber_theta = atan2d(ray_absorber_point(ray,3)-absorber_pipe_z,ray_absorber_point(ray,1)-absorber_pipe_x);
                        
                        
                        
                        if(absorber_theta < 0)
                            absorber_theta = absorber_theta + 360;
                        end
                        absorber_theta_index = ceil(absorber_theta/absorber_pipe_theta_res);
                        
                        if(absorber_theta_index == 0)
                            absorber_theta_index = 1;
                        end
                        absorber_y_index = ceil((ray_absorber_point(ray,2)-(absorber_pipe_y-0.5*absorber_pipe_length))/absorber_pipe_y_res);
                        
                        % see if the absorber point is within the length of the
                        % pipe
                        if (absorber_y_index > 0) && (absorber_y_index < (2 + absorber_pipe_length/absorber_pipe_y_res))
                            absorber_intersect(ray) = true;
                            % determine direct or diffuse for power calculations
                            if ray > number_of_rays_direct
                                % this is a diffuse ray
                                unscaled_absorber_power_diffuse = unscaled_absorber_power_diffuse + ray_intensity(ray);
                                absorber_pipe_array(absorber_y_index,absorber_theta_index) = absorber_pipe_array(absorber_y_index,absorber_theta_index) + ray_intensity(ray)*power_scale_diffuse;
                            else
                                % this is a direct ray
                                unscaled_absorber_power_direct = unscaled_absorber_power_direct + ray_intensity(ray);
                                absorber_pipe_array(absorber_y_index,absorber_theta_index) = absorber_pipe_array(absorber_y_index,absorber_theta_index) + ray_intensity(ray)*power_scale_direct;
                            end
                        else
                            absorber_intersect(ray) = false;
                        end
                    else
                        absorber_intersect(ray) = false;
                    end
                end
            end
        end
        
        if absorber_design_flat_panel == true
            absorber_panel_array = zeros(floor(absorber_panel_size_x/absorber_panel_x_res)+1,floor(absorber_panel_size_y/absorber_panel_y_res)+1);
            for ray = 1:number_of_rays_total
                
                % only continue if the ray intersects a mirror
                if(mirror_intersect(ray) == true)
                    
                    m = dot((absorber_point - ray_intersect_point(ray,:)),absorber_normal) / dot(ray_reflected_direction(ray,:),absorber_normal);
                    ray_absorber_point(ray,:) = ray_intersect_point(ray,:) + m*ray_reflected_direction(ray,:);
                    
                    % calculate if the intersection is on the absorber
                    % there are several ways of doing this
                    % rotate the intersect so it's parallel to x and y axis
                    rot_ray_absorber_point = (ray_absorber_point(ray,:)-absorber_panel_centre_point)*absorber_rotation_matrix+absorber_panel_centre_point;
                    
                    % and set absorber_intersect to true or false
                    if (rot_ray_absorber_point(1) > rot_absorber_min_x) && (rot_ray_absorber_point(1) < rot_absorber_max_x) && (rot_ray_absorber_point(2) > rot_absorber_min_y) && (rot_ray_absorber_point(2) < rot_absorber_max_y)
                        absorber_intersect(ray) = true;
                        % calculate the angle between the reflected ray and the
                        % absorber normal
                        %cos_angle = dot(ray_reflected_direction,absorber_normal)/(norm(ray_reflected_direction)*norm(absorber_normal));
                        %unscaled_absorber_power = unscaled_absorber_power + ray_intensity(ray);
                        absorber_x_index = ceil((absorber_panel_size_x/absorber_panel_x_res)*(rot_ray_absorber_point(1)-rot_absorber_min_x)/(rot_absorber_max_x-rot_absorber_min_x));
                        absorber_y_index = ceil((absorber_panel_size_y/absorber_panel_y_res)*(rot_ray_absorber_point(2)-rot_absorber_min_y)/(rot_absorber_max_y-rot_absorber_min_y));
                        if(absorber_x_index == 0)
                            absorber_x_index = 1;
                        end
                        if(absorber_y_index == 0)
                            absorber_y_index = 1;
                        end
                        if ray > number_of_rays_direct
                            % this is a diffuse ray
                            unscaled_absorber_power_diffuse = unscaled_absorber_power_diffuse + ray_intensity(ray);
                            absorber_panel_array(absorber_x_index,absorber_y_index) = absorber_panel_array(absorber_x_index,absorber_y_index) + ray_intensity(ray)*power_scale_diffuse;
                        else
                            % this is a direct ray
                            unscaled_absorber_power_direct = unscaled_absorber_power_direct + ray_intensity(ray);
                            absorber_panel_array(absorber_x_index,absorber_y_index) = absorber_panel_array(absorber_x_index,absorber_y_index) + ray_intensity(ray)*power_scale_direct;
                        end
                    else
                        absorber_intersect(ray) = false;
                    end
                end
            end
        end
        
        absorber_power_direct = unscaled_absorber_power_direct * power_scale_direct;
        absorber_power_diffuse = unscaled_absorber_power_diffuse * power_scale_diffuse;
        absorber_power_total = absorber_power_direct + absorber_power_diffuse;
    end % end of main ray tracing loop
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % HEAT TRANSFER CODE %
    %%%%%%%%%%%%%%%%%%%%%%
    
    % use Jacobi iteration and energy balance method, for steady state problem
    % use the same absorber array as for the optical part
    
    % make the assumption that the water dominates the heat capacity
    % and thermal equilibrium is established in the pipe on each time step
    
    % see if thermal calculations are required
    if(rem(time_seconds,time_step_thermal_calculation) == 0)
        
        disp(['Thermal calculation at time ',num2str(time_seconds),' s']);
        
        if absorber_design_pipe == true
            for i = 1:thermal_iterations % this is required for convergence
                for theta_index = 1:360/absorber_pipe_theta_res
                    for y_index = 1:absorber_pipe_length/absorber_pipe_y_res+1
                        % calculate power in to the element
                        % radiation exchange with ambient
                        power_in_0 = emissivity_thermal_pipe(y_index,theta_index)*element_area_pipe*SB*(air_ambient_temperature^4 - absorber_temperature_pipe(y_index,theta_index)^4);
                        % energy in from irradiance
                        power_in_1 = emissivity_optical_pipe(y_index,theta_index)*absorber_pipe_array(y_index,theta_index);
                        % convection to the ambient air
                        power_in_2 = h_c_air_copper*element_area_pipe*(air_ambient_temperature - absorber_temperature_pipe(y_index,theta_index));
                        % energy from conduction of neighbouring elements
                        if(theta_index ~= 360/absorber_pipe_theta_res)
                            power_in_3 = element_conductivity_pipe_tan*(absorber_temperature_pipe(y_index,theta_index+1)-absorber_temperature_pipe(y_index,theta_index));
                        else
                            power_in_3 = element_conductivity_pipe_tan*(absorber_temperature_pipe(y_index,1)-absorber_temperature_pipe(y_index,theta_index));
                        end
                        if(theta_index ~= 1)
                            power_in_4 = element_conductivity_pipe_tan*(absorber_temperature_pipe(y_index,theta_index-1)-absorber_temperature_pipe(y_index,theta_index));
                        else
                            power_in_4 = element_conductivity_pipe_tan*(absorber_temperature_pipe(y_index,360/absorber_pipe_theta_res)-absorber_temperature_pipe(y_index,theta_index));
                        end
                        if(y_index ~= 1)
                            power_in_5 = element_conductivity_pipe_y*(absorber_temperature_pipe(y_index-1,theta_index)-absorber_temperature_pipe(y_index,theta_index));
                        else
                            power_in_5 = 0;
                        end
                        if(y_index ~= absorber_pipe_length/absorber_pipe_y_res+1)
                            power_in_6 = element_conductivity_pipe_y*(absorber_temperature_pipe(y_index+1,theta_index)-absorber_temperature_pipe(y_index,theta_index));
                        else
                            power_in_6 = 0;
                        end
                        %power_in_3 = 0;
                        %power_in_4 = 0;
                        %power_in_5 = 0;
                        %power_in_6 = 0;
                       
                       
                       % from conduction of water, assumed a boundary condition
                       if process_continuous == true
                        power_in_7 = h_c_water_copper*element_area_pipe*(water_temperature_pipe(y_index,1)-absorber_temperature_pipe(y_index,theta_index));
                       end
                       if process_batch == true
                        power_in_7 = h_c_water_copper*element_area_pipe*(water_temperature_reactor-absorber_temperature_pipe(y_index,theta_index));
                       end   
                        % total power flux
                        power_in = power_in_0 + power_in_1 + power_in_2+ power_in_3 + power_in_4 + power_in_5 + power_in_6 + power_in_7;
                        absorber_temperature_pipe_2(y_index,theta_index) = absorber_temperature_pipe(y_index,theta_index) + thermal_factor*power_in/c_element_pipe;
                    end
                end
                absorber_temperature_pipe = absorber_temperature_pipe_2;
            end
            % energy balance to calculate new water temperatures
            power_in_total = 0;
            for y_index = 1:absorber_pipe_length/absorber_pipe_y_res+1
                power_in = 0;
                for theta_index = 1:360/absorber_pipe_theta_res
                        % calculate energy in to the element from outside
                        % only
                        % radiation exchange with ambient
                        power_in_0 = emissivity_thermal_pipe(y_index,theta_index)*element_area_pipe*SB*(air_ambient_temperature^4 - absorber_temperature_pipe(y_index,theta_index)^4);
                        % energy in from irradiance
                        power_in_1 = emissivity_optical_pipe(y_index,theta_index)*absorber_pipe_array(y_index,theta_index);
                        % convection to the ambient air
                        power_in_2 = h_c_air_copper*element_area_pipe*(air_ambient_temperature - absorber_temperature_pipe(y_index,theta_index)); 
                        power_in = power_in + power_in_0 + power_in_1 + power_in_2;
                end
                if process_continuous == true
                    water_temperature_pipe(y_index,1) = water_temperature_pipe(y_index,1) + time_step_thermal_calculation*power_in/c_water_element;
                end
                power_in_total = power_in_total + power_in;
            end
            if process_batch == true
                water_temperature_reactor = water_temperature_reactor + time_step_thermal_calculation*power_in_total/c_water_reactor;
            end
          
            % calculate power delivered to water
            % consider final section of the absorber tube
            %energy_in = (water_temperature(absorber_pipe_l/absorber_y_res+1,1) - water_ambient_temperature)*water_c;
            
            %water_power = energy_in/time_step_thermal_calculation; % because one section of water per second
            
            %water_flow = absorber_y_res*3.141*absorber_pipe_r.^2*density_water / time_step_thermal_calculation;
            
            %disp(['Power delivered to water ',num2str(water_power),' W']);
            
           % disp(['Water flow rate ',num2str(water_flow),' kg / s']);
        end
        
       
        if absorber_design_flat_panel == true
                %disp(['Time: ', int2str(second),'s']);
                %disp(['Reactor temperature: ', int2str(water_temperature_reactor),'K']);
                %disp(['Reactor a: ', int2str(reactor_a)]);
                
                
                % iteration
                for i = 1:thermal_iterations % this is required for convergence
                    for x_index = 1:absorber_panel_size_x/absorber_panel_x_res+1
                        for y_index = 1:absorber_panel_size_y/absorber_panel_y_res+1
                            % calculate power in to the element
                            % radiation exchange with ambient
                            power_in_0 = emissivity_thermal_panel(x_index,y_index)*element_area_panel*SB*(air_ambient_temperature^4 - absorber_temperature_panel(x_index,y_index)^4);
                            % energy in from irradiance
                            power_in_1 = emissivity_optical_panel(x_index,y_index)*absorber_panel_array(x_index,y_index);
                            % convection to the ambient air
                            power_in_2 = h_c_air_copper*element_area_panel*(air_ambient_temperature - absorber_temperature_panel(x_index,y_index));                                                    
                            if(x_index ~= 1)
                                power_in_3 = element_conductivity_panel_x*(absorber_temperature_panel(x_index-1,y_index)-absorber_temperature_panel(x_index,y_index));
                            else
                                power_in_3 = 0;
                            end
                            if(x_index ~= absorber_panel_size_x/absorber_panel_x_res+1)
                                power_in_4 = element_conductivity_panel_x*(absorber_temperature_panel(x_index+1,y_index)-absorber_temperature_panel(x_index,y_index));
                            else
                                power_in_4 = 0;
                            end
                            if(y_index ~= 1)
                                power_in_5 = element_conductivity_panel_y*(absorber_temperature_panel(x_index,y_index-1)-absorber_temperature_panel(x_index,y_index));
                            else
                                power_in_5 = 0;
                            end
                            if(y_index ~= absorber_panel_size_y/absorber_panel_y_res+1)
                                power_in_6 = element_conductivity_panel_y*(absorber_temperature_panel(x_index,y_index+1)-absorber_temperature_panel(x_index,y_index));
                            else
                                power_in_6 = 0;
                            end
                            
                            % from conduction of water, assumed a boundary condition
                            power_in_7 = h_c_water_copper*element_area_panel*(water_temperature_reactor-absorber_temperature_panel(x_index,y_index));
                            % total power flux
                            power_in = power_in_0 + power_in_1 + power_in_2+ power_in_3 + power_in_4 + power_in_5 + power_in_6 + power_in_7;
                            absorber_temperature_panel_2(x_index,y_index) = absorber_temperature_panel(x_index,y_index) + thermal_factor*power_in/c_element_panel;
                        end
                    end
                    absorber_temperature_panel = absorber_temperature_panel_2;
                end
                % energy balance to calculate new water temperatures
                power_in_total = 0; % in one second
                for y_index = 1:absorber_panel_size_y/absorber_panel_y_res+1
                    power_in = 0;
                    for x_index = 1:absorber_panel_size_x/absorber_panel_x_res+1
                        % energy in from ambient air
                        power_in_0 = emissivity_thermal_panel(x_index,y_index)*element_area_panel*SB*(air_ambient_temperature^4 - absorber_temperature_panel(x_index,y_index)^4);
                        % energy in from irradiance
                        power_in_1 = emissivity_optical_panel(x_index,y_index)*absorber_panel_array(x_index,y_index);
                        % convection to the ambient air
                        power_in_2 = h_c_air_copper*element_area_panel*(air_ambient_temperature - absorber_temperature_panel(x_index,y_index)); 
                        power_in = power_in + power_in_0 + power_in_1 + power_in_2;
                    end
                    if process_continuous == true
                        water_temperature_panel(y_index,1) = water_temperature_panel(y_index,1) + time_step_thermal_calculation*power_in/c_water_element;
                    end
                    power_in_total = power_in_total + power_in;
                end
                if process_batch == true
                    water_temperature_reactor = water_temperature_reactor + time_step_thermal_calculation*power_in_total/c_water_reactor;
                end
        end % end of flat panel thermal calculation
        
        
            % chemical calculations 2.0 (Olive trimings)
            
        if process_batch == true
            % chemical calculations
            
            % calculate reaction rates
            k_1 = k0_1*exp(-Ea_1/(R*water_temperature_reactor));
            k_2 = k0_2*exp(-Ea_2/(R*water_temperature_reactor));
            k_3 = k0_3*exp(-Ea_3/(R*water_temperature_reactor));
            k_4 = k0_4*exp(-Ea_4/(R*water_temperature_reactor));
            k_5 = k0_5*exp(-Ea_5/(R*water_temperature_reactor));
            
            % calculate carbon concentration rates
            C_B_dt = C_B*(-k_1-k_2-k_3);
            C_L_dt = (C_B*k_1) - (C_L*k_4) - (C_L^n*k_5); 
            C_G1_dt = C_B*k_2;
            C_HC1_dt = C_B*k_3;
            C_G2_dt = C_L*k_4;
            C_HC2_dt = (C_L^n)*k_5;
            
            % calculate new carbon concentrations
            C_B = C_B +(C_B_dt*time_step_thermal_calculation);
            C_L = C_L + (C_L_dt* time_step_thermal_calculation);
            C_G1 = C_G1 + (C_G1_dt*time_step_thermal_calculation);
            C_HC1 = C_HC1 + (C_HC1_dt*time_step_thermal_calculation);
            C_G2 = C_G2 + (C_G2_dt*time_step_thermal_calculation);
            C_HC2 = C_HC2 + (C_HC2_dt*time_step_thermal_calculation);

            if C_B > C_total
                C_B  = 0;
            end            
            if C_B < 0
                C_B  = 0;
            end
            X_B = C_B/C_total;
            
            if C_L > C_total
                C_L = C_total;
            end
            if C_L < 0
                C_L  = 0;
            end
            X_L= C_L/C_total;
                
            if C_G1 > C_total
                C_G1 = C_total;
            end
            if C_G1 < 0
                C_G1  = 0;
            end
            X_G1 =  C_G1/C_total;
            
            if C_HC1 > C_total

                C_HC1 = C_total;
            end
            if C_HC1 < 0
                C_HC1  = 0;
            end
            X_HC1 = C_HC1/C_total;
            
            if C_G2 > C_total
                C_G2 = C_total;
            end
            if C_G2 < 0
                C_G2  = 0;
            end
            X_G2 = C_G2/C_total;
            
            if C_HC2 > C_total
                C_HC2 = C_total;
            end
            if C_HC2 < 0
                C_HC2  = 0;
            end
            X_HC2 = C_HC2/C_total;                
                        
            
        end
        
        
        % still inside thermal calculation loop
        %if process_batch == true
        %    % chemical calculations
        %    reactor_da_dt = reactor_A * (1-reaction_a) * exp(-reactor_E/(R*water_temperature_reactor));
        %    reaction_a = reaction_a + reactor_da_dt*time_step_thermal_calculation;
        %    if reaction_a > 1
        %        reaction_a = 1;
        %    end
        %    if reaction_a < 0
        %        reaction_a = 0;
        %    end
        %end
        % data logging for all absorber types
        % inside thermal calculation loop
        log_index = log_index + 1;
        time_seconds_log(log_index,1) = time_seconds;
        time_minutes_log(log_index,1) = time_minutes;
        time_hours_log(log_index,1) = time_hours;
        reactor_temperature_log(log_index,1) = water_temperature_reactor;
        reactor_power_in_log(log_index,1) = power_in_total;
        
        global_irradiance_log(log_index,1) = direct_normal_irradiance;
        zenith_log(log_index,1) = zenith_radians;
        azimuth_log(log_index,1) = azimuth_radians;
      
        %reactor_a_log(log_index,1) = reaction_a;
        if process_batch == true
            C_B_log(log_index,1) = C_B;
            C_B_dt_log(log_index,1) = C_B_dt;
            C_L_log(log_index,1)= C_L;
            C_HC1_log(log_index,1)= C_HC1;
            C_HC2_log(log_index,1)= C_HC2;
            C_G1_log(log_index,1)= C_G1;
            C_G2_log(log_index,1) = C_G2;
            X_B_log(log_index,1) = X_B;
            X_L_log(log_index,1)= X_L;
            X_HC1_log(log_index,1)= X_HC1;
            X_HC2_log(log_index,1)= X_HC2;
            X_G1_log(log_index,1)= X_G1;
            X_G2_log(log_index,1) = X_G2;
        end
        
        absorber_incident_power_log(log_index,1) = absorber_power_total;
        aperture_incident_power_log(log_index,1) = incident_power_direct + incident_power_diffuse;
    end % end of main thermal calculation loop

    % see if water shuttle is required
    if process_continuous == true
        if(time_left_water_shuttle < 0)
            disp(['Water shuttle at time ',num2str(time_seconds),' s']);
            time_left_water_shuttle = time_left_water_shuttle + time_water_shuttle_pipe;  

            % shuttle water along absorber tube
            if absorber_design_pipe == true
                for y_index = absorber_pipe_length/absorber_pipe_y_res:-1:1
                    water_temperature_pipe(y_index+1,1) = water_temperature_pipe(y_index,1);
                end
                water_temperature_pipe(1,1) = water_ambient_temperature;
            end
            if absorber_design_flat_panel == true
                for y_index = absorber_panel_size_y/absorber_panel_y_res:-1:1
                    water_temperature_panel(y_index+1,1) = water_temperature_panel(y_index,1);
                end
                water_temperature_panel(1,1) = water_ambient_temperature;
            end
        end % end of water shuttle loop
        time_left_water_shuttle = time_left_water_shuttle - time_step_universal;
    end
end % end of main timestepping loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plots and final reporting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate power delivered to water
% consider final section of the absorber tube
%energy_in = (water_temperature_reactor - water_ambient_temperature)*water_c_reactor;
%water_power = energy_in/(time_seconds_stop-time_seconds_start);

%disp(['Final water temperature in reactor: ',num2str(water_temperature_reactor),' K']);
%disp(['Average power delivered to water in reactor: ',num2str(water_power),' W']);
            



% the Matlab line function requires separate vectors for x, y, and z
% set these up
x = zeros(1,2);
y = zeros(1,2);
z = zeros(1,2);

figure(1);
hold on;
axis equal;

% plot the incident rays
if(display_rays == true)
    for ray = 1:number_of_rays_total
        if ray_intensity(ray,1) > 0
        x(1) = ray_incident_point(ray,1);
        y(1) = ray_incident_point(ray,2);
        z(1) = ray_incident_point(ray,3);
        x(2) = ray_intersect_point(ray,1);
        y(2) = ray_intersect_point(ray,2);
        z(2) = ray_intersect_point(ray,3);
        
        if(mirror_intersect(ray) == true && absorber_intersect(ray) == true)
            line(x,y,z,'color','r');
        end
        if(mirror_intersect(ray) == true && absorber_intersect(ray) == false)
            line(x,y,z,'color','b');
        end
        if(mirror_intersect(ray) == false)
            line(x,y,z,'color','g');
        end
        end
    end
    
    % plot the reflected rays
    for ray = 1:number_of_rays_total
        if ray_intensity(ray,1) > 0
        if(mirror_intersect(ray) == true)
            x(1) = ray_intersect_point(ray,1);
            y(1) = ray_intersect_point(ray,2);
            z(1) = ray_intersect_point(ray,3);
            x(2) = ray_absorber_point(ray,1);
            y(2) = ray_absorber_point(ray,2);
            z(2) = ray_absorber_point(ray,3);
            % plot a green line
            if(absorber_intersect(ray) == true)
                line(x,y,z,'color','r');
            end
            %            if(absorber_intersect(ray) == false)
            %                line(x,y,z,'color','b');
            %            end
        end
        end
    end
end

% plot the mirrors
for mirror = 1:number_of_mirrors
    X(1) = mirror_corner_1(mirror,1);
    X(2) = mirror_corner_2(mirror,1);
    X(3) = mirror_corner_3(mirror,1);
    X(4) = mirror_corner_4(mirror,1);
    Y(1) = mirror_corner_1(mirror,2);
    Y(2) = mirror_corner_2(mirror,2);
    Y(3) = mirror_corner_3(mirror,2);
    Y(4) = mirror_corner_4(mirror,2);
    Z(1) = mirror_corner_1(mirror,3);
    Z(2) = mirror_corner_2(mirror,3);
    Z(3) = mirror_corner_3(mirror,3);
    Z(4) = mirror_corner_4(mirror,3);
    % plot a filled plane
    fill3(X,Y,Z,'y');
    hold on;
end

% plot the absorber tile
if absorber_design_flat_panel == true
    X(1) = absorber_corner_1(1);
    X(2) = absorber_corner_2(1);
    X(3) = absorber_corner_3(1);
    X(4) = absorber_corner_4(1);
    Y(1) = absorber_corner_1(2);
    Y(2) = absorber_corner_2(2);
    Y(3) = absorber_corner_3(2);
    Y(4) = absorber_corner_4(2);
    Z(1) = absorber_corner_1(3);
    Z(2) = absorber_corner_2(3);
    Z(3) = absorber_corner_3(3);
    Z(4) = absorber_corner_4(3);
    % plot a filled plane
    fill3(X,Y,Z,'k');
    hold on;
end

% plot the absorber pipe
if absorber_design_pipe == true
    for theta = 0:absorber_pipe_theta_res:360
        X(1) = absorber_pipe_x + absorber_pipe_radius * cosd(theta);
        X(2) = absorber_pipe_x + absorber_pipe_radius * cosd(theta + absorber_pipe_theta_res);
        X(3) = absorber_pipe_x + absorber_pipe_radius * cosd(theta + absorber_pipe_theta_res);
        X(4) = absorber_pipe_x + absorber_pipe_radius * cosd(theta);
        Y(1) = absorber_pipe_y - absorber_pipe_length * 0.5;
        Y(2) = absorber_pipe_y - absorber_pipe_length * 0.5;
        Y(3) = absorber_pipe_y + absorber_pipe_length * 0.5;
        Y(4) = absorber_pipe_y + absorber_pipe_length * 0.5;
        Z(1) = absorber_pipe_z + absorber_pipe_radius * sind(theta);
        Z(2) = absorber_pipe_z + absorber_pipe_radius * sind(theta + absorber_pipe_theta_res);
        Z(3) = absorber_pipe_z + absorber_pipe_radius * sind(theta + absorber_pipe_theta_res);
        Z(4) = absorber_pipe_z + absorber_pipe_radius * sind(theta);
        % plot a filled plane (mirror)
        fill3(X,Y,Z,'k');
        hold on;
    end
    for theta = 180:absorber_pipe_theta_res:360
        X(1) = absorber_pipe_x + absorber_pipe_radius * cosd(theta);
        X(2) = absorber_pipe_x + absorber_pipe_radius * cosd(theta + absorber_pipe_theta_res);
        X(3) = absorber_pipe_x + absorber_pipe_radius * cosd(theta + absorber_pipe_theta_res);
        X(4) = absorber_pipe_x + absorber_pipe_radius * cosd(theta);
        Y(1) = absorber_pipe_y - absorber_pipe_length*0.5;
        Y(2) = absorber_pipe_y - absorber_pipe_length*0.5;
        Y(3) = absorber_pipe_y + absorber_pipe_length*0.5;
        Y(4) = absorber_pipe_y + absorber_pipe_length*0.5;
        Z(1) = absorber_pipe_z + absorber_pipe_radius * sind(theta);
        Z(2) = absorber_pipe_z + absorber_pipe_radius * sind(theta + absorber_pipe_theta_res);
        Z(3) = absorber_pipe_z + absorber_pipe_radius * sind(theta + absorber_pipe_theta_res);
        Z(4) = absorber_pipe_z + absorber_pipe_radius * sind(theta);
       
        % plot a filled plane (mirror)
            fill3(X,Y,Z,'y')
            hold on;
         
    end
end

%Plot the lamps for the solar simulator

if light_source_solar_simulator == true
    
teta=0:0.01:2*pi ;

%lamp 1 
lamp_x_1=lamp_1_centre(1)+lamps_radius*cos(teta);
lamp_y_1=lamp_1_centre(2)+lamps_radius*sin(teta) ;
lamp_z_1 = lamp_1_centre(3)+zeros(size(lamp_x_1)) ;
patch(lamp_x_1,lamp_y_1,lamp_z_1,'k')
hold on
plot3(lamp_1_centre(1),lamp_1_centre(2),lamp_1_centre(3),'*r')

%lamp 2
lamp_x_2=lamp_2_centre(1)+lamps_radius*cos(teta);
lamp_y_2=lamp_2_centre(2)+lamps_radius*sin(teta) ;
lamp_z_2 = lamp_2_centre(3)+zeros(size(lamp_x_2)) ;
patch(lamp_x_2,lamp_y_2,lamp_z_2,'k')
hold on
plot3(lamp_2_centre(1),lamp_2_centre(2),lamp_2_centre(3),'*r')

%lamp 3
lamp_x_3=lamp_3_centre(1)+lamps_radius*cos(teta);
lamp_y_3=lamp_3_centre(2)+lamps_radius*sin(teta) ;
lamp_z_3 = lamp_3_centre(3)+zeros(size(lamp_x_3)) ;
patch(lamp_x_3,lamp_y_3,lamp_z_3,'k')
hold on
plot3(lamp_3_centre(1),lamp_3_centre(2),lamp_3_centre(3),'*r')

%lamp 4
lamp_x_4=lamp_4_centre(1)+lamps_radius*cos(teta);
lamp_y_4=lamp_4_centre(2)+lamps_radius*sin(teta) ;
lamp_z_4 = lamp_4_centre(3)+zeros(size(lamp_x_4)) ;
patch(lamp_x_4,lamp_y_4,lamp_z_4,'k')
hold on
plot3(lamp_4_centre(1),lamp_4_centre(2),lamp_4_centre(3),'*r')

%lamp 5
lamp_x_5=lamp_5_centre(1)+lamps_radius*cos(teta);
lamp_y_5=lamp_5_centre(2)+lamps_radius*sin(teta) ;
lamp_z_5 = lamp_5_centre(3)+zeros(size(lamp_x_5)) ;
patch(lamp_x_5,lamp_y_5,lamp_z_5,'k')
hold on
plot3(lamp_5_centre(1),lamp_5_centre(2),lamp_5_centre(3),'*r')

%lamp 6
lamp_x_6=lamp_6_centre(1)+lamps_radius*cos(teta);
lamp_y_6=lamp_6_centre(2)+lamps_radius*sin(teta) ;
lamp_z_6 = lamp_6_centre(3)+zeros(size(lamp_x_6)) ;
patch(lamp_x_6,lamp_y_6,lamp_z_6,'k')
hold on
plot3(lamp_6_centre(1),lamp_6_centre(2),lamp_6_centre(3),'*r')

%lamp 7
lamp_x_7=lamp_7_centre(1)+lamps_radius*cos(teta);
lamp_y_7=lamp_7_centre(2)+lamps_radius*sin(teta) ;
lamp_z_7 = lamp_7_centre(3)+zeros(size(lamp_x_7)) ;
patch(lamp_x_7,lamp_y_7,lamp_z_7,'k')
hold on
plot3(lamp_7_centre(1),lamp_7_centre(2),lamp_7_centre(3),'*r')

%lamp 8
lamp_x_8=lamp_8_centre(1)+lamps_radius*cos(teta);
lamp_y_8=lamp_8_centre(2)+lamps_radius*sin(teta) ;
lamp_z_8 = lamp_8_centre(3)+zeros(size(lamp_x_8)) ;
patch(lamp_x_8,lamp_y_8,lamp_z_8,'k')
hold on
plot3(lamp_8_centre(1),lamp_8_centre(2),lamp_8_centre(3),'*r')

end

% figure 2 shows incident irradiance

if absorber_design_pipe == true
    figure(2);
    hold on;
    pcolor(absorber_pipe_array'/element_area_pipe);
    shading flat;
    title('Irradiation on absorber surface (W m^-2)');
    xlabel('pipe length (element)');
    ylabel('theta resolution (element)');
    colorbar;
end

if absorber_design_flat_panel == true
    figure(2);
    hold on;
    pcolor(absorber_panel_array'/element_area_panel);
    shading flat;
    axis equal;
%     xticks([0 10 20]);
%     xticklabels({'0','10','20'})
%     yticks([0 50 100]);
%     yticklabels({'0','50','100'})
    title('Irradiation on absorber surface (W m^-2)');
    xlabel('pipe width (element)');
    ylabel('pipe length (element)');
    colorbar;
end

disp(['Direct power incident on aperture ',num2str(incident_power_direct), 'W']);
disp(['Diffuse power incident on aperture ',num2str(incident_power_diffuse), 'W']);
disp(['Direct power incident on absorber ',num2str(absorber_power_direct), 'W']);
disp(['Diffuse power incident on absorber ',num2str(absorber_power_diffuse), 'W']);

% figure 3 plots absorber temperature

if absorber_design_pipe == true
    figure(3);
    hold on;
    pcolor(absorber_temperature_pipe');
    shading flat;
    title('Absorber surface temperature (K)');
    xlabel('pipe length (element)');
    ylabel('theta resolution (element)');
    colorbar
end


if absorber_design_flat_panel == true
    figure(3);
    hold on;
    pcolor(absorber_temperature_panel');
    shading flat;
    title('Absorber surface temperature (K)');
    axis equal;
    xlabel('x (element)');
    ylabel('y (element)');
    colorbar
end

% figure 4 plots optical data against time
figure(4)
subplot(3,1,1);
hold on;
plot(time_hours_log,zenith_log*360/(2*pi));
xlabel('time (h)');
ylabel('zenith angle (degrees)');

subplot(3,1,2);
hold on;
plot(time_hours_log,azimuth_log*360/(2*pi));
xlabel('time (h)');
ylabel('azimuth angle (degrees)');

subplot(3,1,3);
hold on;
plot(time_hours_log,global_irradiance_log);
xlabel('time (h)');
ylabel('DNI (W m^-2)');

% figure 4 power data against time
figure(5)

subplot(3,1,1);
hold on;
plot(time_hours_log,aperture_incident_power_log);
xlabel('time (h)');
ylabel('power incident on aperture (W)');

subplot(3,1,2);
hold on;
plot(time_hours_log,absorber_incident_power_log);
xlabel('time (h)');
ylabel('power incident on absorber (W)');

subplot(3,1,3);
hold on;
plot(time_hours_log,reactor_power_in_log);
xlabel('time (h)');
ylabel('power delivered to reactor or water (W)');


% figure 5 plots reactor and chemical properties

if process_continuous == true
    figure(6)
    subplot(1,1,1);

    hold on;
    plot(1:absorber_pipe_length/absorber_pipe_y_res+1,water_temperature_pipe(:,1));
    title('Water temperature along absorber (at final time)');
    xlabel('pipe length (element)');
    ylabel('T (K)');
end

if process_batch == true
    figure(6)
    
    subplot(2,1,1);
    hold on;
    plot(time_hours_log,reactor_temperature_log);
    xlabel('time (h)');
    ylabel('reactor temperature (K)');
    
    %subplot(2,1,2);
    %hold on;
    %plot(time_hours_log,reactor_a_log);
    %xlabel('time (h)');
    %ylabel('a');

    subplot(2,1,2);
    plot (time_hours_log,C_B_log);
    hold on; plot(time_hours_log,C_L_log,'r'); plot(time_hours_log,C_HC1_log,'g'); plot(time_hours_log,C_HC2_log,'m'); plot(time_hours_log,C_G1_log,'k'); plot(time_hours_log,C_G2_log,'y'); hold off;
    legend (' initial biomass', 'liquid intermediate', 'hydochar', 'solid byproduct',  'gas1',  'gas2')
    xlabel('time (h)');
    ylabel('Carbon concentration  (mol/L)');
end

