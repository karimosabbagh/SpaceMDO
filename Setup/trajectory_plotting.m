%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              EXAMPLES                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXAMPLE 8.7%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global mu
% mu = 1.327124e11;
% deg = pi/180;
% 
% %...Input data
% planet_id = 3;
% year = 2003;
% month = 8;
% day = 27;
% hour = 12;
% minute = 0;
% second = 0;
% %...
% 
% %...Algorithm 8.1:
% [coe, r, v, jd] = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);
% 
% %...Convert the planet_id and month numbers into names for output:
% [month_name, planet_name] = month_planet_names(month, planet_id);
% 
% %...Echo the input data and output the solution to
% % the command window:
% fprintf('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
% fprintf('\n Example 8.7')
% fprintf('\n\n Input data:\n');
% fprintf('\n Planet: %s', planet_name)
% fprintf('\n Year : %g', year)
% fprintf('\n Month : %s', month_name)
% fprintf('\n Day : %g', day)
% fprintf('\n Hour : %g', hour)
% fprintf('\n Minute: %g', minute)
% fprintf('\n Second: %g', second)
% fprintf('\n\n Julian day: %11.3f', jd)
% 
% fprintf('\n\n');
% fprintf(' Orbital elements:')
% fprintf('\n');
% 
% fprintf('\n Angular momentum (km^2/s) = %g', coe(1));
% fprintf('\n Eccentricity = %g', coe(2));
% fprintf('\n Right ascension of the ascending node (deg) = %g', coe(3));
% fprintf('\n Inclination to the ecliptic (deg) = %g', coe(4));
% fprintf('\n Argument of perihelion (deg) = %g', coe(5));
% fprintf('\n True anomaly (deg) = %g', coe(6));
% fprintf('\n Semimajor axis (km) = %g', coe(7));
% fprintf('\n');
% 
% fprintf('\n Longitude of perihelion (deg) = %g', coe(8));
% fprintf('\n Mean longitude (deg) = %g', coe(9));
% fprintf('\n Mean anomaly (deg) = %g', coe(10));
% fprintf('\n Eccentric anomaly (deg) = %g', coe(11));
% 
% fprintf('\n\n');
% fprintf(' State vector:')
% fprintf('\n');
% 
% fprintf('\n Position vector (km) = [%g %g %g]', r(1), r(2), r(3))
% fprintf('\n Magnitude = %g\n', norm(r))
% fprintf('\n Velocity (km/s) = [%g %g %g]', v(1), v(2), v(3))
% fprintf('\n Magnitude = %g', norm(v))
% 
% fprintf('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXAMPLE 8.8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc
global mu
mu = 1.327124e11;
deg = pi/180;

%...Data declaration for Example 8.8:

%...Departure
planet_id_dep = 3;
year = 2024;
month_dep = 11;
day = 22;
hour = 0;
minute = 0;
second = 0;
depart = [planet_id_dep year month_dep day hour minute second];

%...Arrival
planet_id_arr = 4;
year = 2025;
month_arr = 10;
day = 20;
hour = 0;
minute = 0;
second = 0;
arrive = [planet_id_arr year month_arr day hour minute second];

%...

%...Algorithm 8.2:
[planet1, planet2, trajectory] = interplanetary(depart, arrive);

R1 = planet1(1,1:3);
Vp1 = planet1(1,4:6);
jd1 = planet1(1,7);

R2 = planet2(1,1:3);
Vp2 = planet2(1,4:6);
jd2 = planet2(1,7);

V1 = trajectory(1,1:3);
V2 = trajectory(1,4:6);

tof = jd2 - jd1;

%...Use Algorithm 4.2 to find the orbital elements of the spacecraft trajectory based on [Rp1, V1]...
coe = coe_from_sv(R1, V1, mu);
% ... and [R2, V2]
coe2 = coe_from_sv(R2, V2, mu);

%...Equations 8.94 and 8.95:
vinf1 = V1 - Vp1;
vinf2 = V2 - Vp2;

%...Echo the input data and output the solution to the command window:
fprintf('–––––––––––––––––––––––––––––––––––––––––––––––––––––')

fprintf('\n Example 8.8')
fprintf('\n\n Departure:\n');
fprintf('\n Planet: %d', depart(1))
fprintf('\n Year : %g', depart(2))
fprintf('\n Month : %d', depart(3))
fprintf('\n Day : %g', depart(4))
fprintf('\n Hour : %g', depart(5))
fprintf('\n Minute: %g', depart(6))
fprintf('\n Second: %g', depart(7))
fprintf('\n\n Julian day: %11.3f\n', jd1)
fprintf('\n Planet position vector (km) = [%g %g %g]', R1(1),R1(2), R1(3))
fprintf('\n Magnitude = %g\n', norm(R1))
fprintf('\n Planet velocity (km/s) = [%g %g %g]', Vp1(1), Vp1(2), Vp1(3))
fprintf('\n Magnitude = %g\n', norm(Vp1))
fprintf('\n Spacecraft velocity (km/s) = [%g %g %g]', V1(1), V1(2), V1(3))
fprintf('\n Magnitude = %g\n', norm(V1))
fprintf('\n v-infinity at departure (km/s) = [%g %g %g]', vinf1(1), vinf1(2), vinf1(3))
fprintf('\n Magnitude = %g\n', norm(vinf1))
fprintf('\n\n Time of flight = %g days\n', tof)
fprintf('\n\n Arrival:\n');
fprintf('\n Planet: %d', arrive(1))
fprintf('\n Year : %g', arrive(2))
fprintf('\n Month : %d', arrive(3))
fprintf('\n Day : %g', arrive(4))
fprintf('\n Hour : %g', arrive(5))
fprintf('\n Minute: %g', arrive(6))
fprintf('\n Second: %g', arrive(7))
fprintf('\n\n Julian day: %11.3f\n', jd2)
fprintf('\n Planet position vector (km) = [%g %g %g]', R2(1), R2(2), R2(3))
fprintf('\n Magnitude = %g\n', norm(R1))
fprintf('\n Planet velocity (km/s) = [%g %g %g]', Vp2(1), Vp2(2), Vp2(3))
fprintf('\n Magnitude = %g\n', norm(Vp2))
fprintf('\n Spacecraft Velocity (km/s) = [%g %g %g]', V2(1), V2(2), V2(3))
fprintf('\n Magnitude = %g\n', norm(V2))
fprintf('\n v-infinity at arrival (km/s) = [%g %g %g]', vinf2(1), vinf2(2), vinf2(3))
fprintf('\n Magnitude = %g', norm(vinf2))
fprintf('\n\n\n Orbital elements of flight trajectory:\n')
fprintf('\n Angular momentum (km^2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Right ascension of the ascending node (deg) = %g', coe(3)/deg)
fprintf('\n Inclination to the ecliptic (deg) = %g', coe(4)/deg)
fprintf('\n Argument of perihelion (deg) = %g', coe(5)/deg)
fprintf('\n True anomaly at departure (deg) = %g', coe(6)/deg)
fprintf('\n True anomaly at arrival (deg) = %g\n', coe2(6)/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))

% If the orbit is an ellipse, output the period:
if coe(2) < 1
    fprintf('\n Period (days) = %g', 2*pi/sqrt(mu)*coe(7)^1.5/24/3600)
end

fprintf('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       PLOTTING TRAJECTORY                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 8.8

% Extract orbital elements
h = coe(1);          % Angular momentum (km^2/s)
e = coe(2);          % Eccentricity
RA = coe(3);         % Right ascension of the ascending node (rad)
incl = coe(4);       % Inclination (rad)
w = coe(5);          % Argument of perihelion (rad)
theta1 = coe(6);     % True anomaly at departure (rad)
theta2 = coe2(6); % True anomaly at arrival (rad)
a = coe(7);          % Semimajor axis (km)

% Generate the trajectory points
num_points = 500; % Number of points for the trajectory
true_anomalies = linspace(theta1, theta2, num_points); % True anomaly range
r = (h^2 / mu) ./ (1 + e * cos(true_anomalies)); % Orbital radius (km)

% Convert to Cartesian coordinates in the orbital plane
x_orb = r .* cos(true_anomalies);
y_orb = r .* sin(true_anomalies);

% Rotate to the ecliptic plane using orbital elements
R_matrix = rotation_matrix(RA, incl, w); % Rotation matrix based on RA, incl, w
trajectory_points = R_matrix * [x_orb; y_orb; zeros(1, num_points)];

% Plot the solar system
figure; hold on;
title('Interplanetary Trajectory');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on;
axis equal;

% Plot the Sun at the origin
plot3(0, 0, 0, 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');

% Plot departure planet
plot3(R1(1), R1(2), R1(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'blue');
text(R1(1), R1(2), R1(3), ' Earth (Departure)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Plot arrival planet
plot3(R2(1), R2(2), R2(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
text(R2(1), R2(2), R2(3), ' Mars (Arrival)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot the spacecraft trajectory
plot3(trajectory_points(1, :), trajectory_points(2, :), trajectory_points(3, :), 'g-', 'LineWidth', 2);

% Legends
legend('Sun', 'Earth (Departure)', 'Mars (Arrival)', 'Trajectory', 'Location', 'best');

hold off;

% Function to compute the rotation matrix
function R = rotation_matrix(RA, incl, w)
    Rz_RA = [cos(RA), -sin(RA), 0; sin(RA), cos(RA), 0; 0, 0, 1];
    Rx_incl = [1, 0, 0; 0, cos(incl), -sin(incl); 0, sin(incl), cos(incl)];
    Rz_w = [cos(w), -sin(w), 0; sin(w), cos(w), 0; 0, 0, 1];
    R = Rz_RA * Rx_incl * Rz_w; % Combined rotation
end