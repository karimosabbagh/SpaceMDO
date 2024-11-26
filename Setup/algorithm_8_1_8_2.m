%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             FUNCTIONS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% month_planet_names.m

% This function returns the name of the month and the planet corresponding, respectively, to the numbers "month_id" and "planet_id".
function [month, planet] = month_planet_names(month_id, planet_id)
    months = ['January  ' 'February ' 'March    ' 'April    ' 'May      ' 'June     ' 'July     ' 'August   ' 'September' 'October  ' 'November ' 'December '];
    planets = ['Mercury' 'Venus  ' 'Earth  ' 'Mars   ' 'Jupiter' 'Saturn ' 'Uranus ' 'Neptune' 'Pluto  '];

    month = months (month_id, 1:9);
    planet = planets(planet_id, 1:7);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sv_from_coe.m

% This function computes the classical orbital elements (coe) from the state vector (R,V) using Algorithm 4.1.
function [r, v] = sv_from_coe(coe,mu)
    h = coe(1);
    e = coe(2);
    RA = coe(3);
    incl = coe(4);
    w = coe(5);
    TA = coe(6);
    
    %...Equations 4.45 and 4.46 (rp and vp are column vectors):
    rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
    vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);
    
    %...Equation 4.34:
    R3_W = [ cos(RA) sin(RA) 0
    -sin(RA) cos(RA) 0
    0 0 1];
    
    %...Equation 4.32:
    R1_i = [1 0 0
    0 cos(incl) sin(incl)
    0 -sin(incl) cos(incl)];
    
    %...Equation 4.34:
    R3_w = [ cos(w) sin(w) 0
    -sin(w) cos(w) 0
    0 0 1];
    
    %...Equation 4.49:
    Q_pX = (R3_w*R1_i*R3_W)';
    
    %...Equations 4.51 (r and v are column vectors):
    r = Q_pX*rp;
    v = Q_pX*vp;
    
    %...Convert r and v into row vectors:
    r = r';
    v = v';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coe_from_sv.m

% This function computes the state vector (r,v) from the classical orbital elements (coe).
function coe = coe_from_sv(R,V,mu)

eps = 1.e-10;

r = norm(R);
v = norm(V);

vr = dot(R,V)/r;

H = cross(R,V);
h = norm(H);

%...Equation 4.7:
incl = acos(H(3)/h);

%...Equation 4.8:
N = cross([0 0 1],H);
n = norm(N);

%...Equation 4.9:
if n ~= 0
    RA = acos(N(1)/n);
    if N(2) < 0
        RA = 2*pi - RA;
    end
else
    RA = 0;
end

%...Equation 4.10:
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

%...Equation 4.12 (incorporating the case e = 0):
if n ~= 0
    if e > eps
        w = acos(dot(N,E)/n/e);
        if E(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end

%...Equation 4.13a (incorporating the case e = 0):
if e > eps
    TA = acos(dot(E,R)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end

%...Equation 4.62 (a < 0 for a hyperbola):
a=h^2/mu/(1 - e^2);

coe = [h e RA incl w TA a];

end %coe_from_sv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% kepler_E.m

% This function uses Newton’s method to solve Kepler’s equation E - e*sin(E) = M for the eccentric anomaly, given the eccentricity and the mean anomaly.
function E = kepler_E(e, M)
    %...Set an error tolerance:
    error = 1.e-8;
    
    %...Select a starting value for E:
    if M < pi
        E = M + e/2;
    else
        E = M - e/2;
    end
    
    %...Iterate on Equation 3.17 until E is determined to within the error tolerance:
    ratio = 1;
    while abs(ratio) > error
        ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
        E = E - ratio;
    end
end %kepler_E

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JO.m

%This function computes the Julian day number at 0 UT for any year between 1900 and 2100 using Equation 5.48.
function j0 = J0(year, month, day)

    j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) + fix(275*month/9) + day + 1721013.5;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lambert.m

%This function solves Lambert's problem.
function [V1, V2] = lambert(R1, R2, t, string)
    global mu
    
    %...Magnitudes of R1 and R2:
    r1 = norm(R1);
    r2 = norm(R2);
    
    c12 = cross(R1, R2);
    theta = acos(dot(R1,R2)/r1/r2);
    
    %...Determine whether the orbit is prograde or retrograde:
    if nargin < 4 jj (strcmp(string,'retro') & (strcmp(string,'pro')))
        string = 'pro';
        fprintf('\n ** Prograde trajectory assumed.\n')
    end
    
    if strcmp(string,'pro')
        if c12(3) <= 0
            theta = 2*pi - theta;
        end
    elseif strcmp(string,'retro')
        if c12(3) >= 0
            theta = 2*pi - theta;
        end
    end
    
    %...Equation 5.35:
    A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));
    
    %...Determine approximately where F(z,t) changes sign, and use that value of z as the starting value for Equation 5.45:
    z = -100;
    while F(z,t) < 0
        z = z + 0.1;
    end
    
    %...Set an error tolerance and a limit on the number of iterations:
    tol = 1.e-8;
    nmax = 5000;
    
    %...Iterate on Equation 5.45 until z is determined to within the error tolerance:
    ratio = 1;
    n = 0;
    while (abs(ratio) > tol) & (n <= nmax)
        n = n + 1;
        ratio = F(z,t)/dFdz(z);
        z = z - ratio;
    end
    
    %...Report if the maximum number of iterations is exceeded:
    if n >= nmax
        fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax)
    end
    
    %...Equation 5.46a:
    f = 1 - y(z)/r1;
    
    %...Equation 5.46b:
    g = A*sqrt(y(z)/mu);
    
    %...Equation 5.46d:
    gdot = 1 - y(z)/r2;
    
    %...Equation 5.28:
    V1 = 1/g*(R2 - f*R1);
    
    %...Equation 5.29:
    V2 = 1/g*(gdot*R2 - R1);
    
    return
    
    % 
    % Subfunctions used in the main body:
    % 
    
    %...Equation 5.38:
    function dum = y(z)
    dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
    end
    
    %...Equation 5.40:
    function dum = F(z,t)
        dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
    end
    
    %...Equation 5.43:
    function dum = dFdz(z)
        if z == 0
            dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
        else
            dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) + A*sqrt(C(z)/y(z)));
        end
    end
    
    %...Stumpff functions:
    function dum = C(z)
        dum = stumpC(z);
    end
    
    function dum = S(z)
        dum = stumpS(z);
    end
end %lambert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stumpS.m

% This function evaluates the Stumpff function S(z) according to Equation 3.52.
function s = stumpS(z)
    if z > 0
        s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stumpC.m

% This function evaluates the Stumpff function C(z) according to Equation 3.53.
function c = stumpC(z)
    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
    else
        c = 1/2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% planet_elements_and_sv.m

% This function calculates the orbital elements and the state vector of a planet from the date (year, month, day) and universal time (hour, minute, second).
function [coe, r, v, jd] = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)
    global mu
    deg = pi/180;
    
    %...Equation 5.48:
    j0 = J0(year, month, day);
    
    ut = (hour + minute/60 + second/3600)/24;
    
    %...Equation 5.47
    jd = j0 + ut;
    
    %...Obtain the data for the selected planet from Table 8.1:
    [J2000_coe, rates] = planetary_elements(planet_id);
    
    %...Equation 8.93a:
    t0 = (jd - 2451545)/36525;
    
    %...Equation 8.93b:
    elements = J2000_coe + rates*t0;
    
    a = elements(1);
    e = elements(2);
    
    %...Equation 2.71:
    h = sqrt(mu*a*(1 - e^2));
    
    %...Reduce the angular elements to within the range 0 - 360 degrees:
    incl = elements(3);
    RA = zero_to_360(elements(4));
    w_hat = zero_to_360(elements(5));
    L = zero_to_360(elements(6));
    w = zero_to_360(w_hat - RA);
    M = zero_to_360((L - w_hat));
    
    %...Algorithm 3.1 (for which M must be in radians)
    E = kepler_E(e, M*deg);
    
    %...Equation 3.13 (converting the result to degrees):
    TA = zero_to_360(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))/deg);
    
    coe = [h e RA incl w TA a w_hat L M E/deg];
    
    %...Algorithm 4.5 (for which all angles must be in radians):
    [r, v] = sv_from_coe([h e RA*deg incl*deg w*deg TA*deg],mu);

return

% This function extracts a planet’s J2000 orbital elements and centennial rates from Table 8.1.
function [J2000_coe, rates] = planetary_elements(planet_id)
    J2000_elements = [ 0.38709893 0.20563069 7.00487 48.33167 77.45645 252.25084
    0.72333199 0.00677323 3.39471 76.68069 131.53298 181.97973
    1.00000011 0.01671022 0.00005 -11.26064 102.94719 100.46435
    1.52366231 0.09341233 1.85061 49.57854 336.04084 355.45332
    5.20336301 0.04839266 1.30530 100.55615 14.75385 34.40438
    9.53707032 0.05415060 2.48446 113.71504 92.43194 49.94432
    19.19126393 0.04716771 0.76986 74.22988 170.96424 313.23218
    30.06896348 0.00858587 1.76917 131.72169 44.97135 304.88003
    39.48168677 0.24880766 17.14175 110.30347 224.06676 238.92881];
    
    cent_rates = [ 0.00000066 0.00002527 -23.51 -446.30 573.57 538101628.29
    0.00000092 -0.00004938 -2.86 -996.89 -108.80 210664136.06
    -0.00000005 -0.00003804 -46.94 -18228.25 1198.28 129597740.63
    -0.00007221 0.00011902 -25.47 -1020.19 1560.78 68905103.78
    0.00060737 -0.00012880 -4.15 1217.17 839.93 10925078.35
    -0.00301530 -0.00036762 6.11 -1591.05 -1948.89 4401052.95
    0.00152025 -0.00019150 -2.09 -1681.4 1312.56 1542547.79
    -0.00125196 0.00002514 -3.64 -151.25 -844.43 786449.21
    -0.00076912 0.00006465 11.07 -37.33 -132.25 522747.90];
    
    J2000_coe = J2000_elements(planet_id,:);
    rates = cent_rates(planet_id,:);
    
    %...Convert from AU to km:
    au = 149597871;
    J2000_coe(1) = J2000_coe(1)*au;
    rates(1) = rates(1)*au;
    
    %...Convert from arcseconds to fractions of a degree:
    rates(3:6) = rates(3:6)/3600;
    end %planetary_elements
    
    function y = zero_to_360(x)
        if x >= 360
            x = x - fix(x/360)*360;
        elseif x < 0
            x = x - (fix(x/360) - 1)*360;
    end
    
    y = x;
    
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interplanetary.m

% This function determines the spacecraft trajectory from the sphere of influence of planet 1 to that of planet 2 using Algorithm 8.2
function [planet1, planet2, trajectory] = interplanetary(depart, arrive)
    global mu
    
    planet_id = depart(1);
    year = depart(2);
    month = depart(3);
    day = depart(4);
    hour = depart(5);
    minute = depart(6);
    second = depart(7);
    
    %...Use Algorithm 8.1 to obtain planet 1’s state vector (don’t need its orbital elements ["dum"]):
    [dum, Rp1, Vp1, jd1] = planet_elements_and_sv (planet_id, year, month, day, hour, minute, second);
    
    planet_id = arrive(1);
    year = arrive(2);
    month = arrive(3);
    day = arrive(4);
    hour = arrive(5);
    minute = arrive(6);
    second = arrive(7);
    
    %...Likewise use Algorithm 8.1 to obtain planet 2’s state vector:
    [dum, Rp2, Vp2, jd2] = planet_elements_and_sv (planet_id, year, month, day, hour, minute, second);
    
    tof = (jd2 - jd1)*24*3600;
    
    %...Patched conic assumption:
    R1 = Rp1;
    R2 = Rp2;
    
    %...Use Algorithm 5.2 to find the spacecraft’s velocity at departure and arrival, assuming a prograde trajectory:
    [V1, V2] = lambert(R1, R2, tof, 'pro');
    
    planet1 = [Rp1, Vp1, jd1];
    planet2 = [Rp2, Vp2, jd2];
    trajectory = [V1, V2];

end %interplanetary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% planet_state_for_date.m

% This function calculates the state vector of a planet from the date (year, month, day) and universal time (hour, minute, second).
function [planet_state_vector] = planet_state_for_date(planet_and_date)
    global mu
    
    planet_id = planet_and_date(1);
    year = planet_and_date(2);
    month = planet_and_date(3);
    day = planet_and_date(4);
    hour = planet_and_date(5);
    minute = planet_and_date(6);
    second = planet_and_date(7);
    
    %...Use Algorithm 8.1 to obtain the state vector of the planet at the given date
    [dum, Rp, Vp, jd] = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);
    
    % Return the state (position and velocity)
    planet_state_vector = [Rp, Vp, jd];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rotation_matrix.m

% This function computes the rotation matrix for plotting trajectories
function R = rotation_matrix(RA, incl, w)
    Rz_RA = [cos(RA), -sin(RA), 0; sin(RA), cos(RA), 0; 0, 0, 1];
    Rx_incl = [1, 0, 0; 0, cos(incl), -sin(incl); 0, sin(incl), cos(incl)];
    Rz_w = [cos(w), -sin(w), 0; sin(w), cos(w), 0; 0, 0, 1];
    R = Rz_RA * Rx_incl * Rz_w; % Combined rotation
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        PLANET STATE VECTORS                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc
global mu
mu = 1.327124e11;
deg = pi/180;

%same dates for plotting orbits
departure_start_date = datetime(2024, 11, 21); % Nov 21, 2024
departure_end_date = datetime(2027, 6, 11); % June 11, 2027
arrival_start_date = datetime(2024, 11, 21); % Nov 21, 2024
arrival_end_date = datetime(2027, 6, 11); % Jun 11, 2027

% Create date vectors (days)
departure_dates = departure_start_date:departure_end_date;
arrival_dates = arrival_start_date:arrival_end_date;


result1 = {};
result2 = {};

% Loop over each departure date
for dep_idx = 1:length(departure_dates)

    depart_date = departure_dates(dep_idx);
    planet_and_date = [3 depart_date.Year depart_date.Month depart_date.Day 0 0 0]; % Earth (planet 3)
    [planet_state_vector] = planet_state_for_date(planet_and_date);
    Rp = planet_state_vector(1,1:3);
    Vp = planet_state_vector(1,4:6);

    result1{end+1, 1} = depart_date; % Departure date
    result1{end, 2} = Rp(1); % Earth Position X
    result1{end, 3} = Rp(2); % Earth Position Y
    result1{end, 4} = Rp(3); % Earth Position Z
    result1{end, 5} = norm(Rp); % Earth Position Magnitude
    result1{end, 6} = Vp(1); % Earth Velocity X
    result1{end, 7} = Vp(2); % Earth Velocity Y
    result1{end, 8} = Vp(3); % Earth Velocity Z
    result1{end, 9} = norm(Vp); % Earth Velocity Magnitude
end

% Loop over each arrival date
for arr_idx = 1:length(arrival_dates)
    arrival_date = arrival_dates(arr_idx);
    planet_and_date = [4 arrival_date.Year arrival_date.Month arrival_date.Day 0 0 0]; % Mars (planet 4)
    [planet_state_vector] = planet_state_for_date(planet_and_date);
    Rp = planet_state_vector(1,1:3);
    Vp = planet_state_vector(1,4:6);

    result2{end+1, 1} = arrival_date; % Departure date
    result2{end, 2} = Rp(1); % Mars Position X
    result2{end, 3} = Rp(2); % Mars Position Y
    result2{end, 4} = Rp(3); % Mars Position Z
    result2{end, 5} = norm(Rp); % Mars Position Magnitude
    result2{end, 6} = Vp(1); % Mars Velocity X
    result2{end, 7} = Vp(2); % Mars Velocity Y
    result2{end, 8} = Vp(3); % Mars Velocity Z
    result2{end, 9} = norm(Vp); % Mars Velocity Magnitude
end

earth_result_table = cell2table(result1, 'VariableNames', {'DepartureDate', 'Earth_Position_X', 'Earth_Position_Y', 'Earth_Position_Z', 'Earth_Position_Magnitude', 'Earth_Velocity_X', 'Earth_Velocity_Y', 'Earth_Velocity_Z', 'Earth_Velocity_Magnitude'});
mars_result_table = cell2table(result2, 'VariableNames', {'ArrivalDate', 'Mars_Position_X', 'Mars_Position_Y', 'Mars_Position_Z', 'Mars_Position_Magnitude', 'Mars_Velocity_X', 'Mars_Velocity_Y', 'Mars_Velocity_Z', 'Mars_Velocity_Magnitude'});

disp(earth_result_table);
disp(mars_result_table);

writetable(earth_result_table, 'earth_result_table.csv');
writetable(mars_result_table, 'mars_result_table.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          ANIMATE ORBITS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Earth and Mars positions and dates
Earth_X = earth_result_table.Earth_Position_X;
Earth_Y = earth_result_table.Earth_Position_Y;
Earth_Z = earth_result_table.Earth_Position_Z;
Mars_X = mars_result_table.Mars_Position_X;
Mars_Y = mars_result_table.Mars_Position_Y;
Mars_Z = mars_result_table.Mars_Position_Z;
Earth_Dates = earth_result_table.DepartureDate; 
Mars_Dates = mars_result_table.ArrivalDate; 

% Plot
figure;
hold on;
grid on;
axis equal;
xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');
title('Earth and Mars Motion Around the Sun');

sunPlot = plot3(0, 0, 0, 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun'); % Sun
earthTrail = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Earth orbit');
marsTrail = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Mars orbit');
earthPlot = plot3(NaN, NaN, NaN, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');
marsPlot = plot3(NaN, NaN, NaN, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'DisplayName', 'Mars');

dateText = text(0.1, 0.9, 0.9, '', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

legend([earthPlot, marsPlot, earthTrail, marsTrail, sunPlot], 'Location', 'northeastoutside');

% Animate
videoFileName = 'earth_mars_orbit_animation.mp4';
v = VideoWriter(videoFileName, 'MPEG-4'); 
v.FrameRate = 30; 
open(v); 
numFrames = min(length(Earth_X), length(Mars_X));
for i = 1:numFrames
    set(earthPlot, 'XData', Earth_X(i), 'YData', Earth_Y(i), 'ZData', Earth_Z(i));
    set(marsPlot, 'XData', Mars_X(i), 'YData', Mars_Y(i), 'ZData', Mars_Z(i));
    set(earthTrail, 'XData', Earth_X(1:i), 'YData', Earth_Y(1:i), 'ZData', Earth_Z(1:i));
    set(marsTrail, 'XData', Mars_X(1:i), 'YData', Mars_Y(1:i), 'ZData', Mars_Z(1:i));
  
    currentDate = Earth_Dates(i); 
    dateStr = datestr(currentDate, 'yyyy-mm-dd');
    set(dateText, 'String', ['Date: ' dateStr]);
    
    frame = getframe(gcf); 
    writeVideo(v, frame); 
    pause(0.005);
end

close(v);