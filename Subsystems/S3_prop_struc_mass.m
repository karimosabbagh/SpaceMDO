function [m_prop,m_structure,m_SC,Isp,cost,S3_constraints] = propellant_structure_mass(delta_v_escape,delta_v_arrival, departure_date, arrival_date)
    % Add subsystem paths
    currentFilePath = fileparts(mfilename('fullpath'));
    subsytems = fullfile(currentFilePath, '..', 'Setup');
    addpath(subsytems);

    % define constants
    m_payload = 1000;           % kg, similar to Mars Recon Orbiter
    g = 9.81;                   % m/s^2
    c_prop = 50;                % $/kg
    gamma = 0.8;                % 
    c_struct = 70;              % $/kg
    n = 0.75;                    % 
    c_isp = 40;                  % $/s
    lam = 1.1;                  % 
    d_Earth = 1;                  % Earth's orbital radius (AU)
    d_Mars = 1.52;                % Mars' orbital radius (AU)
    a = (d_Earth + d_Mars) / 2;   % Semi-major axis (AU)
    e = (d_Mars - d_Earth) / (d_Mars + d_Earth); % Eccentricity
    mu = 1.327e11;                % Solar gravitational parameter (km^3/s^2) 
    
    % constraint constants
    m_minshield = 100; % minimum shield mass required for all missions
    beta = 0.3; 
    D_avg = 657;                  % Radiation dose at 1 AU (mSv/year), local parameter
    alpha = 100;  % radiation absorbed / kg of material (mSv/kg)
     
    % objective function
    cost = c_prop*(m_prop)^gamma + c_struct*(m_structure)^n + c_isp*(Isp)^lam;
    m_SC = m_structure + m_prop + m_payload;

    % constraints

    % Relation between masses and velocity total delta V
    delta_v = delta_v_escape + delta_v_arrival; 
    
    % calculate time of flight
    tof = determine_tof(departure_date, arrival_date);

    % g1: Minimum propellant mass required
    g1 = m_SC*(1-exp(-1*delta_v/(g*Isp)))- m_prop;     % m_prop >= m_SC*(1-^(-delta_v/g*Isp))
                    
    % g2: Radiation shielding from structural mass
    
    n_points = 1000;
    times = linspace(0, tof, n_points); % Number of times across the transfer period

    M = 2 * pi * times / (2*tof); % Mean anomaly for each time

    E = zeros(size(M));
    for i = 1:length(M)
        E(i) = fzero(@(E) E - e * sin(E) - M(i), M(i));  % Solving for eccentric anomaly (numerically)
    end

    theta = 2 * atan(sqrt((1 + e) / (1 - e)) .* tan(E / 2)); % Finding true anomaly

    r = a * (1 - e^2) ./ (1 + e * cos(theta)); % Distance from Sun across the times

    D = D_avg * (1 ./ r).^2;
    D_total = trapz(times, D); % Integrate to find a total dose applied to spacecraft across flight
    g2 = m_minshield + D_total/(alpha) - beta*m_structure;  % beta*m_structure = m_minshield + alpha*(D_total)
    
    S3_constraints = [g1, g2];

end





