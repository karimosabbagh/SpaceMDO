function [m_prop,cost,S3_constraints] = S3_prop_struc_mass(delta_v_escape, delta_v_arrival, departure_date, arrival_date, Isp, m_SC, eta_FOV, IFOV)
    % Add subsystem paths
    currentFilePath = fileparts(mfilename('fullpath'));
    subsytems = fullfile(currentFilePath, '..', 'Setup');
    addpath(subsytems);

    % global m_SC;

    % departure_date = datetime(departure_date, "ConvertFrom", "datenum", "Format", 'yyyy-MM-dd');
    % arrival_date = datetime(arrival_date, "ConvertFrom", "datenum", "Format", 'yyyy-MM-dd');

    % define constants
    IFOV = IFOV/1000;
    g = 9.81/1000;              % km/s^2

    c_prop = 50;                % $/kg
    gamma = 0.8;                % 
    c_struct = 70;              % $/kg
    n = 1.1;                    % 
    c_isp = 10;                  % $/s
    lam = 1.1;                  % 
    c_FOV = 50;                %
    phi = 1.5;
    
    % constraint constants
    m_minshield = 100; % minimum shield mass required for all missions
    beta = 0.5; 
    D_avg = 657;                  % Radiation dose at 1 AU (mSv/year), local parameter
    alpha = 100;  % radiation absorbed / kg of material (mSv/kg)
    d_Earth = 1;                  % Earth's orbital radius (AU)
    d_Mars = 1.52;                % Mars' orbital radius (AU)
    a = (d_Earth + d_Mars) / 2;   % Semi-major axis (AU)
    e = (d_Mars - d_Earth) / (d_Mars + d_Earth); % Eccentricity
    m_payload = 500; % kg
     
    % objective function
    delta_v = delta_v_escape + delta_v_arrival; 
    m_prop = m_SC*(1-exp(-1*delta_v/(g*Isp)));
    m_structure = m_SC - m_prop;
    disp(m_structure);
    cost = c_prop*(m_prop)^gamma + c_struct*(m_structure)^n + c_isp*(Isp)^lam + c_FOV*(eta_FOV / IFOV)^phi;
   
    % constraints
    
    % calculate time of flight
    tof = determine_tof(departure_date, arrival_date);
    
   
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
    shield = m_minshield + (D_total/alpha) - beta*m_structure;  % beta*m_structure = m_minshield + alpha*(D_total)
    


    if isnan(shield)
        disp('g2 is NaN!')
    end 

   
    S3_constraints = shield;

end





