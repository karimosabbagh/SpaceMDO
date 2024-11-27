function tof = determine_tof(departure_date, arrival_date)
% Compute the Julian day number at 0 UT for departure and arrival dates, and the time of flight

    % Extract year, month, and day
    departure_year = year(departure_date);
    departure_month = month(departure_date);
    departure_day = day(departure_date);

    arrival_year = year(arrival_date);
    arrival_month = month(arrival_date);
    arrival_day = day(arrival_date);

    % Compute J0 values
    jd1 = 367*departure_year - fix(7*(departure_year + fix((departure_month + 9)/12))/4) + ...
        fix(275*departure_month/9) + departure_day + 1721013.5;
    jd2 = 367*arrival_year - fix(7*(arrival_year + fix((arrival_month + 9)/12))/4) + ...
        fix(275*arrival_month/9) + arrival_day + 1721013.5;

    % Determine TOF
    tof = jd2 - jd1;
end