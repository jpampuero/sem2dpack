function energies = sem2d_read_and_compute_energies(dir, x, Trup, Gc_sim, Heat_sim)
% Reads energies from file and computes the derivative of energy with respect to rupture length.
%
% Inputs:
%   dir  - Directory paths for energy output file
%   x    - array of fault coordinates
%   Trup - Time array for rupture front
%   Gc_sim - Rupture velocities
%   Heat_sim - Original time array
%
% Outputs:
%   energies - Struct containing energy derivatives

    % Read energy file
    en = readmatrix(strcat(char(dir),'\energy_sem2d.tab'), 'FileType', 'text');

    % Interpolate rupture front location
    t_en = en(:,1);
    rf_en = smoothdata(interp1(Trup, x, t_en), 'gaussian', 30);

    % Calculate energy derivatives
    dEp = calculate_energy_derivative(en(:,2), rf_en);
    dEk = calculate_energy_derivative(en(:,3), rf_en);
    if size(en,2) == 4
        % backwards compatibility with files before E_W energy term
        dEe = calculate_energy_derivative(en(:,4), rf_en);
    else
        dEe = calculate_energy_derivative(en(:,4) + en(:,5), rf_en);
    end

    dGc_sim = calculate_energy_derivative(Gc_sim, rf_en);
    dHeat_sim = calculate_energy_derivative(Heat_sim, rf_en);

    % Compile results into a struct
    energies = struct('t_en', t_en, 'rf_en', rf_en, 'dEp', dEp, 'dEk', dEk, 'dEe', dEe, ...
                     'dGc_sim', dGc_sim, 'dHeat_sim', dHeat_sim);
end

function derivative = calculate_energy_derivative(energy_values, rf_en)
% Computes the derivative of energy with respect to rupture length.
%
% Inputs:
%   energy_values - Array of energy values
%   rf_en         - Rupture front positions
%
% Output:
%   derivative    - Computed derivative with respect to rupture length

    derivative = zeros(size(energy_values));

    derivative(1) = (energy_values(2) - energy_values(1)) / (rf_en(2) - rf_en(1));
    derivative(2:end-1) = (energy_values(3:end) - energy_values(1:end-2)) ./ (rf_en(3:end) - rf_en(1:end-2));
    derivative(end) = (energy_values(end) - energy_values(end-1)) / (rf_en(end) - rf_en(end-1));

    derivative = abs(derivative);
end
