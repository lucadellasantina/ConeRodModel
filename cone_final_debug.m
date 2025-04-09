intensities = 50;

% Parameters of the cone cell (Kourennyi 2004)
C_cone = 10; % Capacitance in pF
ECa_cone = 40; % Calcium reversal potential in mV
EKv_cone = -80; % Kv reversal potential in mV
Eh_cone = -32.5; % Ih reversal potential in mV
EKCa_cone = -80; % KCa reversal potential in mV
ECl_cone = -45; % Cl reversal potential in
EL_cone = -63; % Leak reversal potential in mV

gCa_cone = 4.92; % Calcium conductance in nS
gKv_cone = 2; % Kv conductance in nS
gh_cone = 3.5; % Ih conductance in nS
gKCa_cone = 0.5; % KCa conductance in nS
gCl_cone = 6.5; % Cl conductance in nS 
gL_cone = 1; % Leak conductance in nS

aoCa = 3.1; % Calcium activation constant (1/s)
SCa = 5.7; % Calcium scaling factor in mV
VhalfCa = -16.6; % Calcium half activation voltage in mV
Cainf_cone = 0.001; % Steady-state calcium concentration in mM (Peak value 0.0005 – 0.001 mM )

% Intracellular calcium concentration (Liu & Kourennyi 2004)
radius_cone = 7.8; % microns
volume_cone = (4/3) * pi * radius_cone^3; % volume in um^3
z = 2; % Valence of calcium ions
tau_Cai_cone = 20; % Time constant for calcium decay in ms

% Conversion factors
fI = 10e-6; % Conversion factor for ICa (already in pA)
fF = 10e9; % Conversion factor for Faraday's constant
F = 96485; % Faraday constant (C/mol)

% Cone photocurrent parameters
Idark_cone = -20; % Dark current in pA

% Photocurrent model parameters
Iamp_max = 30; % pA maximal amplitude of photocurrent
Vamp_max = -41; % Maximum voltage amplitude
n = 1;         % Hill coefficient
R_half = 400;   % flash intensity at half-maximal amplitude
R_sat = 100000;   % saturating flash intensity

% Time parameters
% TODO: Investigate why model breaks if dt>0.01
dt = 0.1; % Time step in ms (same as 1e-5 s)
t_end = 1500; % Total simulation time in ms
t_del = 500; % stimulus trigered offset time (ms)
num_time_steps = t_end / dt;  
 
% Intensities to simulate
num_intensities = length(intensities);

time = linspace(0, t_end - dt, num_time_steps);  % Exclude the final step at t_end

% Preallocate variables for currents and states
V_all = zeros(num_intensities, num_time_steps); % Voltage for all intensities
Iphoto_all = zeros(num_intensities, num_time_steps); % Photocurrent for all intensities
ICa_all = zeros(num_intensities, num_time_steps); % Calcium current
IKv_all = zeros(num_intensities, num_time_steps); % Kv current
ICl_all = zeros(num_intensities, num_time_steps); % Cl current
IL_all = zeros(num_intensities, num_time_steps); % Leak current
I_noise = zeros(1, num_time_steps); % Preallocate noise current array
Ih_all = zeros(num_intensities, num_time_steps); % Ih current for all intensities
IKCa_all = zeros(num_intensities, num_time_steps); % IKCa current for all intensities
Cai_all = zeros(num_intensities, num_time_steps); % Initial intracellular calcium concentration in mM

nh_all = zeros(num_intensities, num_time_steps); % Ih activation variable
mKv_all = zeros(num_intensities, num_time_steps); % Kv activation variable
hKv_all = zeros(num_intensities, num_time_steps); % Kv inactivation variable
nCa_all = zeros(num_intensities, num_time_steps); % Calcium activation variable

function Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c, t_del, num_time_steps, dt)
    % Compute the original Iphoto
    Iphoto = Iamp .* ((1 - exp(-(time - t_del) / tau1))) .* (time > t_del) - ...
             Iamp .* (1 ./ (1 + exp(-(time - t_del - c) / tau2))) .* (time > t_del);

    % Preallocate noise array
    I_noise_photo = zeros(1, num_time_steps);

    sigma_noise = 1;
    tau_noise = 1;

    % Generate noise for Iphoto
    for i = 1:num_time_steps - 1
        noise = sigma_noise * sqrt(2 / tau_noise) * randn() / sqrt(dt);
        I_noise_photo(i + 1) = I_noise_photo(i) + dt * (-I_noise_photo(i) / tau_noise + noise);
    end

    % Add noise to Iphoto
    Iphoto = Iphoto + I_noise_photo;

    % Clamp Iphoto to ensure it does not go below zero
    %Iphoto(Iphoto < 0) = 0;
end

% Main simulation loop
for intensity_index = 1:num_intensities
    R_star = intensities(intensity_index) ;
    
    % Initial membrane parameters
    V = -41 * ones(1, num_time_steps); % Voltage for each intensity  (mV)
    mKv = 0.3709 * ones(1, num_time_steps); % Kv activation variable
    hKv = 0.9998 * ones(1, num_time_steps); % Kv inactivation variable
    nCa = 0.0115 * ones(1, num_time_steps); % Calcium activation variable
    nh = 0.0877 * ones(1, num_time_steps); % Ih activation variable
    Cai = 0.00002 * ones(1, num_time_steps); % Initial intracellular calcium concentration in mM (Rest Value 0.00002 – 0.00005 mM)
    
    
    Ih_cone = zeros(1, num_time_steps); % Reset Ih current for each intensity
    ICa_cone = zeros(1, num_time_steps); % Reset calcium current for each intensity
    IKv_cone = zeros(1, num_time_steps); % Reset Kv current for each intensity
    IKCa_cone = zeros(1, num_time_steps); % Reset IKCa current for each intensity
    ICl_cone = zeros(1, num_time_steps); % Reset Cl current for each intensity
    IL_cone = zeros(1, num_time_steps); % Reset leak current for each intensity

    % Hill equation for Iamp
    R_star = max(R_star, 1e-6); % Ensure R_star_i > 0 to avoid log(0)
    Iamp = Iamp_max * (R_star^n / (R_star^n + R_half^n));

     % Time constants and c based on R_star
     if R_star < 40000
        c = (5.05 * log(R_star / R_half) + 300)/2;
        tau1 = (-5.263 * log(R_star / R_half) + 26.235)/2;
        tau2 = (2.392 * log(R_star / R_half) + 236.09)/3;
    end

    % Compute Iphoto for all time steps
    Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c, t_del, num_time_steps, dt);

    % Add Iphoto to total photocurrent at each time step
    Iphoto_all(intensity_index, :) = Iphoto;

    % Main simulation loop
    for i = 1:num_time_steps-1
        % Update gating variables
        amKv =  (5 * (V(i) - 100) ./ (1 - exp((100 - V(i)) / 42)));
        bmKv =  (9 * exp((20 - V(i)) / 40));

        ahKv =  (0.15 * exp(-V(i) / 22));
        bhKv =  (0.4125 ./ (1 + exp((10 - V(i)) / 7)));

        anh =  (0.001 * 18 ./ (1 + exp((V(i) + 88) / 12)));
        bnh =  (0.001 * 18 ./ (1 + exp((-V(i) - 18) / 19)));

        anCa =  (aoCa * exp((V(i) - VhalfCa) / (2 * SCa)));
        bnCa =  (aoCa * exp((-V(i) + VhalfCa) / (2 * SCa)));

        hKCa = Cai(i) ./ (Cai(i) + 0.3);
        mCl =  1 ./ (1 + exp(0.00037 - Cai(i)) ./ 0.00009);

        % Update currents
        ICa_cone(i) = gCa_cone * nCa(i) * (V(i) - ECa_cone); % Calcium current in pA
        IKv_cone(i) = gKv_cone * mKv(i)^3 * hKv(i) * (V(i) - EKv_cone); % Kv current in pA
        IKCa_cone(i) = gKCa_cone * hKCa * (V(i) - EKCa_cone); % KCa current in pA
        ICl_cone(i) = gCl_cone * mCl * (V(i) - ECl_cone); % Cl current in pA
        IL_cone(i) = gL_cone * (V(i) - EL_cone); % Leak current in pA

        % Update Ih current
        Ih_cone(i) = gh_cone * nh(i) * (V(i) - Eh_cone); % Ih current in pA

        % Update membrane voltage with photocurrent from Iphoto
        % dV_dt = (-IKv_cone(i) ...
        %         - Ih_cone(i) ...
        %         - ICa_cone(i) ...
        %         - IKCa_cone(i) ...
        %         - ICl_cone(i) ...
        %         - IL_cone(i) ...
        %         - Idark_cone - Iphoto(i)) / C_cone;
        dV_dt = (-IKCa_cone(i)-ICa_cone(i)-ICl_cone(i) -IKv_cone(i) - IL_cone(i) - Idark_cone - Iphoto(i)) / C_cone;

        V(i+1) = V(i) + dt * dV_dt;

        % Update channel gating variables

        %mKv(i+1) = mKv(i) + dt * (amKv * (1 - mKv(i)) - bmKv * mKv(i));

        tau_mKv = 1 / (amKv + bmKv);
        mKv_inf = amKv * tau_mKv;
        mKv(i+1) = mKv_inf + (mKv(i) - mKv_inf) * exp(-dt / tau_mKv);

        %hKv(i+1) = hKv(i) + dt * (ahKv * (1 - hKv(i)) - bhKv * hKv(i));

        tau_hKv = 1 / (ahKv + bhKv);
        hKv_inf = ahKv * tau_hKv;
        hKv(i+1) = hKv_inf + (hKv(i) - hKv_inf) * exp(-dt / tau_hKv);

        %nh(i+1) = nh(i) + dt * (anh * (1 - nh(i)) - bnh * nh(i));

        tau_nh = 1 / (anh + bnh);
        nh_inf = anh * tau_nh;
        nh(i+1) = nh_inf + (nh(i) - nh_inf) * exp(-dt / tau_nh);

        %nCa(i+1) = nCa(i) + dt * (anCa * (1 - nCa(i)) - bnCa * nCa(i));

        tau_nCa = 1 / (anCa + bnCa);
        nCa_inf = anCa * tau_nCa;
        nCa(i+1) = nCa_inf + (nCa(i) - nCa_inf) * exp(-dt / tau_nCa);


        % Update calcium concentration
        dCai_dt = -ICa_cone(i) / (volume_cone * z * F * fI * fF) - (Cai(i) - Cainf_cone) / tau_Cai_cone;
        Cai(i+1) = Cai(i) + dt * dCai_dt;

        % disp(['step =' num2str(i)]);
        % disp(['mKv =' num2str(mKv(i))]);
        % disp(['amKv =' num2str(amKv)]);
        % disp(['bmKv =' num2str(bmKv)]);
        % disp(['ahKv =' num2str(ahKv)]);
        % disp(['bhKv =' num2str(bhKv)]);
        % disp(['mKv =' num2str(mKv(i+1))]);
        % disp(' ');
    end
    
    % Store the results 
    V_all(intensity_index, :) = V;
    Ih_all(intensity_index, :) = Ih_cone;
    ICa_all(intensity_index, :) = ICa_cone;
    IKv_all(intensity_index, :) = IKv_cone;
    ICl_all(intensity_index, :) = ICl_cone;
    IL_all(intensity_index, :) = IL_cone;
    IKCa_all(intensity_index, :) = IKCa_cone;
    Cai_all(intensity_index, :) = Cai;

    nh_all(intensity_index, :) = nh;
    mKv_all(intensity_index, :) = mKv;
    hKv_all(intensity_index, :) = hKv;
    nCa_all(intensity_index, :) = nCa;
    Cai_all(intensity_index, :) = Cai;

    
end

% Preallocate the total current matrix
Itotal_all = zeros(num_intensities, num_time_steps);

for i = 1:num_intensities
    Itotal_all(i, :) = IKv_all(i, :) + ...
                       Ih_all(i, :) + ...
                       ICa_all(i, :) + ...
                       IKCa_all(i, :) + ...
                       ICl_all(i, :) + ...
                       IL_all(i, :) + ...
                       Idark_cone + ...
                       Iphoto_all(i, :);
end

% Preallocate an array to store the lowest peak value for each intensity
lowest_peaks = zeros(num_intensities, 1);

% Loop through each intensity and find the minimum value in V_all
for i = 1:num_intensities
    % Ensure V_all contains valid values and matches the dimensions
    lowest_peaks(i) = min(V_all(i, :));
end

% Create strings for the legend that show both intensity and lowest peak values
legend_entries = cell(num_intensities, 1);
for i = 1:num_intensities
    legend_entries{i} = sprintf('R*=%g, Peak=%.2f mV', intensities(i), lowest_peaks(i));
end

% Plot the photovoltage
figure;
hold on;
for intensity_index = 1:num_intensities
    % Ensure V_all and time are correctly aligned
    plot(time, V_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', legend_entries{intensity_index});
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (mV)', 'FontSize', 14, 'FontWeight', 'bold');
title('Cone Photovoltage', 'FontSize', 16, 'FontWeight', 'bold');
%legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'off');
grid on;
xlim([0, t_end]); % Adjust x-axis range based on t_del
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot the total current vs. time for each intensity
figure('Name', 'Cone Total Current vs Time', 'NumberTitle', 'off');
hold on;
for i = 1:num_intensities
    plot(time, Itotal_all(i,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%g', intensities(i)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('Cone Total Current vs. Time', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([t_del * 0.75, t_del * 1.75]);  % Adjust as needed
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot the photocurrent for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Iphoto_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('Cone Photocurrent', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Photocurrent (pA)', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot Ih currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Ih_all(intensity_index, :),'LineWidth', 2.5, ... 
         'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Ih Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('Rod Ih Current', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del * 0.75, t_del * 1.75]);
grid on;
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot IKv currents for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKv_all(intensity_index, :), ...
         'LineWidth', 2.5, ... 
         'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('IKv Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('IKv Current', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del * 0.75, t_del * 1.75]);
grid on;
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot ICa currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICa_all(intensity_index, :), ...
         'LineWidth', 2.5, ... 
         'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('ICa Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('ICa Current', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del * 0.75, t_del * 1.75]);
grid on;
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot the IKCa current for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKCa_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('IKCa Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('IKCa Current', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot the Cl current (ICl) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICl_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Cl Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('ICl Current', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot the leak current (IL) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IL_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Leak Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('IL Current', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Plot gating variables
figure('Name', 'Gating Variables', 'NumberTitle', 'off');

% Kv Activation Variable (mKv)
subplot(2, 2, 1);
hold on;
for intensity_index = 1:num_intensities
    plot(time, mKv_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('mKv', 'FontSize', 14, 'FontWeight', 'bold');
title('Kv Activation Variable (mKv)', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

% Kv Inactivation Variable (hKv)
subplot(2, 2, 2);
hold on;
for intensity_index = 1:num_intensities
    plot(time, hKv_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('hKv', 'FontSize', 14, 'FontWeight', 'bold');
title('Kv Inactivation Variable (hKv)', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

% Ih Activation Variable (nh)
subplot(2, 2, 3);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nh_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('nh', 'FontSize', 14, 'FontWeight', 'bold');
title('Ih Activation Variable (nh)', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

% Calcium Activation Variable (nCa)
subplot(2, 2, 4);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nCa_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('nCa', 'FontSize', 14, 'FontWeight', 'bold');
title('Calcium Activation Variable (nCa)', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location','bestoutside', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([t_del * 0.75 t_del * 1.75]);
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;



%% Plotting Intensity Response for Photocurrent
function plot_intensity_response(Iamp_max, R_half, n)
    % Generate intensities and compute response
    intensities = 0:0.04:50000;
    Iamp_all = Iamp_max * (intensities.^n ./ (intensities.^n + R_half^n));
    
    % Experimental data
    dataX = [0.04, 0.4, 4, 40, 400, 4000, 40000];
    dataY = [0.0, 0.002, 0.01, 0.11, 0.59, 0.98, 1];

    % Plot intensity response
    figure('Name', 'Intensity Response', 'NumberTitle', 'off');
    semilogx(intensities, Iamp_all / Iamp_max, 'LineWidth', 2.5, 'DisplayName', 'Simulation');
    hold on;
    plot(dataX, dataY, 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
        'DisplayName', 'Experiment [Nange Jin, 2024]');
    xlabel('Flash Intensity (R*/cone/flash)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Normalized Current Response (Iamp/Iamp_{max})', 'FontSize', 14, 'FontWeight', 'bold');
    title('Photocurrent Intensity Response', 'FontSize', 16, 'FontWeight', 'bold');
    legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
    grid on;
    set(gca, 'FontWeight', 'bold', 'FontSize', 14);
    hold off;
end


plot_intensity_response(Iamp_max, R_half, n);

%% Plotting Intensity Response for Voltage
function plot_voltage_intensity_response(Vamp_max, R_half, n)
    % Generate intensities and compute voltage response
    intensities = 0:0.04:50000;
    Vamp_all = Vamp_max * (intensities.^n ./ (intensities.^n + R_half^n));

    % Experimental data
    intensity_levels = [0.04, 0.4, 4, 40, 400, 4000, 40000];
    voltage_data = [0.0, 0.0089, 0.01, 0.07, 0.42, 0.88, 1];

    % Plot voltage intensity response
    figure('Name', 'Voltage Intensity Response', 'NumberTitle', 'off');
    semilogx(intensities, Vamp_all / Vamp_max, 'LineWidth', 2.5, 'DisplayName', 'Simulation');
    hold on;
    plot(intensity_levels, voltage_data / max(voltage_data), 'o', ...
        'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Experiment [Nange Jin, 2024]');
    xlabel('Flash Intensity (R*/cone/flash)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Normalized Voltage Response (Vamp/Vamp_{max})', 'FontSize', 14, 'FontWeight', 'bold');
    title('Voltage Intensity Response', 'FontSize', 16, 'FontWeight', 'bold');
    legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
    grid on;
    set(gca, 'FontWeight', 'bold', 'FontSize', 14);
    hold off;
end

% Running the Voltage Intensity Plot Function with Updated Intensity Data
plot_voltage_intensity_response(Vamp_max, R_half,n);


%% Load data from Excel file for cone data
filename = '/Users/harshith/Downloads/rod-Cx36KO rod and cone mV pA from 19n05007-8 29-30 _1K_ 9632ms.xlsx';  
voltage_data = readmatrix(filename, 'Sheet', 1); % Specify the sheet for cone voltage data
current_data = readmatrix(filename, 'Sheet', 2); % Specify the sheet for cone current data

% Define intensity levels for labeling
intensity_levels = [0.04, 0.4, 4, 40, 400, 4000, 40000];  % Provided light intensity levels

% Initialize arrays to store peak values
voltage_peaks = zeros(length(intensity_levels), 1);
current_peaks = zeros(length(intensity_levels), 1);

%% Plot Experimental Photovoltage 
figure('Name', 'Cone Photovoltage (Experimental)', 'NumberTitle', 'off');
hold on;

% Cell array for legend entries
legend_entries = cell(length(intensity_levels), 1);

for i = 1:length(intensity_levels)
    voltage_trace = voltage_data(i, :);  % Each row corresponds to a specific intensity

    % Plot each voltage trace
    plot(voltage_trace, 'LineWidth', 2.5);

    % Find the peak (minimum) voltage
    voltage_peaks(i) = min(voltage_trace);  % For minimum peak

    % Create legend entry with intensity and peak value
    legend_entries{i} = sprintf('Intensity %.2f R*/cone/20ms: Peak = %.2f mV', intensity_levels(i), voltage_peaks(i));
end

xlabel('Time (arbitrary units)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (mV)', 'FontSize', 14, 'FontWeight', 'bold');
title(' Experimental Photovoltage', 'FontSize', 16, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on'); % Improved legend placement
xlim([0 5000]); % Adjust x-axis for better focus
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
grid on;
hold off;



%% Plot Experimental PhotoCurrent 
figure('Name', 'Cone Photocurrent (Experimental)', 'NumberTitle', 'off');
hold on;

% Cell array for legend entries
legend_entries = cell(length(intensity_levels), 1);

for i = 1:length(intensity_levels)
    current_trace = current_data(i, :);  % Each row corresponds to a specific intensity

    % Plot each current trace
    plot(current_trace, 'LineWidth', 2.5);

    % Find the peak (maximum) current
    current_peaks(i) = max(current_trace);  % For maximum peak

    % Create legend entry with intensity and peak value
    legend_entries{i} = sprintf('Intensity %.2f R*/cone/20ms: Peak = %.2f pA', intensity_levels(i), current_peaks(i));
end

xlabel('Time (mS)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('Experimental Photocurrent', 'FontSize', 16, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on'); 
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
xlim([0 5000]); % Adjust x-axis for better focus
hold off;



%% Interactive Plot for Simulated vs Experimental Photovoltage
V_all = V_all(:, 1:length(time));

% Create a figure for the interactive plot
figure('Name', 'Interactive Plot: Simulated vs Experimental Voltage', 'NumberTitle', 'off');
hold on;

% Define distinct colors for simulated and experimental plots
simulated_colors = lines(num_intensities); % Colors for simulated data

% Plot the simulated data (initially visible)
simulated_plots = gobjects(num_intensities, 1);
for intensity_index = 1:num_intensities
    simulated_plots(intensity_index) = plot(time, V_all(intensity_index, :), ...
        'Color', simulated_colors(intensity_index, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Simulated R*=%d', intensities(intensity_index)));
end

% Plot the experimental data with lighter colors
experimental_plots = gobjects(length(intensity_levels), 1);
for i = 1:length(intensity_levels)
    voltage_trace = voltage_data(i, :);
    lighter_color = simulated_colors(i, :);
    experimental_plots(i) = plot(1:length(voltage_trace), voltage_trace, ...
        '-', 'Color', lighter_color, ...
        'LineWidth', 2, ... 
        'DisplayName', sprintf('Experimental Intensity %.2f R*/cone/20ms', intensity_levels(i)));
end

% Add labels, title, and legend
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Combined Plot: Simulated vs Experimental Voltage');
legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on'); 
xlim([t_del * 0.75, t_del * 1.75]); % Adjust x-axis limits for better focus
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
grid on;

% Add toggle buttons to switch between simulated, experimental, or both
uicontrol('Style', 'pushbutton', 'String', 'Show Simulated', ...
    'Position', [20 20 120 40], ...
    'Callback', @(src, event) toggle_visibility('simulated', simulated_plots, experimental_plots));

uicontrol('Style', 'pushbutton', 'String', 'Show Experimental', ...
    'Position', [160 20 120 40], ...
    'Callback', @(src, event) toggle_visibility('experimental', simulated_plots, experimental_plots));

uicontrol('Style', 'pushbutton', 'String', 'Show Both', ...
    'Position', [300 20 120 40], ...
    'Callback', @(src, event) toggle_visibility('both', simulated_plots, experimental_plots));

hold off;



%% Load experimental photocurrent data from Excel file
experimental_photocurrent = readmatrix(filename, 'Sheet', 2); % Replace 'Sheet' with the actual sheet name

% Define intensity levels for labeling
intensity_levels = [0.04, 0.4, 4, 40, 400, 4000, 40000]; % Provided light intensity levels

Iphoto_all = Iphoto_all(:, 1:length(time)); % Trim simulated photocurrent data to match `time`

%% Create Interactive Plot for Simulated vs Experimental Photocurrent
figure('Name', 'Interactive Plot: Simulated vs Experimental Photocurrent', 'NumberTitle', 'off');
hold on;

% Define distinct colors for simulated and experimental plots
simulated_colors = lines(num_intensities); % Generate colors for simulated data

% Plot the simulated photocurrent data
simulated_plots = gobjects(num_intensities, 1); % Preallocate plot objects for simulated data
for intensity_index = 1:num_intensities
    simulated_plots(intensity_index) = plot(time, Iphoto_all(intensity_index, :), ...
        'Color', simulated_colors(intensity_index, :), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('Simulated R*=%d', intensities(intensity_index))); % Use appropriate intensity for labels
end

% Plot the experimental photocurrent data
experimental_plots = gobjects(length(intensity_levels), 1); % Preallocate plot objects for experimental data
for i = 1:length(intensity_levels)
    photocurrent_trace = experimental_photocurrent(i, :); % Extract experimental photocurrent trace
    lighter_color = simulated_colors(i, :); % Use matching color for corresponding intensity
    experimental_plots(i) = plot(1:length(photocurrent_trace), photocurrent_trace, ...
        '-', 'Color', lighter_color, ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('Experimental R*=%d', intensity_levels(i))); % Label with experimental intensities
end

% Add labels, title, and legend
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Photocurrent (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('Simulated vs Experimental Photocurrent', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'Box', 'on'); 
xlim([0,5000]);
grid on;

% Add toggle buttons to switch between simulated, experimental, or both
uicontrol('Style', 'pushbutton', 'String', 'Show Simulated', ...
    'Position', [20 20 120 40], ...
    'Callback', @(src, event) toggle_visibility('simulated', simulated_plots, experimental_plots));

uicontrol('Style', 'pushbutton', 'String', 'Show Experimental', ...
    'Position', [160 20 120 40], ...
    'Callback', @(src, event) toggle_visibility('experimental', simulated_plots, experimental_plots));

uicontrol('Style', 'pushbutton', 'String', 'Show Both', ...
    'Position', [300 20 120 40], ...
    'Callback', @(src, event) toggle_visibility('both', simulated_plots, experimental_plots));

hold off;

%% Helper Function for Toggling Visibility
function toggle_visibility(option, simulated_plots, experimental_plots)
    switch option
        case 'simulated'
            % Show only simulated plots
            for i = 1:length(simulated_plots)
                simulated_plots(i).Visible = 'on'; % Show simulated
            end
            for i = 1:length(experimental_plots)
                experimental_plots(i).Visible = 'off'; % Hide experimental
            end
        case 'experimental'
            % Show only experimental plots
            for i = 1:length(simulated_plots)
                simulated_plots(i).Visible = 'off'; % Hide simulated
            end
            for i = 1:length(experimental_plots)
                experimental_plots(i).Visible = 'on'; % Show experimental
            end
        case 'both'
            % Show both simulated and experimental plots
            for i = 1:length(simulated_plots)
                simulated_plots(i).Visible = 'on'; % Show simulated
            end
            for i = 1:length(experimental_plots)
                experimental_plots(i).Visible = 'on'; % Show experimental
            end
    end
end