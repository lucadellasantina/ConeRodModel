intensities = 50;

% Parameters of the rod cell (Liu & Kourennyi 2004)
C_rod = 10; % Capacitance in pF
ECa_rod = 40; % Calcium reversal potential in mV
EKv_rod = -74; % Kv reversal potential in mV
EKx_rod = -74; % Kx reversal potential in mV
Eh_rod = -32; % Ih reversal potential in mV
EKCa_rod = -74; % KCa reversal potential in mV
ECl_rod = -20; % Cl reversal potential in mV
EL_rod = -74; % Leak reversal potential in mV

gCa_rod = 4; % Calcium conductance in nS
gKv_rod = 10; % Kv conductance in nS
gKx_rod = 1.04; % Kx conductance in nS
gh_rod = 2.5; % Ih conductance in nS
gKCa_rod = 5; % KCa conductance in nS
gCl_rod = 1.3; % Cl conductance in nS
gL_rod = 0.52; % Leak conductance in n

K_KCa_rod = 0.00032; % KCa constant in mM
K_Cl_rod = 0.0015; % Cl constant in mM
Cainf_rod = 5e-5; % Steady-state calcium concentration in mM

% Intracellular calcium concentration parameters
radius_rod = 7.8; % microns
volume_rod = (4/3) * pi * radius_rod^3; % volume in um^3
z = 2; % Valence of calcium ions
tau_Cai_rod = 20; % Time constant for calcium decay in ms

% Conversion factors
fI = 10e-6; % Conversion factor for ICa
fF = 10e9; % Conversion factor for Faraday's constant
F = 96485; % Faraday constant (C/mol)

% Rod photocurrent parameters
Idark_rod = -42; % Dark current in pA

% Photocurrent model parameters
Iamp_max = 30; % Maximal amplitude of photocurrent in pA
Vamp_max = -42;
n = 1; % Hill coefficient
R_half = 4; % Flash intensity at half-maximal amplitude
R_sat = 1000; % Saturating flash intensity

% Time parameters
% TODO: Investigate why model breaks if dt>0.1
dt = 0.1; % Time step in ms
t_end = 1500; % Total simulation time in ms
num_time_steps = t_end / dt;  

% Generate the time vector with exactly 50,000 points from 0 to t_end - dt
time = linspace(0, t_end, num_time_steps);
t_del =500;

% Intensities to simulate
num_intensities = length(intensities);


% Preallocate variables for currents and states
V_all = zeros(num_intensities, num_time_steps);
Iphoto_all = zeros(num_intensities, num_time_steps);
IKv_all = zeros(num_intensities, num_time_steps);
IKx_all = zeros(num_intensities, num_time_steps);
Ih_all = zeros(num_intensities, num_time_steps);
ICa_all = zeros(num_intensities, num_time_steps);
IKCa_all = zeros(num_intensities, num_time_steps);
ICl_all = zeros(num_intensities, num_time_steps);
IL_all = zeros(num_intensities, num_time_steps);
I_noise = zeros(1, num_time_steps);
Cai_all = zeros(num_intensities, num_time_steps);

nKv_all = zeros(num_intensities, num_time_steps);
nKx_all = zeros(num_intensities, num_time_steps);
nh_all = zeros(num_intensities, num_time_steps);
mCa_all = zeros(num_intensities, num_time_steps);
hCa_all = zeros(num_intensities, num_time_steps);

% Function to calculate photocurrent Iphoto for rod
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
    R_star = intensities(intensity_index);
    
    % Initial membrane parameters
    V = -42.8 * ones(1, num_time_steps); % Voltage for each intensity  (mV)
    nKv = 0.1518 * ones(1, num_time_steps); % Kv activation variable
    nKx = 0.7735 * ones(1, num_time_steps); % Kx inactivation variable
    nh = 0.0007 * ones(1, num_time_steps); % Ih activation variable
    mCa = 0.0041 * ones(1, num_time_steps); % Calcium activation variable
    hCa = 0.9975 * ones(1, num_time_steps); % Calcium inactivation variable
    Cai = 0.00000001 * ones(1, num_time_steps); % Initial intracellular calcium concentration in mM
    
    % Preallocate currents
    Ih_rod = zeros(1, num_time_steps); % Reset Ih current for each intensity
    ICa_rod = zeros(1, num_time_steps); % Reset Calcium current for each intensity
    IKv_rod = zeros(1, num_time_steps); % Reset Kv current for each intensity
    IKx_rod =  zeros(1, num_time_steps); % Reset Kx current for each intensity
    IKCa_rod = zeros(1, num_time_steps); % Reset IKCa current for each intensity
    ICl_rod = zeros(1, num_time_steps); % Reset Cl current for each intensity
    IL_rod = zeros(1, num_time_steps); % Reset Leak current for each intensity

    % Hill equation for Iamp
    R_star = max(R_star, 1e-6); % Prevent log(0)
    Iamp = Iamp_max * (R_star^n / (R_star^n + R_half^n));
    
    % Time constants and c based on R_star
    if R_star <= 100000
        c = (366.05 * log(R_star / R_half) + 1347)/3;
        tau1 = (-5.263 * log(R_star / R_half) + 55.235)/3;
        tau2 = (20.392 * log(R_star / R_half) + 236.09)/3;
    end

    % Compute Iphoto for all time steps
    Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c, t_del, num_time_steps, dt);

    % Add Iphoto to total photocurrent at each time step
    Iphoto_all(intensity_index, :) = Iphoto;

    % Main simulation loop
    for i = 1:num_time_steps-1
        % Update gating variables
        anKv = 0.005 * (20 - V(i)) / (exp((20 - V(i)) / 22) - 1);
        bnKv = 0.0625 * exp(-V(i) / 80);

        anKx = 0.001 * 0.66 * exp((V(i) + 50) / (2 * 5.7));
        bnKx = 0.001 * 0.66 * exp(-(V(i) + 50) / (2 * 5.7));

        anh = 0.001 * exp((V(i) + 75) / (-10.6));
        bnh = 0.001 * exp(-(V(i) + 75) / (-10.6));

        amCa = 0.1 * exp((V(i) + 10) / 12);
        bmCa = 0.1 * exp(-(V(i) + 10) / 12);

        ahCa = 0.001 * 0.5 * exp(-(V(i) - 11) / 18);
        bhCa = 0.001 * exp((V(i) - 11) / 18);

        nKCa = 1 / (1 + (K_KCa_rod / Cai(i))^4);
        nCl = 1 / (1 + (K_Cl_rod / Cai(i))^4);

        % Update currents
        IKv_rod(i) = gKv_rod * nKv(i)^4 * (V(i) - EKv_rod); % IKv current in pA
        IKx_rod(i) = gKx_rod * nKx(i) * (V(i) - EKx_rod); % IKx current in pA
        ICa_rod(i) = gCa_rod * mCa(i) * hCa(i) * (V(i) - ECa_rod); % Calcium current in pA
        IKCa_rod(i) = gKCa_rod * nKCa^4 * (V(i) - EKCa_rod); % IKCa current in pA
        ICl_rod(i) = gCl_rod * nCl * (V(i) - ECl_rod); % ICl current in pA
        IL_rod(i) = gL_rod * (V(i) - EL_rod); % IL current in pA
        
        % Update Ih current
        Ih_rod(i) = gh_rod * nh(i) * (V(i) - Eh_rod);

        % Update membrane voltage
        % dV_dt = (-IKv_rod(i) ...
        %     - IKx_rod(i) ...
        %     - Ih_rod(i) ...
        %     - ICa_rod(i) ...
        %     - IKCa_rod(i) ...
        %     - ICl_rod(i) ...
        %     - IL_rod(i) ...
        %     - Idark_rod ...
        %     - Iphoto(i)) / C_rod;
        dV_dt = (-IKv_rod(i) -IKx_rod(i) -Ih_rod(i) -ICa_rod(i) -IKCa_rod(i) - ICl_rod(i) -IL_rod(i) - Idark_rod - Iphoto(i)) / C_rod;

        V(i+1) = V(i) + dt * dV_dt;

        % Update channel gating variables
        %nKv(i+1) = nKv(i) + dt * (anKv * (1 - nKv(i)) - bnKv * nKv(i));
        
        tau_nKv = 1 / (anKv + bnKv);
        nKv_inf = anKv * tau_nKv;
        nKv(i+1) = nKv_inf + (nKv(i) - nKv_inf) * exp(-dt / tau_nKv);

        %nKx(i+1) = nKx(i) + dt * (anKx * (1 - nKx(i)) - bnKx * nKx(i));

        tau_nKx = 1 / (anKx + bnKx);
        nKx_inf = anKx * tau_nKx;
        nKx(i+1) = nKx_inf + (nKx(i) - nKx_inf) * exp(-dt / tau_nKx);


        %nh(i+1) = nh(i) + dt * (anh * (1 - nh(i)) -  bnh * nh(i));

        tau_nh = 1 / (anh + bnh);
        nh_inf = anh * tau_nh;
        nh(i+1) = nh_inf + (nh(i) - nh_inf) * exp(-dt / tau_nh);

        %mCa(i+1) = mCa(i) + dt * (amCa * (1 - mCa(i)) - bmCa * mCa(i));

        tau_mCa = 1 / (amCa + bmCa);
        mCa_inf = amCa * tau_mCa;
        mCa(i+1) = mCa_inf + (mCa(i) - mCa_inf) * exp(-dt / tau_mCa);

        %hCa(i+1) = hCa(i) + dt * (ahCa * (1 - hCa(i)) - bhCa * hCa(i));

        tau_hCa = 1 / (ahCa + bhCa);
        hCa_inf = ahCa * tau_hCa;
        hCa(i+1) = hCa_inf + (hCa(i) - hCa_inf) * exp(-dt / tau_hCa);

        % Update Calcium concentration 
        dCai_dt = -ICa_rod(i) / (volume_rod * z * F * fI * fF) - (Cai(i) - Cainf_rod) / tau_Cai_rod;
        Cai(i+1) = Cai(i) + dt * dCai_dt;
    end

    % Store the results
    V_all(intensity_index, :) = V;
    Ih_all(intensity_index, :) = Ih_rod;
    ICa_all(intensity_index, :) = ICa_rod;
    IKv_all(intensity_index, :) = IKv_rod;
    IKx_all(intensity_index, :) = IKx_rod; 
    IKCa_all(intensity_index, :) = IKCa_rod;
    ICl_all(intensity_index, :) = ICl_rod;
    IL_all(intensity_index, :) = IL_rod;
    Cai_all(intensity_index, :) = Cai;

    nKv_all(intensity_index, :) = nKv;
    nKx_all(intensity_index, :) = nKx;
    nh_all(intensity_index, :) = nh;
    mCa_all(intensity_index, :) = mCa;
    hCa_all(intensity_index, :) = hCa;

    
end

% Preallocate an array to store the lowest peak value for each intensity
lowest_peaks = zeros(num_intensities, 1);

% Loop through each intensity to find the minimum peak voltage
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
    % Ensure alignment between time and V_all
    plot(time, V_all(intensity_index, :),'LineWidth', 2.5, ...
        'DisplayName', legend_entries{intensity_index}); 
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (mV)', 'FontSize', 14, 'FontWeight', 'bold');
title('Rod Photovoltage', 'FontSize', 16, 'FontWeight', 'bold');
%legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
xlim([0, t_end]); % Adjust x-axis range based on t_del
grid on;
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
hold off;

%% Compute total current (Itotal) and plot vs time
% Preallocate total current matrix
Itotal_all = zeros(num_intensities, num_time_steps);

% Sum all the individual ionic currents plus Idark_rod and Iphoto
for intensity_index = 1:num_intensities
    Itotal_all(intensity_index, :) = ...
        IKv_all(intensity_index, :) + ...
        IKx_all(intensity_index, :) + ...
        Ih_all(intensity_index, :) + ...
        ICa_all(intensity_index, :) + ...
        IKCa_all(intensity_index, :) + ...
        ICl_all(intensity_index, :) + ...
        IL_all(intensity_index, :) + ...
        Idark_rod + ...
        Iphoto_all(intensity_index, :);
    %Itotal_all(intensity_index, :) = Iphoto_all(intensity_index, :);

end

% Plot total current for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Itotal_all(intensity_index,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('R*=%g', intensities(intensity_index)));
end
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total Current (pA)', 'FontSize', 14, 'FontWeight', 'bold');
title('Rod Total Current vs. Time', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 10, 'Box', 'on');
xlim([0, t_del*1.75]);  % Adjust the x-axis for your region of interest
grid on;
hold off;

%% Plot the photocurrent for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Iphoto_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('Rod Photocurrent ','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Photocurrent (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
    set(gca, 'FontWeight', 'bold','FontSize', 14);
grid on;
hold off;

%% Plot Ih currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Ih_all(intensity_index, :),'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('Rod Ih Current ','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Ih (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del* 0.20, t_del*3.5]);
grid on;
set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

%% Plot IKv currents for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKv_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('IKv Current','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('IKv (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
    set(gca, 'FontWeight', 'bold','FontSize', 14);
grid on;
hold off;

%% Plot IKx currents for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKx_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('IKx Current ','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('IKx (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
    set(gca, 'FontWeight', 'bold','FontSize', 14);
grid on;
hold off;

%% Plot ICa currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICa_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('ICa Current','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('ICa (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
    set(gca, 'FontWeight', 'bold','FontSize', 14);
grid on;
hold off;

%% Plot IKCa currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKCa_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('IKCa Current ','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('IKCa (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
    set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

%% Plot ICl currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICl_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('ICl Current ','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('ICl (pA)','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
    set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

%% Plot the leak current (IL) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IL_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Leak Current close','FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
title('IL Current','FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

%% Plot gating variables for different intensities
figure;
subplot(2,2,1);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nKv_all(intensity_index, :),'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('nKv Gating Variable','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('nKv','FontSize', 14, 'FontWeight', 'bold');
xlim([t_del*0.20 t_del*3.5]);
set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

subplot(2,2,2);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nKx_all(intensity_index, :),'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('nKx Gating Variable','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('nKx','FontSize', 14, 'FontWeight', 'bold');
xlim([t_del*0.20 t_del*3.5]);
set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

subplot(2,2,3);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nh_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('nh Gating Variable','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('nh','FontSize', 14, 'FontWeight', 'bold');
xlim([t_del*0.20 t_del*3.5]);
set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

subplot(2,2,4);
hold on;
for intensity_index = 1:num_intensities
    plot(time, mCa_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('mCa Gating Variable','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('mCa','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
xlim([t_del*0.20 t_del*3.5]);
set(gca, 'FontWeight', 'bold','FontSize', 14);
hold off;

%% Plot for hCa
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, hCa_all(intensity_index, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('hCa Gating Variable for All Intensities','FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (ms)','FontSize', 14, 'FontWeight', 'bold');
ylabel('hCa','FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 10, 'Box', 'on');
    set(gca, 'FontWeight', 'bold','FontSize', 14);
xlim([t_del*0.20 t_del*3.5]);
hold off;


%% Plotting Intensity Response for Voltage 
function plot_intensity_response(Iamp_max, R_half, n)

    intensities = 0:0.04:10000;
    Iamp_all = Iamp_max * (intensities.^n ./ (intensities.^n + R_half^n));
    
% Updated intensity levels and averaged voltage data from Excel
dataX = [0.04, 0.4, 4, 40, 400, 4000];  % Using provided intensity levels
dataY = [0.0, 0.07, 0.59, 0.88, 0.96, 1];

    % Plot intensity response
    figure;
    semilogx(intensities, Iamp_all / Iamp_max,'LineWidth', 2.5, 'DisplayName', 'Simulation');
    hold on;
    plot(dataX, dataY, 'o', 'MarkerSize', 8, 'LineWidth', 2,'DisplayName', 'Experiment [Nange Jin, 2024]');
    xlabel('Flash intensity (in R*/rod/flash)','FontSize', 14, 'FontWeight', 'bold');
    ylabel('Normalized response (Iamp/Iamp_{max})','FontSize', 14, 'FontWeight', 'bold');
    legend( 'FontSize', 10, 'Box', 'on');
    title('Photocurrent Intensity Response','FontSize', 16, 'FontWeight', 'bold');
        set(gca, 'FontWeight', 'bold','FontSize', 14);
    hold off;
end

plot_intensity_response(Iamp_max, R_half, n);


%% Plotting Intensity Response for Voltage 
function plot_voltage_intensity_response(Vamp_max, R_half, n)

    % Define the range of intensities for simulation
    intensities = 0:0.04:10000;

    % Calculate the simulated voltage response curve
    Vamp_all = Vamp_max * (intensities.^n ./ (intensities.^n + R_half^n));

    % Updated intensity levels and averaged voltage data from Excel
    intensity_levels = [0.04, 0.4, 4, 40, 400, 4000];  % Using provided intensity levels
    voltage_data = [0.0, 0.05, 0.61, 0.92, 0.93, 1];

    % Plot intensity response on a semi-logarithmic scale
    figure;
    semilogx(intensities, Vamp_all / Vamp_max,'LineWidth', 2.5, 'DisplayName', 'Simulation');
    hold on;
    plot(intensity_levels, voltage_data / max(voltage_data), 'o','MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Experiment [Nange Jin, 2024]');
    xlabel('Flash intensity (in R*/rod/flash)','FontSize', 14, 'FontWeight', 'bold');
    ylabel('Normalized Voltage Response (Vamp/Vamp_{max})','FontSize', 14, 'FontWeight', 'bold');
    legend( 'FontSize', 10, 'Box', 'on');
    title('Voltage Intensity Response','FontSize', 16, 'FontWeight', 'bold');

    % Set plot properties for clarity
    set(gca, 'FontWeight', 'bold','FontSize', 14);
    hold off;
end

%Plot the voltage intensity response
plot_voltage_intensity_response(Vamp_max, R_half, n);



%% Load data from Excel file
filename = '/Users/harshith/downloads/rod-Cx36KO rod and cone mV pA from 19n05007-8 29-30 _1K_ 9632ms.xlsx';  % Replace with your actual file path
voltage_data = readmatrix(filename, 'Sheet', 3);
current_data = readmatrix(filename, 'Sheet', 4);

% Define intensity levels for labeling
intensity_levels = [0.04, 0.4, 4, 40, 400, 4000];  % Provided light intensity levels

%% Plot Experimental Photovoltage for Different Intensities 
figure;
hold on;
for i = 1:length(intensity_levels)
    voltage_trace = voltage_data(i, :);  % Each row corresponds to a specific intensity
    plot(voltage_trace, 'LineWidth', 2,'DisplayName', sprintf('Intensity: %.2f R*/rod/20ms', intensity_levels(i)));
end
xlabel('Time (mS)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage (mV)','FontSize', 14, 'FontWeight', 'bold');
title('Experimental Photovoltage','FontSize', 16, 'FontWeight', 'bold');
legend('show', 'FontSize', 12, 'Box', 'off');
grid on;
xlim([0 5000]);
hold off;

%% Plot Experimental Photocurrent Different Intensities 
figure;
hold on;
for i = 1:length(intensity_levels)
    current_trace = current_data(i, :);  % Each row corresponds to a specific intensity
    plot(current_trace, 'LineWidth', 2,'DisplayName', sprintf('Intensity: %.2f R*/rod/20ms', intensity_levels(i)));
end
xlabel('Time (mS)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Current (pA)','FontSize', 14, 'FontWeight', 'bold');
title('Experimental photocurrent','FontSize', 14, 'FontWeight', 'bold');
legend('show', 'FontSize', 12, 'Box', 'off');
grid on;
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
        'LineWidth', 2, ...
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
legend('show', 'Location', 'bestoutside'); % Place legend outside the plot for clarity
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

%% Interactive Plot for Simulated vs Experimental Photocurrent

experimental_photocurrent = readmatrix(filename, 'Sheet', 4); % Replace 'Sheet' with the actual sheet name

% Define intensity levels for labeling
intensity_levels = [0.04, 0.4, 4, 40, 400, 4000]; % Provided light intensity levels

% Ensure `time` matches simulation dimensions
Iphoto_all = Iphoto_all(:, 1:length(time)); % Trim simulated photocurrent data to match `time`

% Ensure there are enough colors for both simulated and experimental data
max_colors = max(num_intensities, length(intensity_levels));
simulated_colors = lines(max_colors); % Generate enough colors for all plots

%% Create Interactive Plot for Simulated vs Experimental Photocurrent
figure('Name', 'Interactive Plot: Simulated vs Experimental Photocurrent', 'NumberTitle', 'off');
hold on;

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
legend('show', 'Location', 'bestoutside', 'FontSize', 12, 'Box', 'off'); % Place legend outside the plot for clarity
xlim([0, 5000]);
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