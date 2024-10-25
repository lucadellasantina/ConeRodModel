% Parameters of the rod cell (Liu & Kourennyi 2004)
C_rod = 30; % Capacitance in pF
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
gL_rod = 0.52; % Leak conductance in nS

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

% Noise parameters
sigma_noise_rod = 0.5; % Noise amplitude in pA
tau_noise_rod = 1; % Noise time constant in ms

% Rod photocurrent parameters
Idark_rod = -40; % Dark current in pA

% Photocurrent model parameters
Iamp_max = 16; % Maximal amplitude of photocurrent in pA
n = 1; % Hill coefficient
R_half = 20; % Flash intensity at half-maximal amplitude
R_sat = 200; % Saturating flash intensity

%% Time parameters
dt = 0.1; % Time step in ms
t_end = 4000; % Total simulation time in ms
num_time_steps = t_end / dt;  
% Generate the time vector with exactly 40,000 points from 0 to t_end - dt
time = linspace(0, t_end - dt, num_time_steps);  % Exclude the final step at t_end
t_del =1000;

% Intensities to simulate
num_intensities = length(intensities);

%% Preallocate variables for currents and states
V_all = zeros(num_intensities, num_time_steps);
Iphoto_all = zeros(num_intensities, num_time_steps);
IKv_all = zeros(num_intensities, num_time_steps);
IKx_all = zeros(num_intensities, num_time_steps);
Ih_all = zeros(num_intensities, num_time_steps);
ICa_all = zeros(num_intensities, num_time_steps);
IKCa_all = zeros(num_intensities, num_time_steps);
ICl_all = zeros(num_intensities, num_time_steps);
IL_all = zeros(num_intensities, num_time_steps);
Cai_all = zeros(num_intensities, num_time_steps);

nKv_all = zeros(num_intensities, num_time_steps);
nKx_all = zeros(num_intensities, num_time_steps);
nh_all = zeros(num_intensities, num_time_steps);
mCa_all = zeros(num_intensities, num_time_steps);
hCa_all = zeros(num_intensities, num_time_steps);

%% Function to calculate photocurrent Iphoto for rod
function Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c, t_del)
    Iphoto = Iamp .* ((1 - exp(-(time - t_del) / tau1))) .* (time > t_del) ...
           - Iamp ./ (1 + exp(-(time - t_del - c) / tau2)) .* (time > t_del);
end

%% Main simulation loop over each intensity
for intensity_index = 1:num_intensities
    R_star = intensities(intensity_index);
    
    % Initial membrane parameters
    V = -42.8 * ones(1, num_time_steps); % Initial voltage (mV)
    nKv = 0.1518 * ones(1, num_time_steps);
    nKx = 0.7735 * ones(1, num_time_steps);
    nh = 0.0007 * ones(1, num_time_steps);
    mCa = 0.0041 * ones(1, num_time_steps);
    hCa = 0.9975 * ones(1, num_time_steps);
    Cai = 0.00001 * ones(1, num_time_steps); % Initial calcium concentration in mM
    
    % Preallocate currents
    Ih_rod = zeros(1, num_time_steps);
    ICa_rod = zeros(1, num_time_steps);
    IKv_rod = zeros(1, num_time_steps);
    IKx_rod =  zeros(1, num_time_steps);
    IKCa_rod = zeros(1, num_time_steps);
    ICl_rod = zeros(1, num_time_steps);
    IL_rod = zeros(1, num_time_steps);
    I_noise = zeros(1, num_time_steps);

    % Hill equation for Iamp
    R_star = max(R_star, 1e-6); % Prevent log(0)
    Iamp = Iamp_max * (R_star^n / (R_star^n + R_half^n));
    
    % Time constants and c based on R_star
    if R_star <= R_sat
        c = (66.05 * log(R_star / R_half) + 647);
        tau1 = -5.263 * log(R_star / R_half) + 34.235;
        tau2 = 20.392 * log(R_star / R_half) + 106.09;
    else
        c = 0.33 * R_star + 836;
        tau1 = -5.263 * log(R_star / R_half) + 34.235;
        tau2 = 20.392 * log(R_star / R_half) + 106.09;
    end

    % Compute Iphoto for all time steps
    Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c, t_del);
    Iphoto_all(intensity_index, :) = Iphoto;

    %% Time loop for the simulation
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
        IKv_rod(i) = gKv_rod * nKv(i)^4 * (V(i) - EKv_rod);
        IKx_rod(i) = gKx_rod * nKx(i) * (V(i) - EKx_rod);
        ICa_rod(i) = gCa_rod * mCa(i) * hCa(i) * (V(i) - ECa_rod);
        IKCa_rod(i) = gKCa_rod * nKCa^4 * (V(i) - EKCa_rod);
        ICl_rod(i) = gCl_rod * nCl * (V(i) - ECl_rod);
        IL_rod(i) = gL_rod * (V(i) - EL_rod);
        Ih_rod(i) = gh_rod * nh(i) * (V(i) - Eh_rod);

        % Noise generation
        noise = sigma_noise_rod * sqrt(2 / tau_noise_rod) * randn() / sqrt(dt);
        I_noise(i+1) = I_noise(i) + dt * (-I_noise(i) / tau_noise_rod + noise);

        % Membrane voltage update
        dV_dt = (-IKv_rod(i) - IKx_rod(i) - Ih_rod(i) - ICa_rod(i) - IKCa_rod(i) ...
                - ICl_rod(i) - IL_rod(i) - Idark_rod - Iphoto(i) - I_noise(i)) / C_rod;
        V(i+1) = V(i) + dt * dV_dt;

        % Update channel gating variables
        nKv(i+1) = nKv(i) + dt * (anKv * (1 - nKv(i)) - bnKv * nKv(i));
        nKx(i+1) = nKx(i) + dt * (anKx * (1 - nKx(i)) - bnKx * nKx(i));
        nh(i+1) = nh(i) + dt * (anh * (1 - nh(i)) - bnh * nh(i));
        mCa(i+1) = mCa(i) + dt * (amCa * (1 - mCa(i)) - bmCa * mCa(i));
        hCa(i+1) = hCa(i) + dt * (ahCa * (1 - hCa(i)) - bhCa * hCa(i));

        % Calcium concentration update
        dCai_dt = -ICa_rod(i) / (volume_rod * z * F * fI * fF) - (Cai(i) - Cainf_rod) / tau_Cai_rod;
        Cai(i+1) = Cai(i) + dt * dCai_dt;
    end

    % Store results
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


%% Preallocate an array to store the lowest peak value for each intensity
lowest_peaks = zeros(length(intensities), 1);

colors = lines(num_intensities); % You can use other colormaps like 'jet', 'hsv', etc.


% Loop through each intensity and find the minimum value in V_all
for i = 1:length(intensities)
    % Find the minimum (lowest peak) value for this intensity curve
    lowest_peaks(i) = min(V_all(i, :));
end

% Create strings for the legend that show both intensity and lowest peak values
legend_entries = cell(length(intensities), 1);
for i = 1:length(intensities)
    legend_entries{i} = sprintf('R*=%d, Peak=%.2f mV', intensities(i), lowest_peaks(i));
end

% Plot the photovoltage
figure;
hold on;

% Initialize an array to hold plot handles for the legend
plot_handles = [];
for idx = 1:length(intensities)
    plot(time, V_all(idx, :)); % Plot each curve without a specific legend entry
end

% Add a legend that includes both intensity and lowest peak values
lgd = legend(plot_handles, legend_entries, 'Location', 'best');
set(lgd, 'Units', 'normalized'); % Allows the legend to be repositioned

% Adjust plot labels and limits
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Cone Photovoltage');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot Ih currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Ih_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('Ih Current (pA)');
xlabel('Time (ms)');
ylabel('Ih (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot photocurrents for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Iphoto_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('Rod Photocurrent (pA)');
xlabel('Time (ms)');
ylabel('Photocurrent (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot IKv currents for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKv_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('IKv Current (pA)');
xlabel('Time (ms)');
ylabel('IKv (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot IKx currents for all intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKx_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('IKx Current (pA)');
xlabel('Time (ms)');
ylabel('IKx (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;



% Plot ICa currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICa_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('ICa Current (pA)');
xlabel('Time (ms)');
ylabel('ICa (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot IKCa currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKCa_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('IKCa Current (pA)');
xlabel('Time (ms)');
ylabel('IKCa (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot ICl currents for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICl_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('ICl Current (pA)');
xlabel('Time (ms)');
ylabel('ICl (pA)');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Plot the leak current (IL) for different intensities
figure;
hold on;
for idx = 1:num_intensities
    plot(time, IL_all(idx, :), 'DisplayName', sprintf('R*=%d', intensities(idx)));
end
xlabel('Time (ms)');
ylabel('Leak Current (pA)');
legend('show', 'Location', 'best');
xlim([t_del*0.20 t_del*3.5]);
title('IL Current');
hold off;


% Plot gating variables for different intensities
figure;
subplot(2,2,1);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nKv_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('nKv Gating Variable');
xlabel('Time (ms)');
ylabel('nKv');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

subplot(2,2,2);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nKx_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('nKx Gating Variable');
xlabel('Time (ms)');
ylabel('nKx');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

subplot(2,2,3);
hold on;
for intensity_index = 1:num_intensities
    plot(time, nh_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('nh Gating Variable');
xlabel('Time (ms)');
ylabel('nh');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

subplot(2,2,4);
hold on;
for intensity_index = 1:num_intensities
    plot(time, mCa_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('mCa Gating Variable');
xlabel('Time (ms)');
ylabel('mCa');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;

% Additional plot for hCa
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, hCa_all(intensity_index, :), 'Color', colors(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
title('hCa Gating Variable for All Intensities');
xlabel('Time (ms)');
ylabel('hCa');
legend('Location', 'eastoutside');
xlim([t_del*0.20 t_del*3.5]);
hold off;


%% Plotting Intensity Response
function plot_intensity_response(Iamp_max, R_half, n)

    intensities = 0:0.1:10000;
    Iamp_all = Iamp_max * (intensities.^n ./ (intensities.^n + R_half^n));
    
    % Experimental data from Ingram et al. (2017)
    dataX = [1, 2, 5, 10, 30, 100, 200, 500, 1000, 2000];
    dataY = [0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 0.95, 1, 1];

    % Plot intensity response
    figure;
    semilogx(intensities, Iamp_all / Iamp_max, 'DisplayName', 'Simulation');
    hold on;
    plot(dataX, dataY, 'o', 'DisplayName', 'Experiment [Ingram et al., 2017]');
    xlabel('Flash intensity (in R*/rod/flash)');
    ylabel('Normalized response (Iamp/Iamp_{max})');
    legend();
    title('Intensity Response');
    hold off;
end

%% Run the simulations
intensities = [1, 2, 5, 10, 30, 100, 200, 500, 1000, 2000];
time = 0:dt:t_end;

plot_intensity_response(Iamp_max, R_half, n);
