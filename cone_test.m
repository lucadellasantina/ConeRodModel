intensities = [1, 2, 3, 5, 10,20, 30, 50,100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000]; % Array of intensities

% Parameters of the cone cell (Kourennyi 2004)
C_cone = 16; % Capacitance in pF
ECa_cone = 40; % Calcium reversal potential in mV
EKv_cone = -80; % Kv reversal potential in mV
Eh_cone = -32.5; % Ih reversal potential in mV
EKCa_cone = -80; % KCa reversal potential in mV
ECl_cone = -45; % Cl reversal potential in mV
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
Cainf_cone = 0.0001; % Steady-state calcium concentration in mM


% Intracellular calcium concentration (Liu & Kourennyi 2004)
radius_cone = 7.8; % microns
volume_cone = (4/3) * pi * radius_cone^3; % volume in um^3
z = 2; % Valence of calcium ions
tau_Cai_cone = 20; % Time constant for calcium decay in ms

% Conversion factors
fI = 10e-6; % Conversion factor for ICa (already in pA)
fF = 10e9; % Conversion factor for Faraday's constant
F = 96485; % Faraday constant (C/mol)

% Noise parameters
sigma_noise_cone = 0.5; % Noise amplitude in pA
tau_noise_cone = 1; % Noise time constant in ms

% Cone photocurrent parameters
Idark_cone = -20; % Dark current in pA

% Photocurrent model parameters
Iamp_max = 30; % pA maximal amplitude of photocurrent
n = 1;         % Hill coefficient
R_half = 950;   % flash intensity at half-maximal amplitude
R_sat = 2;   % saturating flash intensity


%% Time parameters
dt = 0.01; % Time step in ms (same as 1e-5 s)
T = 4000; % Total simulation time in ms
time = 0:dt:T-dt; % Time vector
t_del = 1000; % stimulus trigered offset time (ms)
num_time_steps = length(time);


%% Preallocate variables for currents and states
num_intensities = length(intensities);
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



%% Function to calculate photocurrent Iphoto for cone
function Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c, t_del)
    Iphoto = Iamp .* ((1 - exp(-(time - t_del) / tau1))).* (time > t_del)  - Iamp .* ((1 ./ (1 + exp(-(time - t_del - c) / tau2)))) .* (time > t_del);
end

%% Main simulation loop
for intensity_index = 1:num_intensities
    R_star = intensities(intensity_index) ; % Current intensity
    

    % Preallocate variables 
    V = -41 * ones(1, num_time_steps); % Reset voltage for each intensity  (mV)
    mKv = 0.3709 * ones(1, num_time_steps); % Kv activation variable
    hKv = 0.9998 * ones(1, num_time_steps); % Kv inactivation variable
    nCa = 0.0115 * ones(1, num_time_steps); % Calcium activation variable
    nh = 0.0877 * ones(1, num_time_steps); % Ih activation variable
    Cai = 0.00001 * ones(1, num_time_steps); % Initial intracellular calcium concentration in mM
    
    
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
    if R_star <= R_sat                % for intensities below r_sat
        c = (5.05 * log(R_star / R_half) + 230);
        tau1 = -5.263 * log(R_star / R_half) + 15.235;
        tau2 = 5.392 * log(R_star / R_half) + 40.09;
    else
        c = (5.05 * log(R_star / R_half) + 150);
        tau1 = -5.263 * log(R_star / R_half) + 15.235;
        tau2 = 5.392 * log(R_star / R_half) + 30.09;
    end
    

    % Compute Iphoto for the current intensity
    Iphoto = calc_Iphoto(time, Iamp, tau1, tau2, c ,t_del);
     
    % Add Iphoto to total photocurrent at each time step
    Iphoto_all(intensity_index, :) = Iphoto;
    

    %% Main simulation loop

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
        Ih_cone(i) = gh_cone * (1 - (1 + 3 * nh(i)) * (1 - nh(i)).^3) * (V(i) - Eh_cone); % Ih current in pA

        % Noise generation
        noise = sigma_noise_cone * sqrt(2 / tau_noise_cone) * randn() / sqrt(dt);
        I_noise(i+1) = I_noise(i) + dt * (-I_noise(i) / tau_noise_cone + noise);

        % Update membrane voltage with photocurrent from Iphoto
        dV_dt = (-IKv_cone(i) ...
                - Ih_cone(i) ...
                - ICa_cone(i) ...
                - IKCa_cone(i) ...
                - ICl_cone(i) ...
                - IL_cone(i) ...
                - Idark_cone - Iphoto(i)  ...
                - I_noise(i)) / C_cone;

        V(i+1) = V(i) + dt * dV_dt;

        % Update channel gating variables
        mKv(i+1) = mKv(i) + dt * (amKv * (1 - mKv(i)) - bmKv * mKv(i));
        hKv(i+1) = hKv(i) + dt * (ahKv * (1 - hKv(i)) - bhKv * hKv(i));
        nh(i+1) = nh(i) + dt * (anh * (1 - nh(i)) - bnh * nh(i));
        nCa(i+1) = nCa(i) + dt * (anCa * (1 - nCa(i)) - bnCa * nCa(i));

        % Update calcium concentration
        dCai_dt = -ICa_cone(i) / (volume_cone * z * F * fI * fF) - (Cai(i) - Cainf_cone) / tau_Cai_cone;
        Cai(i+1) = Cai(i) + dt * dCai_dt;
    end
    
    % Store the results for plotting
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
%% Preallocate an array to store the lowest peak value for each intensity
lowest_peaks = zeros(num_intensities, 1);

% Loop through each intensity and find the minimum value in V_all
for i = 1:num_intensities
    % Find the minimum (lowest peak) value for this intensity curve
    lowest_peaks(i) = min(V_all(i, :));
end

% Create strings for the legend that show both intensity and lowest peak values
legend_entries = cell(num_intensities, 1);
for i = 1:num_intensities
    legend_entries{i} = sprintf('R*=%d, Peak=%.2f mV', intensities(i), lowest_peaks(i));
end

% Plot the photovoltage
figure;
hold on;

% Initialize an array to hold plot handles for the legend
plot_handles = [];
for intensity_index = 1:num_intensities
    plot(time, V_all(intensity_index, :)); % Plot each curve without a specific legend entry
end

% Add a legend that includes both intensity and lowest peak values
lgd = legend(plot_handles, legend_entries, 'Location', 'best');
set(lgd, 'Units', 'normalized'); % Allows the legend to be repositioned

% Adjust plot labels and limits
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Cone Photovoltage');
xlim([t_del*0.75 t_del*1.75]);

hold off;


% Plot the Ih current for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Ih_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('Ih Current (pA)');
legend('show', 'Location', 'best');
title('Cone Ih Current ');
xlim([t_del*0.75 t_del*1.75]);
hold off;

% Plot the photocurrent (Iphoto) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, Iphoto_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('Photocurrent (pA)');
legend('show', 'Location', 'best');
title('Cone Photocurrent ');
xlim([t_del*0.75 t_del*1.75]);
hold off;



% Plot the Kv current (IKv) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKv_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('Kv Current (pA)');
legend('show', 'Location', 'best');
title('IKv Current for Different Intensities');
xlim([t_del*0.75 t_del*1.75]);
hold off;


% Plot the calcium current (ICa) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICa_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('Calcium Current (pA)');
legend('show', 'Location', 'best');
title('ICa Current ');
xlim([t_del*0.75 t_del*1.75]);
hold off;

% Plot the IKCa current for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IKCa_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('IKCa Current (pA)');
legend( 'Location', 'best');
title('IKCa Current');
xlim([t_del*0.75 t_del*1.75]);
hold off;


% Plot the Cl current (ICl) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, ICl_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('Cl Current (pA)');
legend('show', 'Location', 'best');
title('ICl Current');
xlim([t_del*0.75 t_del*1.75]);
hold off;

% Plot the leak current (IL) for different intensities
figure;
hold on;
for intensity_index = 1:num_intensities
    plot(time, IL_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('Leak Current (pA)');
legend('show', 'Location', 'best');
title('IL Current');
xlim([t_del*0.75 t_del*1.75]);
hold off;

% Plot gating variables
figure;
subplot(2, 2, 1);
hold on;
plot_handles = [];
for intensity_index = 1:num_intensities
    plot(time, mKv_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('mKv');
legend(plot_handles, 'Location', 'best');
title('Kv Activation Variable (mKv)');
xlim([t_del*0.75 t_del*1.75]);
hold off;

subplot(2, 2, 2);
hold on;
plot_handles = [];
for intensity_index = 1:num_intensities
   plot(time, hKv_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('hKv');
legend(plot_handles, 'Location', 'best');
title('Kv Inactivation Variable (hKv)');
xlim([t_del*0.75 t_del*1.75]);
hold off;

subplot(2, 2, 3);
hold on;
plot_handles = [];
for intensity_index = 1:num_intensities
    plot(time, nh_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('nh');
legend(plot_handles, 'Location', 'best');
title('Ih Activation Variable (nh)');
xlim([t_del*0.75 t_del*1.75]);
hold off;

subplot(2, 2, 4);
hold on;
plot_handles = [];
for intensity_index = 1:num_intensities
    plot(time, nCa_all(intensity_index, :), 'DisplayName', sprintf('R*=%d', intensities(intensity_index)));
end
xlabel('Time (ms)');
ylabel('nCa');
legend(plot_handles, 'Location', 'best');
title('Ca Activation Variable (nCa)');
xlim([t_del*0.75 t_del*1.75]);
hold off;


%% Plotting Intensity Response
function plot_intensity_response(Iamp_max, R_half, n)

    intensities = 0:0.1:10000;
    Iamp_all = Iamp_max * (intensities.^n ./ (intensities.^n + R_half^n));
    
    % Experimental data from Ingram et al. (2019)
    dataX = [10, 80, 90, 300, 1000, 3500, 9300];
    dataY = [0.08, 0.097, 0.1, 0.28, 0.6, 0.8, 0.9];

    % Plot intensity response
    figure;
    semilogx(intensities, Iamp_all / Iamp_max, 'DisplayName', 'Simulation');
    hold on;
    plot(dataX, dataY, 'o', 'DisplayName', 'Experiment [Ingram et al., 2019]');
    xlabel('Flash intensity (in R*/rod/flash)');
    ylabel('Normalized response (Iamp/Iamp_{max})');
    legend();
    title('Intensity Response');
    hold off;
end

% Run the simulations
intensities = [1, 2, 5, 10, 30, 100, 200, 500, 1000, 2000];
time = 0:dt:T;

plot_intensity_response(Iamp_max, R_half, n);
