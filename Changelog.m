%% Change Log for ConeRodModel
%
% TODO: Allow exactly 50% of rods stimulated below each cone but with
% random identity
% TODO: Test connectivity with 1rod 1cone bug
%
% Version 3.4 - Created on 2024-09-09 by Luca Della Santina
%
% + Ported to MATLAB R2024a
%
% Version 3.3 - Created on 2024-09-02 by Luca Della Santina
%
% + Cleaned code for sharing
% + Added GPL v3 License
%
% Version 3.2 - Created on 2024-02-06 by Luca Della Santina
%
% + Simplified S_Rod graph creation by manually building diagonal matrix
% + New button to assign current prefs to selected batch item
% + RodsIllRandom now randomized position of rods illuminated within block
%
% Version 3.1 - Created on 2024-01-31 by Luca Della Santina
%
% + User can define a custon number of illuminated rods per cone
% + Illuminated rods are in brighter blue on the graph plot
% + User can choose random chance 50% of a rod being illuminated
% + User can specify Mean and SD of normal distribution for random RCjc
% + Ensure range validity from/to/step before adding to batch list
% + Ensure reference cone number exists before solving or plotting
% + Save table with all rods response and whether connected to ref cone
% + Info button with description of each parameter
%
% Version 3.0 - Created on 2023-12-27 by Luca Della Santina
%
% + Rewrote Cone_Graph and Rod_Graph to allow arbitrary # of rods/cones
% + Fixed legend color of rod spheres
% + Settings ConeRows, ConeCols, RodRows, RodCols to specify amount of PRs
% + Right-click in the figure area allows toggle between 3D/2D graph plot
%
% Version 2.2 - Created on 2023-12-26 by Luca Della Santina
%
% + Reference Cone # can be specified in the model settings
% + Reference cone is plotted in orange in the connectivity graph
%
% Version 2.1 - Created on 2023-12-26 by Luca Della Santina
%
% + Fixed error in filenames containing ":" not allowed in windows OS
% + Preallocated some of the matrices for speed during initialization
%
% Version 2.0 - Created on 2023-12-21 by Luca Della Santina
%
% + Ported to MATLAB R2023b
% + Migrated graph plot from discontinued biograph() to graph()
% + Colored cones in Red, central cone in Magenta, Rods in blue
% + All panels are saved as single PDF
% + Saved results and plots filenames are timestamped
% + exposed t1,t2 SPR (single photon response) parameters
% + all parameters in an editable table

% Version 1.0 - Created on 2023-12-19 by Luca Della Santina
%
% Initial implementation from Christophe scripts
% Requires MATLAB R2022a with Bioinformatics toolbox for biograph()

%% TODO implement or discuss
%
% Initialization code (old Constants.m):
%   Allow variable number of rods coupled for each cone with constant coupling, 
%
% Model solution code (old cone_rod_solver.m)
%   Resulting voltages are corrected by Ec and Er (membrane resting potential cone/rod), why those are both zero? Need to test
%
% Connectivity graph
%
% Saving results table
%

%% General thoughts

% myode() defines the differential equation function
% dVdt = (S_total*V-I)/C 
%
% where:
% S_total is a matrix calculated in cone_graph and rod_graph
% C is the membrane capacitance in each cone and rod (1..9=cones, 10..369=rods)
% V is the voltage at time t
% I is the current in each cone or rod at time different times (interpolated from Itc,IC for cones or Itr,IR for rods)
