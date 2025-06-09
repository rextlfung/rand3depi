%% Scan parameters for echo-planar imaging (EPI) sequence

% Define scanner
run('MRIsystem.m')

% Spatial parameters
res = [1.8 1.8 1.8]*1e-3; % resolution (m)
N = [120 120 80]; % acquisition tensor size
fov = N .* res; % field of view (m)
Nx = N(1); Ny = N(2); Nz = N(3);

% Random undersampling parameters. Total acceleration = Ry*Rz*caipi_z
Ry = 2; Rz = 3; % Acceleration/undersampling factors in each direction
caipi_z = 2; % Number of kz locations to acquire per shot. Must be positive integer
R = [Ry Rz];
acs = [0.1 0.1]; % Central portion of ky-kz space to fully sample
max_ky_step = round(Ny/16); % Maximum gap in fast PE direction

% Temporal parameters
Nshots = ceil(length(1:caipi_z:(Nz - caipi_z + 1))/Rz); % Number of shots per volume
minNframesPerLoop = lcm(40,Nshots)/Nshots; % number of temporal frames to complete one RF spoil cycle
NframesPerLoop = minNframesPerLoop; % 19.2 seconds ~= 1 task cycle

% Decay parameters
TE = 30e-3;                         % echo time (s)
volumeTR = 0.87;                     % temporal frame rate (s)
TR = volumeTR / Nshots;             % repetition time (s)
T1 = 1500e-3;                       % T1 (s)

% Dummy parameters
Ndummyframes = round(3*T1/TR); % dummy frames to reach steady state for calibration

% Exciting stuff
alpha = 180/pi * acos(exp(-TR/T1)); % Ernst angle (degrees)
rfDur = 2e-3;                       % RF pulse duration (s)
rfTB  = 6;                          % RF pulse time-bandwidth product
rf_phase_0 = 117;                   % RF spoiling initial phase (degrees)
NcyclesSpoil = 2;                   % number of Gx and Gz spoiler cycles

% Fat saturation
fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3;       % time-bandwidth product
fatsat.dur     = 4.0;     % pulse duration (ms)