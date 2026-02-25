% =========================================================================
% Canopy-SIP Model (Optical BRF Simulation)
% 
% DESCRIPTION:
% This script demonstrates the pure optical bidirectional reflectance factor 
% (BRF) simulation for a discrete vegetation canopy. It integrates the 
% Geometric-Optical (GO) theory (for four scene components) and spectral 
% invariants theory (p-theory) for multiple scattering.
%
% AUTHORS:
% Yelu Zeng, Min Chen, Dalei Hao, Yachang He
% Updated by: Yachang He (Email: akhyc13@gmail.com)
% Date: 2024-08-01 / Refined for GitHub Open Source: 2026-02
%
% REFERENCE:
% Please cite the corresponding publications related to the SIP model 
% modifications for discrete canopies.
%
%
% =========================================================================

clear all;
clc;

%% 1. Prepare for Simulation

% Load structural data (Gap fractions and Clumping Index)
% Note: Forward slashes used for Mac/Linux/Windows cross-platform compatibility
load('./CI_2/gap_tot.mat');     % Total gap fraction (within + between)
load('./CI_2/gap_within.mat');  % Within-crown gap fraction
load('./CI_2/gap_betw.mat');    % Between-crown gap fraction
load('./CI_2/CI_within.mat');   % Within-crown clumping index

% 1.1 Sun-Sensor Geometry
% Sensor geometry: 13 angles in the principal plane
va = zeros(13, 4);
for t = 1:7
    va(t, 1) = 10 * (7 - t); % Forward view zenith angles
    va(t, 2) = 0;            % Forward view azimuth angle
end
for t = 1:6
    va(t + 7, 1) = 10 * t;   % Backward view zenith angles
    va(t + 7, 2) = 180;      % Backward view azimuth angle
end

SZA = 0; % Sun Zenith Angle [degrees]
SAA = 0; % Sun Azimuth Angle [degrees]

% 1.2 Vegetation Type and Structural Parameters
pathAngle = 1; % 0 = homogeneous canopy, 1 = discrete canopy

if pathAngle ~= 0 
    % Gap fractions for discrete canopy at Nadir (Index 7)
    gap_H        = gap_betw(7, 3);      % Nadir between-crown gap fraction
    gap_H_within = gap_within(7, 3);    % Nadir within-crown gap fraction
    gap_H_tot    = gap_H + gap_H_within;% Nadir total gap fraction
    CI_H_within  = CI_within(7, 3);     
    
    gap_S        = gap_betw(7, 3);      % Solar direction between-crown gap fraction
    gap_S_within = gap_within(7, 3);    
    gap_S_tot    = gap_S + gap_S_within;% Solar direction total gap fraction
end

Crowndeepth = 12.8675;   % Average crown depth [m]
Height      = 20;        % Canopy height [m]
Height_c    = 6.634;     % Crown center height [m]
dthr        = 0.41234;   % Diameter To Height Ratio
bl          = 0.1;       % Leaf width [m]
HotSpotPar  = bl / Height;

c1 = CI_H_within;        % Nadir within-crown clumping index
c2 = 1 - gap_H; 
c  = c1 * c2;            % Canopy-level continuity parameter

iD   = 0.58073;          % Scene hemispherical interceptance
LAI  = 5;                % Single crown sphere LAI [m2/m2]
FAVD = 0.375;            % Foliage Area Volume Density
D    = 0;                % Ratio of diffuse to incoming irradiance

% 1.3 Soil and Canopy Optical Properties
TypeLidf = 2; % 1 = Two-parameter, 2 = Single-parameter (Campbell)

if (TypeLidf == 1)    
    LIDFa = -0.35;
    LIDFb = -0.15;
    [lidf, litab] = dladgen(LIDFa, LIDFb);   
elseif (TypeLidf == 2)
    LIDFa = 57.3; % Average leaf angle (Spherical distribution)
    LIDFb = 0;
    [lidf, litab] = campbell(LIDFa);        
end
lidf = lidf';

% Single band optics (Typical NIR band values)
rho = 0.4957; % Leaf reflectance
tau = 0.4409; % Leaf transmittance 
w   = rho + tau; % Leaf single scattering albedo
rg  = 0.159;  % Soil reflectance  

go_par = dthr * Crowndeepth; % Hotspot parameter at crown scale

%% 2. Start BRF Simulation
BRF3 = ones(13, 1);
deg2rad = pi / 180;

for t = 1:13   
    
    tts = SZA;      
    tto = va(t, 1); 
    psi = va(t, 2); 
    
    % Ensure relative azimuth angle is symmetric
    if psi > 180
        psi = psi - 360; 
    end
    psi = abs(psi);    
    psi = abs(psi - 360 * round(psi / 360));  
    
    CIy1 = 1;
    CIy2 = 1;
    [CIs] = CIxy(CIy1, CIy2, tts);
    [CIo] = CIxy(CIy1, CIy2, tto);
    
    % 2.1 Directional gap extraction
    gap_V_tot    = gap_tot(t, 3);
    gap_V_within = gap_within(t, 3);
    gap_V_betw   = gap_betw(t, 3);
    
    Ps_dir_go = gap_S;       % Between-crown gap fraction in solar direction
    Pv_dir_go = gap_V_betw;  % Between-crown gap fraction in view direction
    
    % 2.2 Calculate four GO components (Kc, Kt, Kg, Kz)
    Kg = Ps_dir_go * Pv_dir_go + get_HSF_go(go_par, SZA, SAA, tto, psi, Ps_dir_go, Pv_dir_go, Height_c);
    Kz = Pv_dir_go - Kg;   
    Kct = 1 - Pv_dir_go;
    
    delta_angle = cosd(SZA) * cosd(tto) + sind(SZA) * sind(tto) * cosd(psi - SAA); 
    phi = acosd(delta_angle);
    delta_val = cosd(phi .* (1 - sin(pi .* c ./ 2)));   
    
    if ((Height - Crowndeepth) < Crowndeepth) && (tto > SZA) && (SAA == psi)   
        Kc = Kct; % Continuous canopy
    else
        Kc = 0.5 * (1 + delta_val) * Kct; % Discrete canopy
    end
    Kt = Kct - Kc;
    
    % 2.3 Calculate BRF_L: Vegetation single scattering contribution
    Ps_dir_inKz = gap_S_within;     
    
    [Gs, Go, k, K, sob, sof] = PHASE(tts, tto, psi, lidf);
    
    [kc, kg]       = sunshade_H(tts, tto, psi, Gs, Go, CIs, CIo, LAI, HotSpotPar);   
    [kc_kt, kg_kt] = sunshade_Kt_He(tts, tto, psi, Gs, Go, CIs, CIo, LAI);   
    
    wso = sob * rho + sof * tau; % Bidirectional scattering coefficient
    
    BRF_v1    = wso .* kc / K;                            % Sunlit crown
    BRF_v1_kt = sqrt(Ps_dir_inKz) .* wso .* kc_kt / K;    % Shaded crown
    
    % 2.4 Calculate BRF_S: Soil contribution
    BRFsc   = kg * rg;       % Sunlit crown visible sunlit soil
    BRFs_kt = kg_kt * rg;    % Shaded crown visible sunlit soil
    
    % 2.5 Calculate BRF_M: Multiple scattering contribution (p-theory)
    i0 = 1 - gap_S_tot;
    i0 = D * iD + (1 - D) * i0;
    iv = 1 - gap_V_tot;
    t0 = 1 - i0;
    tv = 1 - iv;
    
    p         = 1 - iD / LAI;
    rho2      = iv / 2 / LAI;
    rho_hemi2 = iD / 2 / LAI;
    
    Tdn   = t0 + i0 * w * rho_hemi2 ./ (1 - p * w);
    Tup_o = tv + iD * w * rho2 ./ (1 - p * w);
    Rdn   = iD * w * rho_hemi2 ./ (1 - p * w);
    
    BRF_vm = i0 * w.^2 * p * rho2 ./ (1 - p * w);          % Vegetation multiple scattering
    BRFm   = rg .* Tdn .* Tup_o ./ (1 - rg .* Rdn) - t0 * rg * tv; % Vegetation-Soil multiple interactions

    %% 3. Final: Calculate total BRF (Discrete Canopy)
    BRF3(t, :) = Kc .* (BRFsc + BRF_v1) ...         
               + Kt .* (BRFs_kt + BRF_v1_kt) ...
               + (Kg + Kz .* Ps_dir_inKz) .* rg ... 
               + BRFm + BRF_vm;
end

% Save the BRF result for the NIR band at SZA = 0
save('BRF_SIP_SZA00_Nir2.mat', 'BRF3'); 
disp('Simulation completed successfully. Results saved to BRF_SIP_SZA00_Nir.mat');


%% 4. Plot Results
disp('Generating BRF plot in the principal plane...');

signed_vza = zeros(13, 1);
signed_vza(1:7) = -va(1:7, 1);      % 前向: -60, -50, -40, -30, -20, -10, 0
signed_vza(8:13) = va(8:13, 1);     % 后向: 10, 20, 30, 40, 50, 60

figure('Name', 'Canopy BRF (NIR)', 'Color', 'w');
plot(signed_vza, BRF3, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');


xlabel('View Zenith Angle (\circ)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bidirectional Reflectance Factor (BRF)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Canopy BRF in Principal Plane (NIR, SZA = %d\\circ)', SZA), 'FontSize', 14);



grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1, 'XTick', -60:20:60);
xlim([-65 65]);
ylim([0.1 0.5]);

disp('Plot generated successfully.');