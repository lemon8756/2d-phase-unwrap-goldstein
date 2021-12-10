% PhaseUnwrap2D.m

% Created on May 28, 2015 11:05
% Created by Aaron James Lemmer


close all;    
% READ IN AND INITIALIZE RAW DATA
load phaseAng.mat;

% sub_phaseAng = phaseAng(1900:2000,1900:2000);
% phaseAng = sub_phaseAng;

% intrfrgram = intrfrgram_norm;
[num_row, num_col] = size(phaseAng);
x = (-(num_col-1):2:(num_col-1))/num_col;
y = (-(num_row-1):2:(num_row-1))/num_row;  %-1 to 1
[X,Y] = meshgrid(x,y);

mask = ones(num_row, num_col);  %sets the border to null
border = ~mask;

figure();
imagesc(phaseAng);
colorbar;

%% GOLDSTEIN'S PHASE UNWRAPPING ALGORITHM
% STEP 1: LOCATE PHASE RESIDUES
[residues, num_residues] = LocateResidues(phaseAng, border);

figure();
imagesc(residues);

%%
% STEP 2: REMOVE DIPOLE PAIRS (PRE-PROCESSING)
branch_cuts = zeros(num_row, num_col);
% [residues, branch_cuts, num_dipoles] = Dipoles(num_row, num_col, branch_cuts, residues);
% num_dipoles = 0;

% STEP 3: PLACE BRANCH CUTS
% [branch_cuts] = BranchCuts(branch_cuts, residues, num_residues - 2*num_dipoles, border);  % pass it the number of remaining (unbalanced) residues in the map
[branch_cuts] = BranchCuts(branch_cuts, residues, num_residues, border);  % pass it the number of remaining (unbalanced) residues in the map

figure();
imagesc(branch_cuts);

%%
% STEP 3: UNWRAP AROUND BRANCH CUTS
phase_soln = nan(size(branch_cuts));

[num_regions, phase_soln, flag_unwrapped] = UnwrapAroundCuts(phaseAng, phase_soln, branch_cuts, border);

figure();
imagesc(flag_unwrapped);  %Debug

figure();
imagesc(phase_soln);