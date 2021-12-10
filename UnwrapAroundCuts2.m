% UnwrapAroundCuts2.m
% 
% Aaron James Lemmer
% April 9, 2014
%
% Unwraps the phase data by Itoh's method without crossing any branch cuts.
% Returns the number of disconnected pieces.
%

function [num_regions, phase_soln, flag_unwrapped] = UnwrapAroundCuts2(wrapped_phase, branch_cuts, border)

%allocate memory
[num_row, num_col] = size(branch_cuts);
phase_soln = nan(num_row, num_col);
flag_unwrapped = zeros(num_row, num_col);

%determine the regions which are isolated by the branch cuts
CC = bwconncomp(1-branch_cuts,4);  %CC is a structure, use 4-connected test
num_regions = CC.NumObjects;  %number of isolated regions

% num_adjoin = 0;  %initially, there are no pixels in the adjoin list

for reg_index = 1:num_regions  %for each of the regions that are isolated by the branch cut pixels
    
    %find a starting pixel
%     start_pix = round(numel(CC.PixelIdxList{reg_index})/2);  %take a pixel somewhere in the middle of that region as the starting pixel (allows expansion of adjoin list in four directions)
    start_pix = 1;
    [r_start, c_start] = ind2sub([num_row, num_col], CC.PixelIdxList{reg_index}(start_pix));  %set the starting pixel
    
    %if the selected starting pixel is a branch cut, border, or previously unwrapped pixel
    while (branch_cuts(r_start, c_start) == 1 || flag_unwrapped(r_start, c_start) == 1 || border(r_start, c_start) == 1)
        [r_start, c_start] = ind2sub([num_row, num_col], CC.PixelIdxList{reg_index}(start_pix+1));  %take the next pixel in the list for that region as the starting pixel
    end  %if the selected starting pixel is a branch cut, border, or previously unwrapped pixel take the next pixel in the list for that region as the starting pixel
    
    phase_soln(r_start, c_start) = wrapped_phase(r_start, c_start);  %store the wrapped phase at (r_start, c_start) in the solution matrix
    flag_unwrapped(r_start, c_start) = 1;  %mark the initial pixel in the region unwrapped

    %update the adjoin list from the starting pixel (first time)
    [adjoin] = UpdateAdjoinList2(num_row, num_col, flag_unwrapped, branch_cuts);

    while (adjoin.count > 0)  %while the adjoin list is not empty (i.e., while there are more pixels left to unwrap)
        
        % unwrap left
        phase_soln(adjoin.left_ind) = ItohUnwrap(phase_soln(adjoin.left_ind + num_row), wrapped_phase(adjoin.left_ind));
        
        % unwrap right
        phase_soln(adjoin.right_ind) = ItohUnwrap(phase_soln(adjoin.right_ind - num_row), wrapped_phase(adjoin.right_ind));
        
        % unwrap up
        phase_soln(adjoin.above_ind) = ItohUnwrap(phase_soln(adjoin.above_ind + 1), wrapped_phase(adjoin.above_ind));
        
        % unwrap down
        phase_soln(adjoin.below_ind) = ItohUnwrap(phase_soln(adjoin.below_ind - 1), wrapped_phase(adjoin.below_ind));

        % mark these pixels unwrapped
        flag_unwrapped = flag_unwrapped + adjoin.list;  %mark the freshly unwrapped pixels unwrapped

        %phase_r_c = phase_soln(r,c);  %this is the phase value from the solution matrix to pass to UpdateAdjoinList (which uses it to unwrap neighboring pixels)
        [adjoin] = UpdateAdjoinList2(num_row, num_col, flag_unwrapped, branch_cuts);

    end  %while the adjoin list is not empty
end  %for each of the regions that are isolated by the branch cut pixels

%unwrap the branch cut pixels----------------------------------------------

bc_ind = find(branch_cuts);  %linear indices of the branch cut pixels
allow_ind = find(~branch_cuts & flag_unwrapped);  %linear indices of non-branch cut unwrapped pixels (cannot unwrap a branch cut pixel with another branch cut or a wrapped pixel)

%Only the pixels to the left or above matter.  Recalling the 4-square
%definition of a residue, the branch cuts lie between pixels, either to
%the right or below a given branch cut pixel.

% unwrap from left
bc_left_ind = intersect(bc_ind - num_row, allow_ind);  %linear indices of the allowable pixels to the left of branch cut pixels

phase_soln(bc_left_ind + num_row) = ItohUnwrap(phase_soln(bc_left_ind), wrapped_phase(bc_left_ind + num_row));

flag_unwrapped(bc_left_ind + num_row) = 1;  %mark these pixels unwrapped

% re-evaluate remaining unwrapped branch cut pixels
bc_ind = find(branch_cuts & ~flag_unwrapped);  %linear indices of remaining wrapped branch cut pixels
allow_ind = find(~branch_cuts & flag_unwrapped);  %updated linear indices of non-branch cut unwrapped pixels (cannot unwrap a branch cut pixel with another branch cut or a wrapped pixel)

% unwrap from above
bc_above_ind = intersect(bc_ind - 1, allow_ind);  %linear indices of the allowable pixels above branch cut pixels
[I,~] = ind2sub([num_row,num_col],bc_above_ind);  %subscript indices of the pixels above unwrapped pixels
bc_above_ind(I == num_row) = [];  %remove indices for pixels outside range of [1:num_pix], or above top edge?-these will appear as bottom-row indices because the linear indexing wraps around

phase_soln(bc_above_ind + 1) = ItohUnwrap(phase_soln(bc_above_ind), wrapped_phase(bc_above_ind + 1));

flag_unwrapped(bc_above_ind + 1) = 1;  %mark these pixels unwrapped


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UpdateAdjoinList2
%
% Aaron James Lemmer
% April 14, 2014
%
% Modifications:
% April 24, 2014
%   1.  Eliminated all calls to the function intersect, for improved speed.
%   2.  Eliminated the line to find the valid pixel indices, valid_ind (no 
%           longer needed).

function [adjoin] = UpdateAdjoinList2(num_row, num_col, flag_unwrapped, branch_cuts)

% returns a structure with 6 fields:
    % adjoin.list: an array of 4-connected pixels neighboring unwrapped
        % pixels which are allowable for unwrapping (not branch cut pixels,
        % already unwrapped, or border pixels)
        % MAY BE AN UNNECESSARY VARIABLE!!!
    % adjoin.left_ind: vector of linear indices for pixels to the left of
        % unwrapped pixels
    % adjoin.right_ind: vector of linear indices for pixels to the right of
        % unwrapped pixels
    % adjoin.above_ind: vector of linear indices for pixels above unwrapped
        % pixels
    % adjoin.below_ind: vector of linear indices for pixels below unwrapped
        % pixels
    % adjoin.count: total number of pixels in the adjoin list, adjoin.list

% allocate memory
num_pix = num_row*num_col;
adjoin.list = zeros(num_row, num_col);
adjoin.list_left = zeros(num_row, num_col);
adjoin.list_right = zeros(num_row, num_col);
adjoin.list_above = zeros(num_row, num_col);
adjoin.list_below = zeros(num_row, num_col);

valid_pixels = ~flag_unwrapped & ~branch_cuts;  %adjoin list cannot contain branch cut, border, or unwrapped pixels (border is already contained in branch_cuts)

uw_ind = find(flag_unwrapped);  %linear indices of already-unwrapped pixels

% handle the pixels to the left--------------------------------------------
adjoin.left_ind = uw_ind - num_row;  %linear indices of ALL the pixels to the left of unwrapped pixels
adjoin.left_ind(adjoin.left_ind < 1) = [];  %remove indices for pixels outside range of [1:num_pix], or off the left edge
adjoin.list_left(adjoin.left_ind) = 1;  %flag ALL of the pixels to the left of unwrapped pixels
adjoin.list_left = adjoin.list_left & valid_pixels;  %VALID pixels to the left of unwrapped pixels (logical intersection of valid pixels and all the pixels to the left of unwrapped pixels)
adjoin.left_ind = find(adjoin.list_left);  %linear indices of the VALID pixels to the left of unwrapped pixels
    
% handle the pixels to the right-------------------------------------------
adjoin.right_ind = uw_ind + num_row;  %linear indices of ALL the pixels to the right of unwrapped pixels
adjoin.right_ind(adjoin.right_ind > num_pix) = [];  %remove indices for pixels outside range of [1:num_pix], or off the right edge
adjoin.list_right(adjoin.right_ind) = 1;  %flag ALL of the pixels to the right of unwrapped pixels
adjoin.list_right = adjoin.list_right & valid_pixels;  %VALID pixels to the right of unwrapped pixels (logical intersection of valid pixels and all the pixels to the right of unwrapped pixels)
adjoin.right_ind = find(adjoin.list_right);  %linear indices of the VALID pixels to the right of unwrapped pixels

% handle the pixels above--------------------------------------------------
adjoin.above_ind = uw_ind - 1;  %linear indices of ALL the pixels above unwrapped pixels
adjoin.above_ind(adjoin.above_ind < 1) = [];  %remove indices for pixels outside range of [1:num_pix], or off the left edge
adjoin.list_above(adjoin.above_ind) = 1;  %flag ALL of the pixels above unwrapped pixels
adjoin.list_above = adjoin.list_above & valid_pixels;  %VALID pixels above unwrapped pixels (logical intersection of valid pixels and all the pixels above unwrapped pixels)
adjoin.above_ind = find(adjoin.list_above);  %linear indices of the VALID pixels above unwrapped pixels
[I,~] = ind2sub([num_row,num_col],adjoin.above_ind);  %subscript indices of the pixels above unwrapped pixels
adjoin.above_ind(I == num_row) = [];  %remove indices for pixels outside range of [1:num_pix], or above top edge?-these will appear as bottom-row indices because the index wraps around

% handle the pixels below--------------------------------------------------
adjoin.below_ind = uw_ind + 1;  %linear indices of ALL the pixels below unwrapped pixels
adjoin.below_ind(adjoin.below_ind > num_pix) = [];  %remove indices for pixels outside range of [1:num_pix], or off the left edge
adjoin.list_below(adjoin.below_ind) = 1;  %flag ALL of the pixels below unwrapped pixels
adjoin.list_below = adjoin.list_below & valid_pixels;  %VALID pixels below unwrapped pixels (logical intersection of valid pixels and all the pixels below unwrapped pixels)
adjoin.below_ind = find(adjoin.list_below);  %linear indices of the VALID pixels below unwrapped pixels
[I,~] = ind2sub([num_row,num_col],adjoin.below_ind);  %subscript indices of the pixels below unwrapped pixels
adjoin.below_ind(I == 1) = [];  %remove indices for pixels outside range of [1:num_pix], or below bottom edge?-these will appear as top-row indices because the index wraps around

% NOTE: adjoin.left_ind, adjoin.right_ind, adjoin.above_ind, and
% adjoin.below_ind may contain overlapping pixels.  This is okay, because
% we don't have a preference as to which direction we unwrap (left, right,
% up, or down) and the overlapping pixel location in the adjoin list will
% be used to specify a value from the WRAPPED phase in the Itoh unwrapper
% (as opposed to a value from the phase solution), so no significant error
% should be introduced.  Say we unwrap in the order left, right, up, down:
% for the adjoin pixels that are unwrapped multiple times due to the
% overlap, we will have final unwrapped values referenced from the pixels
% directly above.  Previous calculations for those cells will be
% overwritten.  I assume that the time spent overwriting the n overlapping
% pixels at most 3 times per loop will be insignificant compared to the
% time to locate these pixels and resolve the discrepancy (less lines of
% code this way).

adjoin.list = adjoin.list_left | adjoin.list_right | adjoin.list_above | adjoin.list_below;

adjoin.count = sum(sum(adjoin.list));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ItohUnwrap: Itoh's Method for One-Dimensional Phase Unwrapping
%
% Aaron James Lemmer
% April 14, 2014
%
%
% Given the phase value of an unwrapped starting pixel and the (wrapped) 
% phase value of a neighboring pixel, ItohUnwrapNeighbor applies Itoh's 1D 
% method to unwrap the phase.

function unwrapped_out = ItohUnwrap(unwrapped_pixels, wrapped_pixels)

%unwrap by adding the wrapped phase difference
unwrapped_out = unwrapped_pixels + Wrap(wrapped_pixels - unwrapped_pixels);