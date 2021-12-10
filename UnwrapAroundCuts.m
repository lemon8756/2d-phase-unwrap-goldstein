% UnwrapAroundCuts.m
% 
% Aaron James Lemmer
% November 19, 2013
%
% Unwraps the phase data by Itoh's method without crossing any branch cuts.
% Returns the number of disconnected pieces.
%
% Modifications:
%     December 12, 2013: Added capability to account for the image
%     border/mask.

function [num_regions, phase_soln, flag_unwrapped] = UnwrapAroundCuts(wrapped_phase, phase_soln, branch_cuts, border)

%allocate memory
[num_row, num_col] = size(branch_cuts);
adjoin_list = zeros(num_row, num_col);
flag_unwrapped = zeros(num_row, num_col);

%determine the regions which are isolated by the branch cuts
CC = bwconncomp(1-branch_cuts,4);  %CC is a structure, use 4-connected test
num_regions = CC.NumObjects;  %number of isolated regions

num_adjoin = 0;  %initially, there are no pixels in the adjoin list

for reg_index = 1:num_regions  %for each of the regions that are isolated by the branch cut pixels
    
    %find a starting pixel
    [r_start, c_start] = ind2sub(size(1-branch_cuts), CC.PixelIdxList{reg_index}(1));  %just take the first pixel in the list for that region as the starting pixel

    if (branch_cuts(r_start, c_start) ~= 1 && flag_unwrapped(r_start, c_start) ~= 1 && border(r_start, c_start) ~= 1)  %if the selected starting pixel is not a branch cut, border, or previously unwrapped pixel

        phase_soln(r_start, c_start) = wrapped_phase(r_start, c_start);  %store the wrapped phase at (r_start, c_start) in the solution matrix
        phase_r_c = phase_soln(r_start, c_start);  %pass this initial solution point to UpdateAdjoinList for the first time
        flag_unwrapped(r_start, c_start) = 1;  %mark the initial pixel in the region unwrapped
        %disp(['(', int2str(r_start), ',', int2str(c_start), ') just unwrapped and passed to UpdateAdjoinList.']);  %DEBUG: current pixel being passed to UpdateAdjoinList

        %update the adjoin list from the starting pixel
        [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = UpdateAdjoinList(r_start, c_start, phase_r_c, wrapped_phase, phase_soln, adjoin_list, num_adjoin, flag_unwrapped, branch_cuts, border);
        
        %whilecount = 0;  %DEBUG: counter for the while loop
        while (num_adjoin > 0)  %while the adjoin list is not empty
            
            %whilecount = whilecount + 1;  %DEBUG: increment the counter
            %disp(['While loop iteration number ', int2str(whilecount), ':']);  %DEBUG: print the while loop counter
            
            [empty_flag, r, c, num_adjoin] = GetNext2Unwrap(adjoin_list, num_adjoin);  %fetch the next pixel from the adjoin list and decrement the number of pixels in the list
            if (empty_flag == 1); break; end;  %if the adjoin list is empty, there are no more pixels to unwrap

            %flag_unwrapped(r,c) = 1;  %mark the freshly fetched pixel (r,c) unwrapped
            %disp(['(', int2str(r), ',', int2str(c), ') just passed to UpdateAdjoinList.']);  %DEBUG: current pixel being passed to UpdateAdjoinList

            phase_r_c = phase_soln(r,c);  %this is the phase value from the solution matrix to pass to UpdateAdjoinList (which uses it to unwrap neighboring pixels)
            [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = UpdateAdjoinList(r, c, phase_r_c, wrapped_phase, phase_soln, adjoin_list, num_adjoin, flag_unwrapped, branch_cuts, border);

        end  %while the adjoin list is not empty
    end  %if the selected starting pixel is not a branch cut or previously unwrapped pixel
end  %for each of the regions that are isolated by the branch cut pixels

%unwrap the branch cut pixels
[row_bc, col_bc] = find(branch_cuts);  %find the row and column indices of the branch cuts
num_bc = length(row_bc);  %count the branch cut pixels

for b = 1:num_bc  %for each branch cut pixel, unwrap
    
    r_bc_current = row_bc(b);  %the row coordinate of the current branch cut pixel
    c_bc_current = col_bc(b);  %the row coordinate of the current branch cut pixel
    
    %Only the pixels to the left or above matter.  Recalling the 4-square
    %definition of a residue, the branch cuts lie between pixels, either to
    %the right or below a given branch cut pixel.
    
    %the coordinates of the pixel to the left of (r_bc_current,c_bc_current)
    r_left = r_bc_current;
    c_left = c_bc_current - 1;
    
    %the coordinates of the pixel above (r_bc_current,c_bc_current)
    r_up = r_bc_current - 1;
    c_up = c_bc_current;
    
    if (c_left >= 1 && branch_cuts(r_left,c_left) ~= 1)  %if the pixel to the left is valid for unwrapping
        
        phase_soln(r_bc_current,c_bc_current) = ItohUnwrapNeighbor(phase_soln(r_left,c_left), wrapped_phase(r_bc_current,c_bc_current));
        
    elseif (r_up >= 1 && branch_cuts(r_up,c_up) ~= 1)  %if the pixel above is valid for unwrapping
        
        phase_soln(r_bc_current,c_bc_current) = ItohUnwrapNeighbor(phase_soln(r_up,c_up), wrapped_phase(r_bc_current,c_bc_current));
        
    end  %if there is a valid neighboring pixel for unwrapping
    
end  %for each branch cut pixel 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GetNext2Unwrap
%
% Aaron James Lemmer
% November 19, 2013

function [empty_flag, r, c, num_adjoin] = GetNext2Unwrap(adjoin_list, num_adjoin)

empty_flag = 0;  %default is that the adjoin list is not empty

if (num_adjoin < 1)  %if the adjoin list is empty
    empty_flag = 1;
    return;  %invoke an early return
end  %if the adjoin list is empty

%find the indices of the nonzero elements of adjoin_list
[row_list, col_list] = find(adjoin_list);

%the next pixel to unwrap is the last one in the list
r = row_list(num_adjoin);  %row coordinate of the last pixel in the adjoin list
c = col_list(num_adjoin);  %col coordinate of the last pixel in the adjoin list

num_adjoin = num_adjoin - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UpdateAdjoinList
%
% Aaron James Lemmer
% November 19, 2013
%
%
% UpdateAdjoinList is passed the coordinates to a pixel, (r,c).  If the 
% following conditions are true:
%   1. (r,c) is not already on the adjoin list,
%   2. (r,c) is not yet unwrapped,
%   3. (r,c) is not a branch cut
%   4. (r,c) is not a border pixel (it MAY be an edge pixel),
% then UpdateAdjoinList does the following:
%   1. unwraps (r,c)
%   2. adds (r,c) to the adjoin list via Insert2AdjoinList.
%   Insert2AdjoinList, in addition to adding (r,c) to the list,
%       1. stores the phase value at (r,c) in the solution matrix, and
%       2. marks (r,c) as unwrapped.
%
% Below is a diagram of the pixels adjoining (r,c):
% 
%            -----(r-1,c)-----
%            |       |       |
%            |       |       |
%         (r,c-1)--(r,c)--(r,c+1)
%            |       |       |
%            |       |       |
%            -----(r+1,c)-----
%
% Modifications:
%     December 12, 2013: Added capability to account for the image
%     border/mask.  Also fixed an error where only the the first of the 
%     four neighboring pixels was queried to see if it was already on the
%     adjoin list (all are now checked).

function [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = UpdateAdjoinList(r, c, phase_r_c, wrapped_phase, phase_soln, adjoin_list, num_adjoin, flag_unwrapped, branch_cuts, border)

[num_row, num_col] = size(branch_cuts);

%handle the pixel to the left of (r,c)-------------------------------------
m = r;
n = c - 1;
%if (m,n) is not off the left edge of the image and not a branch cut pixel 
%and not a border pixel nor already unwrapped
if (n >= 1 && adjoin_list(m,n) ~= 1 && flag_unwrapped(m,n) ~= 1 && branch_cuts(m,n) ~= 1 && border(m,n) ~= 1)
    unwrapped_m_n = ItohUnwrapNeighbor(phase_r_c, wrapped_phase(m,n));  %unwrap the pixel (m,n)
    
    %insert the pixel into the adjoin list and mark unwrapped
    [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = Insert2AdjoinList(phase_soln, m, n, unwrapped_m_n, adjoin_list, num_adjoin, flag_unwrapped);
% elseif (n < 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is outside the image.']);
% elseif (flag_unwrapped(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already unwrapped.']);
% elseif (branch_cuts(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is a branch cut.']);
% elseif (adjoin_list(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already on the adjoin list.']);
end
    
%handle the pixel to the right of (r,c)------------------------------------
m = r;
n = c + 1;
%if (m,n) is not off the right edge of the image and not a branch cut pixel
%and not a border pixel nor already unwrapped
if (n <= num_col && adjoin_list(m,n) ~= 1 && flag_unwrapped(m,n) ~= 1 && branch_cuts(m,n) ~= 1  && border(m,n) ~= 1)
    unwrapped_m_n = ItohUnwrapNeighbor(phase_r_c, wrapped_phase(m,n));  %unwrap the pixel (m,n)
    
    %insert the pixel into the adjoin list and mark unwrapped
    [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = Insert2AdjoinList(phase_soln, m, n, unwrapped_m_n, adjoin_list, num_adjoin, flag_unwrapped);
% elseif (n > num_col)
%     disp(['(', int2str(m), ',', int2str(n), ') is outside the image.']);
% elseif (flag_unwrapped(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already unwrapped.']);
% elseif (branch_cuts(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is a branch cut.']);
% elseif (adjoin_list(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already on the adjoin list.']);
end

%handle the pixel above (r,c)----------------------------------------------
m = r - 1;
n = c;
%if (m,n) is not off the top edge of the image and not a branch cut pixel 
%and not a border pixel nor already unwrapped
if (m >= 1 && adjoin_list(m,n) ~= 1 && flag_unwrapped(m,n) ~= 1 && branch_cuts(m,n) ~= 1 && border(m,n) ~= 1)
    unwrapped_m_n = ItohUnwrapNeighbor(phase_r_c, wrapped_phase(m,n));  %unwrap the pixel (m,n)
    
    %insert the pixel into the adjoin list and mark unwrapped
    [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = Insert2AdjoinList(phase_soln, m, n, unwrapped_m_n, adjoin_list, num_adjoin, flag_unwrapped);
% elseif (m < 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is outside the image.']);
% elseif (flag_unwrapped(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already unwrapped.']);
% elseif (branch_cuts(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is a branch cut.']);
% elseif (adjoin_list(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already on the adjoin list.']);
end

%handle the pixel below (r,c)----------------------------------------------
m = r + 1;
n = c;
%if (m,n) is not off the bottom edge of the image and not a branch cut 
%pixel and not a border pixel nor already unwrapped
if (m <= num_row && adjoin_list(m,n) ~= 1 && flag_unwrapped(m,n) ~= 1 && branch_cuts(m,n) ~= 1 && border(m,n) ~= 1)
    unwrapped_m_n = ItohUnwrapNeighbor(phase_r_c, wrapped_phase(m,n));  %unwrap the pixel (m,n)
    
    %insert the pixel into the adjoin list and mark unwrapped
    [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = Insert2AdjoinList(phase_soln, m, n, unwrapped_m_n, adjoin_list, num_adjoin, flag_unwrapped);
% elseif (m > num_row)
%     disp(['(', int2str(m), ',', int2str(n), ') is outside the image.']);
% elseif (flag_unwrapped(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already unwrapped.']);
% elseif (branch_cuts(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is a branch cut.']);
% elseif (adjoin_list(m,n) == 1)
%     disp(['(', int2str(m), ',', int2str(n), ') is already on the adjoin list.']);
end

adjoin_list(r,c) = 0;  %remove the center pixel from the adjoin list


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ItohUnwrapNeighbor: Itoh's Method for One-Dimensional Phase Unwrapping
%
% Aaron James Lemmer
% November 20, 2013
%
%
% Given the phase value of an unwrapped starting pixel and the (wrapped) 
% phase value of a neighboring pixel, ItohUnwrapNeighbor applies Itoh's 1D 
% method to unwrap the phase.

function unwrapped = ItohUnwrapNeighbor(unwrapped_pixel, wrapped_pixel)

%unwrap by adding the wrapped phase difference
unwrapped = unwrapped_pixel + Wrap(wrapped_pixel - unwrapped_pixel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Insert2AdjoinList
%
% Aaron James Lemmer
% November 19, 2013

function [phase_soln, adjoin_list, num_adjoin, flag_unwrapped] = Insert2AdjoinList(phase_soln, r, c, unwr_val, adjoin_list, num_adjoin, flag_unwrapped)

phase_soln(r,c) = unwr_val;  %store the unwrapped value in the solution matrix
% disp(['(', int2str(r), ',', int2str(c), ') was unwrapped.']);

adjoin_list(r,c) = 1;  %add the pixel to the adjoin list
% disp(['(', int2str(r), ',', int2str(c), ') was added to the adjoin list.']);
num_adjoin = num_adjoin + 1;  %recount the number of pixels in the list
% disp(['There are ', int2str(num_adjoin), ' pixels in the adjoin list.']);

flag_unwrapped(r,c) = 1;  %mark the pixel unwrapped