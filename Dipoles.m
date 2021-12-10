% Dipoles.m
%
% Aaron James Lemmer
% March 22, 2014
%
% Find the dipoles (adjoining residues of opposite sign) in the residues
% array, remove them, and place branch cut.  The positive and negative
% residues are defined in the residues array as either +1 or -1,
% respectively, while the branch cuts are defined as binary 1s in the
% branch_cuts array.

function [residues, branch_cuts, num_dipoles] = Dipoles(num_row, num_col, branch_cuts, residues)

disp('Dipole pre-processing...');

num_res_start = nnz(residues);  % number of residues before removing dipoles equals number of nonzero elements

pos_res_ind = find(residues == 1);  % linear indices locating positive residues
neg_res_ind = find(residues == -1);  % linear indices locating negative residues

% look left
left_ind = intersect(pos_res_ind - num_row, neg_res_ind);  % linear indices of negative residues to the left of positive residues (eliminating pixels outside left edge)
branch_cuts = PlaceBranchCut2(branch_cuts, left_ind, left_ind + num_row);
residues([left_ind, left_ind + num_row]) = 0;  % remove flags for the pairs from the residues array
pos_res_ind = setdiff(pos_res_ind, left_ind + num_row);  % remove positive residues from list
neg_res_ind = setdiff(neg_res_ind, left_ind);  % remove negative residues from list

% look right
right_ind = intersect(pos_res_ind + num_row, neg_res_ind);  % linear indices of negative residues to the right of remaining positive residues (eliminating pixels outside right edge)
branch_cuts = PlaceBranchCut2(branch_cuts, right_ind, right_ind - num_row);
residues([right_ind, right_ind - num_row]) = 0;  % remove flags for the pairs from the residues array
pos_res_ind = setdiff(pos_res_ind, right_ind - num_row);  % remove positive residues from list
neg_res_ind = setdiff(neg_res_ind, right_ind);  % remove negative residues from list

% look above
above_ind = intersect(pos_res_ind - 1, neg_res_ind);  % linear indices of negative residues above remaining positive residues
[I,~] = ind2sub([num_row,num_col],above_ind);  %subscript indices of the pixels above positive residues
above_ind(I == num_row) = [];  %remove indices for pixels outside range of [1:num_pix], or above top edge?-these will appear as bottom-row indices because the index wraps around
branch_cuts = PlaceBranchCut2(branch_cuts, above_ind, above_ind + 1);
residues([above_ind, above_ind + 1]) = 0;  % remove flags for the pairs from the residues array
pos_res_ind = setdiff(pos_res_ind, above_ind + 1);  % remove positive residues from list
neg_res_ind = setdiff(neg_res_ind, above_ind);  % remove negative residues from list

% look below
below_ind = intersect(pos_res_ind + 1, neg_res_ind);  % linear indices of negative residues below remaining positive residues
[I,~] = ind2sub([num_row,num_col],below_ind);  %subscript indices of the pixels below unwrapped pixels
below_ind(I == 1) = [];  %remove indices for pixels outside range of [1:num_pix], or below bottom edge?-these will appear as top-row indices because the index wraps around
branch_cuts = PlaceBranchCut2(branch_cuts, below_ind, above_ind - 1);
residues([below_ind, below_ind - 1]) = 0;  % remove flags for the pairs from the residues array

num_dipoles = (num_res_start - nnz(residues))/2;  % number of pairs removed
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PlaceBranchCut2
%
% Aaron James Lemmer
% April 24, 2014
%
% Place a branch cut from each pixel in ind1 to the corresponding pixel in
% ind2 by assigning a binary value to the matrix branch_cuts. 

function branch_cuts = PlaceBranchCut2(branch_cuts, ind1, ind2)

%Recall that in MATLAB, when an ordered pair is used to reference a
%location in an array, it lists (row, col) which is like (y,x).  (0,0) is
%undefined and (1,1) refers to the upper left-hand corner.

%The residue is located at the intersection of a four-pixel square, and its
%is labelled as the upper-left corner of this defining 4-square.  Therefore
%the branch cut need not include the uppermost (or leftmost) pixel because 
%it needs only join the center of the four pixels that define each residue.
%See Figure 4.4, p. 106 of Ghiglia and Pritt for further information.
%
%Therefore one does not alter a pixel location if it is found to be on the
%border (if r1 | c1 | r2 | c2 == 1); in this case, the pixel is connected
%to the edge instead of another residue.

[num_row, num_col] = size(branch_cuts);

for step = 1:length(ind1)
    
    % disp(['Step = ',num2str(step),'.  Num_cuts = ',num2str(length(ind1)),'.']);  % DEBUG
    [r1, c1] = ind2sub([num_row, num_col],ind1(step));
    [r2, c2] = ind2sub([num_row, num_col],ind2(step));
    % disp(['Connecting ind1 = ',num2str(ind1(step)),' to ind2 = ',num2str(ind2(step)),'...']);  % DEBUG

    if (r2 > r1 && r1 > 1); r1 = r1 + 1;  %(r1,c1) is above (r2,c2)
    elseif (r2 < r1 && r2 > 1); r2 = r2 + 1; end;  %(r1,c1) is below (r2,c2)
    if (c2 > c1 && c1 > 1); c1 = c1 + 1;  %(r1,c1) is to the left of (r2,c2)
    elseif (c2 < c1 && c2 > 1); c2 = c2 + 1; end;  %(r1,c1) is to the right of (r2,c2)

    %disp(['(r1,c1) = (',int2str(r1),', ',int2str(c1),'); (r2,c2) = (',int2str(r2),', ',int2str(c2),');']);

    %If (r1,c1) = (r2,c2), make the single pixel a cut.  This scenario can
    %occur after modifications are made to the endpoint pixels in the previous 
    %step (see above if/elseif statements).
    if (r1 == r2 && c1 == c2); branch_cuts(r1, c1) = 1; continue; end;

    %calculate the rise and run of the line connecting the centers of (r1,c1) and (r2,c2)
    rise = abs(r2 - r1);  %if (r1 < r2); rise = r2 - r1; else rise = r1 - r2; end;
    run = abs(c2 - c1);  %if (c1 < c2); run = c2 - c1; else run = c1 - c2; end;

    %disp(['rise = ', int2str(rise), '; run = ', int2str(run),';']);

    if (run > rise)  %if true, loop through the columns (along 'x' direction)

        %decide whether to go left or right
        if (c1 < c2)  %if column of start pixel is left of column of end pixel
            c_step = 1;  %increment left to right
        else  %if column of start pixel is even with or right of column of end pixel
            c_step = -1;  %increment right to left
        end  %if (c1 < c2)

        %disp(['c_step = ', int2str(c_step), ';']);

        %calculate the slope where the columns (x) are the independent variable
        slope = (r2 - r1)/(c2 - c1);  %(delta row)/(delta col)

        %disp(['slope = ', int2str(slope), ';']);

        %locate and mark the coordinates of the branch cut
        for col = c1:c_step:c2  %for each pixel along the 'x' direction between c1 and c2

            row = round(r1 + (col - c1)*slope);  %calculate the corresponding pixel along the 'y' direction
            branch_cuts(row, col) = 1;  %mark the pixel as part of the branch cut
        end  %for col = c1:c_step:c2

    else  %if rise >= run; if true, loop through the rows (along 'y' direction)

        %decide whether to go up or down
        if (r1 < r2)  %if row of start pixel is above row of end pixel
            r_step = 1;  %increment downward
        else  %if row of start pixel is below or even with row of end pixel
            r_step = -1;  %increment upward
        end  %if (r1 < r2)

        %disp(['r_step = ', int2str(r_step), ';']);

        %calculate the slope where the rows (y) are the independent variable
        slope = (c2 - c1)/(r2 - r1);  %(delta col)/(delta row)

        %disp(['slope = ', int2str(slope), ';']);

        %locate and mark the coordinates of the branch cut
        for row = r1:r_step:r2  %for each pixel along the 'y' direction between r1 and r2

            col = round(c1 + (row - r1)*slope);  %calculate the corresponding pixel along the 'x' direction
            branch_cuts(row, col) = 1;  %mark the pixel as part of the branch cut
        end  %for row = r1:r_step:r2
        
    end  %if (run > rise)
    
end  %for step...

return
