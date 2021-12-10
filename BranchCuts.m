% BranchCuts.m
%
% Aaron James Lemmer
% November 12, 2013
%
% Modifications:
%     November 29, 2013: Logic was added to check if a residue found in the
%     search of the box perimeter has already been connected to the image 
%     edge.  This prevents additional cuts to the edge being placed that 
%     may isolate regions of the image unnecessarily.
%
%     December 12, 2013: Added capability to account for the image
%     border/mask.
%
%     April 25, 2014: branch_cuts is now initialized in a way to accomodate
%     dipole pre-processing.

function [branch_cuts] = BranchCuts(branch_cuts, residues, num_residues, border)

branch_cuts = branch_cuts + border;  %matrix storing the branch cuts, initially contains only dipole branch cuts from pre-processing step

[row_coord_res, col_coord_res] = find(residues ~= 0);  %find the residue coordinates
%This step essentially covers the role of the outer two for loops (in the x
%and y coordinates of the image) in Goldstein's C implementation.  The for
%loop method covers every pixel of the image, alternating between locating
%and then balancing each residue.  Here, all residues are located first,
%then all are balanced sequentially.

if (isempty(row_coord_res))  %if there are no residues, display message
    disp(['BranchCuts.m: No residues!; length(row_coord_res)=', int2str(length(row_coord_res)),...
        '; num_residues =', int2str(num_residues)])  %change this to fprintf
    return;  %invoke an early return; no residues, so no branch cuts required
end  %if (isempty(row_coord_res))

if (num_residues ~= length(row_coord_res))
    disp(['BranchCuts.m: Residue counts contradict!; num_residues =', int2str(num_residues),...
        '; length(row_coord_res) =', int2str(length(row_coord_res))])
    error('BranchCuts:ResCount','Residue counts contradict.');  %invoke an early return; error in counting residues
end  %if (num_residues ~= length(row_coord_res))

%allocate memory
[num_row, num_col] = size(residues);  %same as the size of wrapped_phase
residues_balanced = zeros(num_row, num_col);  %initially, all residues unbalanced (false = 0)
residues_active = zeros(num_row, num_col);  %initially, all residues inactive (false = 0)
connected2edge = zeros(num_row, num_col);  %initially, no residues connected to border (false = 0)

max_searchbox_rad = floor(length(residues)/2);

for current_residue = 1:num_residues  %loop through all the (initially) unbalanced residues
    
%     disp(['current_residue = ', int2str(current_residue)])
    
    row_current = row_coord_res(current_residue);  %row coordinate of the current residue
    col_current = col_coord_res(current_residue);  %column coordinate of the current residue
    
    %disp(['(row_current, col_current) = (', int2str(row_current), ',', int2str(col_current), ');'])
    
    %evaluate the box center
    if (residues_balanced(row_current, col_current) > 0)  %if current residue is already balanced (possibly by a previous residue)
        continue;  %pass to next iteration of for loop; i.e., go to next unbalanced residue
    end  %if active residue is already balanced
    if (row_current <= 1 || row_current >= num_row || col_current <= 1 || col_current >= num_col || border(row_current, col_current) == 1)  %if active residue is on image edge/border
        %If the current residue is on the image edge/border, place a branch
        %cut between the current pixel and the edge/border.
        branch_cuts(row_current, col_current) = 1;  %make this point a branchcut to the edge
        connected2edge(row_current, col_current) = 1;  %flag that the residue is connected to the edge (helps eliminate extra cuts to edge)
        residues_balanced(row_current, col_current) = 1;  %mark the current residue as balanced
        continue;  %pass to next iteration of for loop; i.e., go to next unbalanced residue
    end  %if active residue is on image edge/border
    
    %designate the current residue as balanced and active
    residues_active(row_current, col_current) = 1;  %mark this residue as active
    residues_balanced(row_current, col_current) = 1;  %mark this residue as balanced
    
    %add current pixel coordinates to active list by updating the entire list
    row_coord_act = row_current;  %add the row coordinate of the current residue to the row-coordinate active list
    col_coord_act = col_current;  %add the col coordinate of the current residue to the col-coordinate active list
    num_active = length(row_coord_act);  %count the active residues
    %num_active = sum(sum(residues_active));  %count the active residues (should = 1 at this point)
    
    charge = residues(row_current, col_current);  %store the initial residue charge for the search box
    
    %set the search box and search the box perimeter
    for radius = 1:max_searchbox_rad  %loop through search box radii
        
%         disp(['radius = ', int2str(radius),';'])
        
        ka = 1;
        while num_active ~= 0  %for each residue in the active list
            
            row_active = row_coord_act(ka);  %row coordinate of pixel corresponding to the active residue
            col_active = col_coord_act(ka);  %column coordinate of the pixel corresponding to the active residue
            
%             disp(['The search box is centered at the active pixel, (', int2str(row_active), ',', int2str(col_active),').'])
            
            %Define the n x n (n = 2*radius + 1) search box such that it is
            %centered on the active pixel (m -> rows, n-> columns, as usual).
            m_min = max(row_active - radius, 1);  %if the active row is the first row, (row_active - radius) will return an invalid (too small) row index
            m_max = min(row_active + radius, num_row);  %if the active row is the last row, (row_active + radius) will return an invalid (too large) row index
            n_min = max(col_active - radius, 1);  %if the active column is the first, (col_active - radius) will return an invalid (too small) column index
            n_max = min(col_active + radius, num_col);  %if the active column is the last, (col_active + radius) will return an invalid (too large) column index
            
            %Note: it is only necessary to search the box perimeter, as the
            %center has already been evaluated.  Therefore, search all elements
            %of the top and bottom rows, and only end elements of intermediate
            %rows.
            for m = m_min:m_max  %loop through search box rows, m
                if (m == m_min || m == m_max)  %if in top or bottom rows of search box
                    n_step = 1;  %search all elements/columns
                else  %if in intermediate rows of search box
                    n_step = n_max - n_min;  %search only edge elements/columns
                end  %if (m == m_min || m == m_max)
                for n = n_min:n_step:n_max  %loop through search box columns, n
                    
%                     disp(['Searching box pixel (m,n) = (', int2str(m), ',', int2str(n), ').'])
                    
                    if (m < 1 || m > num_row || n < 1 || n > num_col)  %if current box pixel is outside the image range, skip this box pixel
                        continue;
                    else
                        if (m == 1 || m == num_row || n == 1 || n == num_col)  %if current search box pixel is on image edge/border
                            
%                             disp(['(m,n) = (', int2str(m), ',', int2str(n), ') is on the image edge/border.'])
                            
                            if (connected2edge(row_active, col_active) == 0)  %if the current center pixel has not been connected to the image edge/border
                                [~, r_nearedge, c_nearedge] = Dist2Edge(row_active, col_active, num_row, num_col, border);
                                branch_cuts = PlaceBranchCut(branch_cuts, row_active, col_active, r_nearedge, c_nearedge);  %place a branch cut between the active pixel (box center) and the nearest image edge/border
%                                 disp(['A branch cut was placed between (', int2str(row_active), ',', int2str(col_active), ') and (', int2str(r_nearedge), ',', int2str(c_nearedge), ').'])
                                connected2edge(row_active, col_active) = 1;  %flag that the residue is connected to the edge (helps eliminate extra cuts to edge)
%                                 disp(['(', int2str(row_active), ',', int2str(col_active), ') is connected to the image edge/border.'])
                            end  %if the current center pixel has not been connected to the image border
                            
                            charge = 0;  %neutralize the search box charge
                            
%                             disp(['The total charge = ', int2str(charge), '.'])
                            
                        elseif (residues(m,n) ~= 0 && residues_active(m,n) == 0)  %else if current search box pixel is an inactive residue (first logical test is for residue, second for inactive)
                            
%                             disp(['(m,n) = (', int2str(m), ',', int2str(n), ') is an inactive residue.'])
                            
                            if (residues_balanced(m,n) == 0)  %if inactive residue is unbalanced
%                                 disp(['(m,n) = (', int2str(m), ',', int2str(n), ') is unbalanced.'])
                                charge = charge + residues(m,n);  %add charge of residue (+1 or -1) to charge counter
                                residues_balanced(m,n) = 1;  %mark the box pixel residue as balanced
                            end %end if inactive residue is unbalanced
                            
%                             disp(['(m,n) = (', int2str(m), ',', int2str(n), ') is balanced.'])
                            
                            residues_active(m,n) = 1;  %mark this inactive residue as active
                            %MAY BE A BETTER ALTERNATIVE WHERE I
                            %PREALLOCATE THE MEMORY FOR THESE NEXT TWO
                            %LINES OF CODE, IMPROVING SPEED.
                            row_coord_act = [row_coord_act, m];  %add m to the row-coordinate active list
                            col_coord_act = [col_coord_act, n];  %add n to the col-coordinate active list
                            num_active = length(row_coord_act);  %recount the active residues
%                             disp(['There are ', int2str(num_active), ' active residues.'])
                            branch_cuts = PlaceBranchCut(branch_cuts, row_active, col_active, m, n);  %place a branch cut between the search box pixel and the box center
%                             disp(['A branch cut was placed between (', int2str(row_active), ',', int2str(col_active), ') and (', int2str(m), ',', int2str(n), ').'])
                            if (connected2edge(m,n) == 1)  %if the box pixel (m,n) is already connected to the edge, the box center is now balanced upon connecting to (m,n)
                                connected2edge(row_active, col_active) = 1;  %(row_active, col_active) is now connected to the edge
                                charge = 0;
                            end  %if (m,n) already connected to the edge
%                             disp(['The total charge = ', int2str(charge), '.'])
                            
                        end %if current search box pixel is on image edge/border, else if it is an inactive residue
                        if (charge == 0); break; end;  %once the charge is
                        %neutralized, the search of the box is complete; move to
                        %next unbalanced residue
                    end %if current box pixel is outside the image range, skip this box pixel; else...
                end  %for n
                if (charge == 0); break; end;  %once the charge is neutralized,
                %the search of the box is complete; move to next unbalanced
                %residue (this second conditional break is necessary because
                %the first only exits the inner (n) loop)
            end  %for m

            if (charge == 0); break; end;  %once the charge is neutralized, the
            %search of the box is complete; move to next unbalanced residue
            %(this third conditional break is necessary in order to exit to the
            %for radius loop)
            
            ka = ka + 1;  %go to the next active residue
            if (ka > num_active); break; end;  %if there is not another active residue, go to the next box size
            
        end  %while num_active ~= 0 (for each residue in active list)
        
        if (charge == 0); break; end;  %once the charge is neutralized, the
            %search of the box is complete; move to next unbalanced residue
            %(this fourth conditional break is necessary in order to exit the
            %for radius loop, allowing iteration to the next unbalanced residue)
    end  %for radius = 1:max_searchbox_rad
    
    if (charge ~= 0)  %if, after iterating through box sizes, still haven't neutralized charge
%         disp(['After iterating through box sizes, the charge was not neutralized; the charge = ', int2str(charge), '.'])
        for ka = 1:num_active  %for each residue in the active list
            row_active = row_coord_act(ka);  %row coordinate of pixel corresponding to the active residue
            col_active = col_coord_act(ka);  %column coordinate of the pixel corresponding to the active residue
            [~, r_nearedge, c_nearedge] = Dist2Edge(row_active, col_active, num_row, num_col, border);
        	branch_cuts = PlaceBranchCut(branch_cuts, row_active, col_active, r_nearedge, c_nearedge);  %place branch cut to nearest pixel on border
        end  %for each residue in the active list
    end  %if (charge ~= 0)
    
    %mark all active pixels inactive
    for ka = 1:num_active
        residues_active(row_coord_act(ka), col_coord_act(ka)) = 0;  %mark all residues inactive
    end  %for ka = 1:num_active
    if (sum(residues_balanced(:)) >= num_residues); break; end;  %if all residues are balanced, no need to keep iterating through them!
end  %for i = 1:num_residues


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dist2Edge
%
% Aaron James Lemmer
% November 18, 2013
%
% Passed the coordinates of a pixel (r,c), Dist2Edge returns the
% coordinates of the nearest edge/border pixel (r_nearedge, c_nearedge) and 
% the square of the distance between (r,c) and (r_nearedge, c_nearedge).
%
% Modifications:
%     December 12, 2013: Added capability to account for the image
%     border/mask.

function [min_dist2, r_nearedge, c_nearedge] = Dist2Edge(r, c, num_row, num_col, border)

max_boxrad = num_row + num_col;  %this is merely any large number that puts an upper bound on the search box radius

for boxrad = 0:max_boxrad  %loop through search box sizes
    found = 0;  %flag that keeps track of whether an edge/border pixel has been found (false := 0)
    min_dist2 = max_boxrad^2;  %initialize the minumum found squared distance to a large number that exceeds maximum value that could be possible in given image
    
    %Define the n x n (n = 2*radius + 1) search box such that it is
    %centered on the active pixel (m -> rows, n-> columns, as usual).
    m_min = max(r - boxrad, 1);  %if r is the first row, (r - boxrad) will return an invalid (too small) row index
    m_max = min(r + boxrad, num_row);  %if r is the last row, (r + boxrad) will return an invalid (too large) row index
    n_min = max(c - boxrad, 1);  %if c is the first column, (c - boxrad) will return an invalid (too small) column index
    n_max = min(c + boxrad, num_col);  %if c is the last column, (c + boxrad) will return an invalid (too large) column index
    
    %Note: it is only necessary to search the box perimeter, as the
    %center has already been evaluated.  Therefore, search all elements
    %of the top and bottom rows, and only end elements of intermediate
    %rows.
    for m = m_min:m_max  %loop through search box rows, m
        if (m == m_min || m == m_max)  %if in top or bottom rows of search box
            n_step = 1;  %search all elements/columns
        else  %if in intermediate rows of search box
            n_step = n_max - n_min;  %search only edge elements/columns
        end  %if (m == m_min || m == m_max)
        for n = n_min:n_step:n_max  %loop through search box columns, n
            
            if (m <= 1 || m >= num_row || n <= 1 || n >= num_col || border(m,n) == 1)  %if (m,n) is on the edge/border of the image
                found = 1;  %change flag to true
                dist2 = (m - r)*(m - r) + (n - c)*(n - c);  %calculate the squared of the distance between (m,n) and (r,c)
                if (dist2 < min_dist2)
                    min_dist2 = dist2;
                    r_nearedge = m;
                    c_nearedge = n;
                end  %if (dist2 < min_dist2)
            end  %if (m <= 1 || m >= num_row || n <= 1 || n >= num_col)
        end  %for n
    end  %for m
    
    if (found == 1); break; end;  %if an edge/pixel is found at the current box size, no need to keep searching 
    
end  %for boxrad = 0:max_boxrad
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PlaceBranchCut
%
% Aaron James Lemmer
% November 13, 2013
%
% Place a branch cut from the pixel (r1,c1) to the pixel (r2,c2) by
% assigning a binary value to the matrix branch_cuts. 

function branch_cuts = PlaceBranchCut(branch_cuts, r1, c1, r2, c2)

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

if (r2 > r1 && r1 > 1); r1 = r1 + 1;  %(r1,c1) is above (r2,c2)
elseif (r2 < r1 && r2 > 1); r2 = r2 + 1; end;  %(r1,c1) is below (r2,c2)
if (c2 > c1 && c1 > 1); c1 = c1 + 1;  %(r1,c1) is to the left of (r2,c2)
elseif (c2 < c1 && c2 > 1); c2 = c2 + 1; end;  %(r1,c1) is to the right of (r2,c2)

%disp(['(r1,c1) = (',int2str(r1),', ',int2str(c1),'); (r2,c2) = (',int2str(r2),', ',int2str(c2),');']);

%If (r1,c1) = (r2,c2), make the single pixel a cut.  This scenario can
%occur after modifications are made to the endpoint pixels in the previous 
%step (see above if/elseif statements).
if (r1 == r2 && c1 == c2); branch_cuts(r1, c1) = 1; return; end;

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

return