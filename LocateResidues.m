% LocateResidues.m
%
% Aaron James Lemmer
% November 8, 2013
%
% Detect residues in wrapped phase data and mark them as positive or 
% negative residues.  Ignore the pixels marked as part of the image border.
% 
% Below is a diagram of the pixels adjoining the pixel (m,n):
% 
%            (m,n)------(m,n+1)
%              |          |
%              |          |
%              |          |
%           (m+1,n)---(m+1,n+1)
%
%
% Modifications:
%     December 12, 2013: Added capability to detect and avoid image border.
%     April 25, 2014: Changed (line 57) to count residues using nnz().
%

function [residues, num_residues] = LocateResidues(wrapped_phase, border)

matrix_size = size(wrapped_phase);
% nx = -(matrix_size(1)-1)/2:(matrix_size(1)-1)/2;
% ny = -(matrix_size(2)-1)/2:(matrix_size(2)-1)/2;

residues = zeros(matrix_size);
sum_wrapped_phase_diffs = zeros(matrix_size);

for n = 1:(matrix_size(2)-1)
    for m = 1:(matrix_size(1)-1)
        %if (m,n) or the neighboring pixels as shown above are border
        %pixels, do not evaluate the gradients
        if (border(m,n) == 1 || border(m+1,n) == 1 || ...
                border(m,n+1) == 1 || border(m+1,n+1) == 1)
            continue;
        end  %if border pixel
        
        delta1 = Wrap(wrapped_phase(m+1,n) - wrapped_phase(m,n));
        delta2 = Wrap(wrapped_phase(m+1,n+1) - wrapped_phase(m+1,n));
        delta3 = Wrap(wrapped_phase(m,n+1) - wrapped_phase(m+1,n+1));
        delta4 = Wrap(wrapped_phase(m,n) - wrapped_phase(m,n+1));
        
        sum_wrapped_phase_diffs(m,n) = delta1 + delta2 + delta3 + delta4;
        
        if (sum_wrapped_phase_diffs(m,n) >= 6)  %+2pi
            residues(m,n) = 1;  %Positive residue is denoted +1
        end
        if (sum_wrapped_phase_diffs(m,n) <= -6)  %-2pi
            residues(m,n) = -1;  %Negative residue is dentoed -1
        end
    end  %for m
end  %for n

num_residues = nnz(residues);  %count the phase residues