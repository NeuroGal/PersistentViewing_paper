function limits = limer(array, round_c)
% This code was written as part of visualization of Vishne et al., biorxiv 2022 https://doi.org/10.1101/2022.08.02.502469
%   'Representing experience over time: sustained sensory patterns and transient frontoparietal patterns'
%   So please cite (-:
% Written by Gal Vishne, bug reports \ requests: gal.vishne@gmail.com
if ~exist('round_c','var'); round_c = 10^(-floor(log10(abs(mean(array(:)))))); end
limits = [floor(min(round_c*array(:))), ceil(max(round_c*array(:)))]/round_c;
end