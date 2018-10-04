% 2018/10/02 Imported from scilab.

function [noi] = mean_chisq_v2(data, df, lnseg)
    % MEAN_CHISQ  Mean of a set assuming Chi-Square dist.
    %
    %  data is flattened along the first dimension.
    %  
    %  MEAN_CHISQ(data, df, lnseg) calculates the noise floor level of the
    %  given set.
    %
    %  data : A set of random variables following the Chi-Square dist.
    %  df   : The degrees of freedom in the given set.
    %         Default is 2 (assuming real and imag).
    %  lnseg: Length of each block into which the given set is divided.
    %         Default is 32.
    %
    % Example:
    %
    %   mean_chisq(data, 2)
    
    if nargin < 3
        lnseg = 32;
    end
    if nargin < 2
        df = 2;
    end
    
    % The number of blocks.
    ndiv = fix(numel(data) / lnseg);
    if ndiv == 0
        error("lnseg (%d) must be shorter than %d.", lnseg, numel(data))
    end
    % Blockwise minimum.
    bmin = nanmin(nanmean(reshape(data(1:lnseg*ndiv), lnseg, ndiv), 1));
    % Correction factor: gaussian to chi-square.
    a = sqrt(2) * erfcinv(2 / ndiv);
%     m = ndiv / sqrt(2*pi)*exp(-0.5*a^2);
    c = 1 / (1 - a / sqrt(lnseg*df));
    noi = c * bmin;
end