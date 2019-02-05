function [ err1, err2 ] = CalcErr(beam1, beamref)
%CalcErr Calculate an errorsignal for an X photodiode
%   Create the mask to isolate the individual quadrants, then sum/subtract
[row,col] = size(beam1);
mask1 = zeros(row,col);
for n=1:row
    for m=1:col
        if n>m
            mask1(n,m) = 1;
        end
    end
end
mask2 = rot90(mask1);

mask_upper = mask1.*mask2;
mask_right = rot90(mask_upper);
mask_lower = rot90(mask_right);
mask_left = rot90(mask_lower);

pwr_upper = sum(sum(mask_upper.*imag(beam1.*conj(beamref))));
pwr_right = sum(sum(mask_right.*imag(beam1.*conj(beamref))));
pwr_lower = sum(sum(mask_lower.*imag(beam1.*conj(beamref))));
pwr_left = sum(sum(mask_left.*imag(beam1.*conj(beamref))));

err1 = pwr_left + pwr_right - pwr_upper - pwr_lower;
err2 = -pwr_left + pwr_right + pwr_upper - pwr_lower;
end

