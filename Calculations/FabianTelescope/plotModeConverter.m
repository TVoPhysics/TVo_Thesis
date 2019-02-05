function [ ] = plotModeConverter(inFabian, outFabian )
%Takes three images where each element of the matrix may have a complex
%form and plots the complex and real parts seperately.

%Find max
mavoFr = max(max(real(outFabian)));
mavoFi = max(max(imag(outFabian)));
mivoFr = min(min(real(outFabian)));
mivoFi = min(min(imag(outFabian)));
macc = sqrt(max(max(conj(outFabian).*outFabian)));
oFr = real(outFabian);
oFi = imag(outFabian);
noFr = oFr-macc;%(oFr -mivoFr)./(mavoFr - mivoFr).*100;
noFi = oFi-macc;%(oFi -mivoFi)./(mavoFi - mivoFi).*100;

figure

subplot(3,3,1)
imagesc(conj(inFabian).*inFabian);
title('| inFabian | ^2')
colormap('gray')

subplot(3,3,2)
[imY,imZ]=imread('Mode_converter','png');
imshow(imY);

subplot(3,3,3)
imagesc(conj(outFabian).*outFabian);
title('| outFabian | ^2')
colormap('gray')

subplot(3,3,4)
imagesc(inFabian)
title('inFabian')
colormap('gray')

subplot(3,3,5)
[imY,imZ]=imread('Mode_converter','png');
imshow(imY);

subplot(3,3,9)
imagesc(noFi);
title('Imaginary part of outFabian')
colormap('gray')

subplot(3,3,6)
imagesc(noFr)
title('Real part of outFabian')

end

