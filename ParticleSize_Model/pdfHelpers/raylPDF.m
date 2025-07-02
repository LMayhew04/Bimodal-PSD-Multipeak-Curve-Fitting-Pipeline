function y = raylPDF(x, sigma, A)
%RAYLPDF  Rayleigh PDF with amplitude scaling.
y = A * (x/sigma^2) .* exp(-x.^2/(2*sigma^2));
end