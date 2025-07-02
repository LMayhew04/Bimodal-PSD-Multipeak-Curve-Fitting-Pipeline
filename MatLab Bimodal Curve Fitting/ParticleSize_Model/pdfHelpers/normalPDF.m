function y = normalPDF(x, mu, sigma, A)
%NORMALPDF  Normal distribution PDF with amplitude scaling.
y = A * exp(-0.5 * ((x-mu)/sigma).^2) ./ (sigma * sqrt(2*pi));
end