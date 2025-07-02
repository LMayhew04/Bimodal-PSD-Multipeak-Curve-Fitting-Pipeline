function y = lognPDF(x, mu, sigma, A)
%LOGNPDF  Lognormal PDF with amplitude scaling.
y = A * exp(-0.5 * ((log(x)-log(mu))/sigma).^2) ./ (x * sigma * sqrt(2*pi));
y(x <= 0) = 0;
end