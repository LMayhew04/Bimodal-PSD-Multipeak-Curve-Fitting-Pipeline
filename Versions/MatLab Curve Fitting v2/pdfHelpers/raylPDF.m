function y = raylPDF(x,b,gain)
    if nargin<3, gain = 1; end
    y = gain * (x./b.^2) .* exp(-(x.^2)./(2*b.^2));
end