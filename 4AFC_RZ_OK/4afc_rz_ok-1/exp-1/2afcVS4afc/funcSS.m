function PredictedCritVals = funcSS( params, xs )

%FUNCSS Summary of this function goes here
%   Detailed explanation goes here

%PredictedCritVals = params{1} .* xs.^ params{2} + params{3} .* xs .^ params{4};
PredictedCritVals = params{1} .* xs.^ (-2) + params{2} .* xs .^ params{3};