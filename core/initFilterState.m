function state = initFilterState(birth)
%INITFILTERSTATE Initialize the undetected PPP and an empty detected-object MBM.
state = struct();
state.poissonUndetected = birth;
state.trackBern         = {};
state.globTable         = uint32([]);
state.globLogW          = [];
end
