%% Script to read in the NCRMP sampling grids, apply a new grid at desired
%   habitat resolution that is lined up with the NCRMP grids, and from
%   there pare down to just polygon squares that contain some amount of
%   reef within them. this is essentially automating what I did in GIS, and
%   I may not have to actually recreate ALL of it
% 12 June 2024

% oh also - Haibo's FVCOM notes its Cartesian projection and CRS (thank
% god). Double check if Sonaljit's did that, may be yet another thing to
% deal with and email him about but we'll see. as I make grids, can
% continue to think about if they can/should line up explicitly with the
% actual ROMS & FVCOM grids. maybe think about emailing Dubravko for
% particle tracking scripts, etc.