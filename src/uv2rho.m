function [ arrout ] = uv2rho( arr, type )
%   Function to map roms u and v grid to rho grid 
%   usage: 
%        arrout=uv2rho(arr,type)
%              arr   :  u or v array
%              type  :  valid type are 'u' or 'v'
%   Note: data should be in the form F(lat,lon,depth) or F(lat,lon) 
%  
%   Author: vihang Bhatt
%   email: bhatt.vihang@gmail.com

if nargin < 2
    error (" usage: uv2rho(uarr,'u')")
end

switch type
    case 'u'
        if length(size(arr))==3 %3d variable
            arrout=0.5*(arr(:,1:end-1,:)+arr(:,2:end,:));
            arrout=cat(1,arrout(1,:,:),arrout);
            arrout=cat(2,arrout(:,1,:),arrout);
        else  %2D variable
            arrout=0.5*(arr(:,1:end-1)+arr(:,2:end));
            arrout=cat(1,arrout(1,:),arrout);
            arrout=cat(2,arrout(:,1),arrout);
        end
        
    case 'v'
        if length(size(arr))==3 %3d variable
            arrout=0.5*(arr(1:end-1,:,:)+arr(2:end,:,:));
            arrout=cat(1,arrout(1,:,:),arrout);
            arrout=cat(2,arrout(:,1,:),arrout);
        else %2D variable
            arrout=0.5*(arr(1:end-1,:)+arr(2:end,:));
            arrout=cat(1,arrout(1,:),arrout);
            arrout=cat(2,arrout(:,1),arrout);
        end
end


end

