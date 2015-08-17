function [ A ] = read_xyz( fileName, N )

%A=[];

fid = fopen(fileName);
%tline = fgetl(fid); %Get N line

A = fscanf( fid, '%*s\t%f\t%f\t%f', [3,20000] )';
% while ischar(tline)
%     
%     tline = fgetl(fid); % Get comment line
%     
%     for i=1:N
%         tline = fgetl(fid);
%         
%         A = [A; sscanf( tline, '%*s\t%f\t%f\t%f' )'];
%         
%     end
%     tline = fgetl(fid); % Get N line
% 
% end

end

