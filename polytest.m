% Read file 
A = dlmread('dat/R_020_03.dat');

N=20;
t = size(A,1)/N;

monomers = [1:N];

%F(t) = struct('cdata',[],'colormap',[]);
close all
f = figure;
%F(1:length(1:100:t)) = struct('cdata',[], 'colormap',[]);
R2 = zeros(1, t);
for i=1:1:t
    X =  A((i-1)*N+monomers, 1);
    Y =  A((i-1)*N+monomers, 2);
    Z =  A((i-1)*N+monomers, 3);
    
    COM = [mean(X),mean(Y),mean(Z)];
    
    X_ = X - COM(1);
    Y_ = Y - COM(2);
    Z_ = Z - COM(3);
    
    R2(i) = mean( X_.^2 + Y_.^2 + Z_.^2 );
    
    %plot3( X, Y, Z, 'r' )
    %hold on
    %scatter3( X, Y, Z,'filled')
    
    
%scatter3( mean(X),mean(Y),mean(Z),'filled','r' )
 %   hold off
    


    axis equal
    %xlim([COM(1)-10 COM(1)+10])
    %ylim([COM(2)-10 COM(2)+10])
    %zlim([COM(3)-10 COM(3)+10])
    xlim([-5 10]);
    ylim([-5 10]);
    zlim([-5 10]);
    drawnow
    
    %k=waitforbuttonpress;
    %F(i) = getframe(gcf);
    
end

plot(sqrt(R2));

%movie2avi(F,'poly.avi');


