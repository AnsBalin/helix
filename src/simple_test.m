files = dir('dat');

Rg_N = [];

for f=1:length(files)
    
    if files(f).name(1) =='R';
        disp(files(f).name)
        N = str2double(files(f).name(3:5));
        sim = str2double(files(f).name(7:8));
        
        A = dlmread(strcat('dat/',files(f).name));
        %A = read_xyz( strcat('dat/',files(f).name), N );
        t = size(A,1)/N;
        monomers = [1:N];
        
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
            
            
            
        end
        
        Rg_N = [Rg_N; N,sqrt(mean(R2))];
        
        %plot(sqrt(R2));
        
    end
    
end

Narr = [5, 10, 20, 40, 80, 160];
Rarr = zeros(1,length(Narr));
Rerr = zeros(1,length(Narr));
for n=1:length(Narr)

    ninds = find(Rg_N(:,1)==Narr(n));
    Rarr(n) = mean(Rg_N(ninds,2));
    Rerr(n) = std(Rg_N(ninds,2)/sqrt(length(ninds)));
    
    
end

errorbar(Narr(1:(end-1)),Rarr(1:(end-1)),Rerr(1:(end-1)),'o');
x = 5:160;
y = x.^(1/2);
z = x.^(3/5);
hold on
plot(x,y,x,z);
ylabel('Rg');
xlabel('N');
ax = gca;
ax.YScale='log';
ax.XScale='log';
xlim([3 200])

