clear all
close all

bInf = -10.0;
bSup = 10.0;

file_height = sprintf('%d_height.txt',1);
file_discharge = sprintf('%d_discharge.txt',1);
time = 'time.txt';
topo = 'topo.txt';

h = load(file_height);
hu = load(file_discharge);
T = load(time);
Z = load(topo);

% We must print H = h + Z
h = h + Z;
nbCell = length(h);
dx = (bSup-bInf)/(nbCell+1);

X = ones(nbCell,1);

X(1) = bInf;
for i=2:nbCell-1
	X(i) = (bInf + dx/2.0) + (i-1)*dx;
end
X(nbCell) = bSup;

for i=1:length(T)
    titreGraph = sprintf('T = %f',T(i));
    
    subplot(2,1,1); 
    plot(X,h)
    hold on
    plot(X,Z)
    hold off
    axis([X(1) X(nbCell) -0.001 1.0])
    title(titreGraph);
    hleg1 = legend('Height');
    
    subplot(2,1,2);
    plot(X,hu)
    axis([X(1) X(nbCell) -0.5 0.5])
    hleg2 = legend('Discharge');
    pause(0.2)
    
	file_height = sprintf('%d_height.txt',i);
    file_discharge = sprintf('%d_discharge.txt',i);
    h = load(file_height);
    hu = load(file_discharge);
    h = h + Z;
end
