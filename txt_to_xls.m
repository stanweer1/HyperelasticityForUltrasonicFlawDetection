clear

%Declare displacement U matrix
U = zeros(1201,65);
A = zeros(1201,1);

%First column for time (edit based on time step)
U(:,1) = [0:5e-6:0.006]';

%Loop to store sensor values in 2 to 65
for i=2:65
    name = ['U3_1_' num2str(i-1)];
    A = dlmread(name);
    U(:,i) = A;
end

file = 'Merged.xls';
dlmwrite(file,U,'delimiter','\t')

%Plot one graph between time and specified sensor
n = input('Enter sensor number (1-64): ');
plot(U(:,1),U(:,n+1));