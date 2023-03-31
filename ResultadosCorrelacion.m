
pathFiles = '/mnt/md0/ricardo/NoiseletProject/ResultsNoiselets_Original_Complete_noOverlap_Monuseg_ecDistArtic/*review.mat';


ResultFiles =dir(pathFiles);
Resultados={};
for x=1:length(ResultFiles)
    fileName = ResultFiles(x).name;
    Tejido = find(fileName=='_');
    NombreTejido = fileName(Tejido(1)+1:Tejido(2)-1)
    load(fileName);

celda1=Resultados_Nucleos;
celda2=Resultados_NucleosMetodo;
% Encontrar la fila con información en celda1
fila_con_info = find(~cellfun(@isempty, celda1(:,1)), 1, 'first');
% Extraer los datos de las columnas correspondientes en celda1
datos1 = celda1(fila_con_info, :);
nombres1 = datos1(1:2:end);
valores1 = [datos1{2:2:end}];
% Extraer los datos de las columnas correspondientes en celda2
datos2 = celda2(fila_con_info, :);
nombres2 = datos2(1:2:end);
valores2 = [datos2{2:2:end}];
% Calcular la correlación entre los valores de las columnas correspondientes
correlacion = corr(valores1', valores2');
Nimagenes = length(valores1);
Resultados(x,:)={NombreTejido,correlacion,Nimagenes};
end