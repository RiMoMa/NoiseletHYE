% Encontrar la fila con información en la celda
fila_con_info = find(~cellfun(@isempty, celda(:,1)), 1, 'first');

% Extraer los datos de las columnas correspondientes en la celda
datos = celda(fila_con_info, :);
valores = [datos{2:2:end}];

% Calcular la media de los valores
media = mean(valores);
fprintf('%d\n',media)
