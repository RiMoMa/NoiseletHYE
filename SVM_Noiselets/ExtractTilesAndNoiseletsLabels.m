function [pila_patches_orig,coords,lienzo,pila_patches_gray,labels] = ExtractTilesAndNoiseletsLabels(win_size,ime,aa,bb,GroundT)
    cell_colums = floor(aa/win_size);
    cell_rows = floor(bb/win_size);

    FuturesPatches = rand(win_size,win_size);
    oneTransform = fnt2d(FuturesPatches);
    pila_patches_orig = zeros(cell_rows*cell_colums,size(oneTransform(:)',2) );
    pila_patches_gray = zeros(cell_rows*cell_colums,win_size*win_size );
    lienzo = zeros(size(ime));
    t = 0; %Contador t, que permitira Apilacion...
    coords = zeros(cell_rows*cell_colums,2);   
    labels = zeros(cell_rows*cell_colums,1);    

        % se definen rangos de Corte(Ventaneo) de filas y columnas para sacar los Parches..
    for i = 1:cell_colums
        for j = 1:cell_rows
            Patches_orig = imcrop(ime,[(j-1)*win_size+1 (i-1)*win_size+1 win_size-1 win_size-1]);         %se recortan PARCHES  
            Patches_Mask = imcrop(GroundT,[(j-1)*win_size+1 (i-1)*win_size+1 win_size-1 win_size-1]);         %se recortan PARCHES Mascara                [a,b,p] = size(Patches_orig); %se especifican los tamaños en filas(m) y columnas(n) que deben tener los parches para no ser Descartados..
            [a,b,p] = size(Patches_orig); %se especifican los tamaños en filas(m) y columnas(n) que deben tener los parches para no ser Descartados..

            if a==win_size && b==win_size %APILAR TODOS LOS PARCHES DE BG...
                
                t = t+1;
                pila_patches_gray(t,:) = reshape(Patches_orig,[1,win_size*win_size]); 
                NoiseletTransform = fnt2d(Patches_orig);
                %NoiseletTransform = featureCurvelet(Patches_orig);
                pila_patches_orig(t,:) = NoiseletTransform(:)';
                %% lienzo final
                lienzo(1+win_size*(i-1):win_size*i,1+win_size*(j-1):win_size*j)=ones(win_size,win_size);
                coords(t,:) = [i,j];
                labels(t,:) = any(Patches_Mask(:));
                
            end              
        end
    end
pila_patches_orig = pila_patches_orig(1:t,:);

end