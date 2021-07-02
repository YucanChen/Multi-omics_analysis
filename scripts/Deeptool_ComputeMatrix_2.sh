#!/bin/bash

cd ./Deeptools_Results
plotProfile -m Combined_matrix_regionclass.gz --refPointLabel "0" --samplesLabel "WT1" "WT2" "Setd2KO" --yMax 1.5 -out ./Combined_matrix_regionclass_row.png --numPlotsPerRow 2
plotProfile -m Combined_matrix_regionclass.gz --refPointLabel "0" --samplesLabel "WT1" "WT2" "Setd2KO" --yMax 1.5 -out ./Combined_matrix_regionclass_fill.png --plotType=fill --perGroup --colors blue purple red
plotProfile -m Combined_matrix_regionclass.gz --refPointLabel "0" --samplesLabel "WT1" "WT2" "Setd2KO" --yMax 1.5 -out ./Combined_matrix_regionclass.png --plotType=lines --perGroup --colors blue purple red
plotHeatmap -m Combined_matrix_sampleclass.gz -out Combined_matrix_sampleclass.png --colorMap RdYlBu --whatToShow 'heatmap and colorbar' --zMin -4 --zMax 4 --dpi 100
exit 0
