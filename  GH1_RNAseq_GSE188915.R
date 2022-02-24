#ANALISIS DE EXPRESION DIFERENCIAL
#RNA-seq data analysis con DESeq2 de la Comparación de los  organoides provenientes de  Fluido Menstrual, Endometrio y Endometrio Hormonalmente-tratado "GSE188915"

#Llamar las librerias, DESeq2 ya había sido previamente instalado 
library(dplyr)
library(DESeq2)

#Importar countdata
library(readr)

#Buscar la ruta del archivo csv
file.choose()

#Copiar la ruta en la consola
ruta_csv<- "/Users/montse/Documents/Estancia Estudiantil/TRABAJOS/GSE188915(Comparison of organoids)/GSE188915_Organoid_TMM_Normalised_MTremoved_CPM_Counts_ENSEMBL.csv.gz"

#Importar datos
countData<-read.csv(ruta_csv, row.names = 1) %>% 
  as.matrix()
#Ordenar columnas
countData [, c(5,6,2,4,7,1,3)]
countDataordenado<- countData [, c(5,6,2,4,7,1,3)]

#Filtrar los datos donde solo se tenga 0 o 1 en todas las lecturas de las muestras
countData= countDataordenado[rowSums(countData) >1, ]
head(countData)

#Importar metadata 

#Importar paquete para leer archivos de excel
#Llamar a la biblioteca
library(readxl)

#Buscar el archivo de excel
file.choose()

#Ruta del archivo de excel
ruta_excel<- "/Users/montse/Documents/Estancia Estudiantil/TRABAJOS/GSE188915(Comparison of organoids)/MetadataGSE188915.xlsx"


#Importar desde excel
colData<- read_excel(ruta_excel)
colData

#Configurar el objeto DESeqDataSet y correr el pipeline de DESeq
dds <- DESeqDataSetFromMatrix(countData=round(countData), 
                                  colData=colData, 
                                  design=~Muestra)
dds<- DESeq(dds)
dds

##Ahora obtener los resultados de la expresión de genes de Endometrial Organoid vs. Hormone treated endometrial organoid, y reordenarlos por valor. 
##Llamar a la función "summary" en los resultados del objeto para obtener noción de cuantos genes estan up o down-regulated en FDR 0.1.

res<- results(dds, contrast<- c("Muestra","Hormone_treated_Endometrial_Organoid", "Endometrial_Organoid"))
res<- res[order(res$pvalue),]
summary(res)

##Ahora obtener los resultados de del knockdown de Menstrual_Fluid_Organoid vs. Hormone treated endometrial organoid, y reordenarlos por valor. 
res2<- results(dds, contrast<- c("Muestra","Hormone_treated_Endometrial_Organoid", "Menstrual_Fluid_Organoid"))
res2<- res2[order(res2$pvalue),]
summary(res2)

#Como solo realizamos un mapeo y cuenta vs. anotacion de Ensembl, nuestros resultados solo tienen informacion acerca de los IDs de genes de Esembl.
#La forma canocida del Bioconductor para hacer esto es con la AnnotationDbi y paquetes de anotación del organismo. 
#Aqui estamos usando organism package ("org") para Homo sapiens ("Hs"), organizado como Annotation Dbi data base ("db") usando Entrez Gene IDs ("eg") como primary keys.
#Para ver todas las llaves, usar función de columnas

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

#Ahora usemos la funcion mapIds para agregar mas columnas a los resultados
#Los row.names de los resultados de nuestra tabla tienen el geneID de Esembl, asi que necesitamos especificar key=ENSEMBL.
#El argumento columna le dice  a la funcion mapIdz que informacion queremos, y el argumento multiVals le dice a la funcion que debe hacer si hay múltiples posibles valores para un solo valor de entrada.
#Aqui le pedimos que solo nos regrese la primera que ocurra en el database. 
#Vamos a obtener los IDs de Entrez, simbolos de genes y el nombre completo de los genes.

res$symbol<- mapIds(org.Hs.eg.db,
                    keys = row.names(res),
                    column ="SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")

res$entrez<- mapIds(org.Hs.eg.db,
                    keys = row.names(res),
                    column ="ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals = "first")

res$name<- mapIds(org.Hs.eg.db,
                  keys = row.names(res),
                  column ="GENENAME",
                  keytype = "ENSEMBL",
                  multiVals = "first")
head(res, 10)

#ANALISIS DE VIAS

#Vamos a usar el paquete de gage (Generally Applicable Gene-set Enrichment for Pathway Analysis) para analisis de vias.
#Ver también el gage package workflow vignette para analisis de vias de RNA-seq 
#Una vez tengamos una lista de vias enriquecidas, vamos a usar el paquete pathview para dibujar diagramas de vias, sombreando las moleculas de la via por su grado de expresion a la alta o a la baja (up/down-regulation)

#KEGG pathways
#El paquete gageData ha pre-compilado bases de datos de genes mapeados a KEGG pathways y terminos GO para organismos comunes
#Kegg.sets.hs es una lista nombrada de 299 elementos. Cada elemento es un vector de caracter del gen miembro de ENTREZ IDs para una sola via KEGG
#(tambien ver go.sets.hs). sigmet.idx.hs es un index de numeros de señalizacion y vias metabolicas en kegg.set.gs.
#En otras palabras, KEGG pathway incluye otro tipo de definiciones de vias, como "Global Map" y "Enfermedades Humanas", que podrian ser indeseables en el analisis de vias
#Por tanto, kegg.sets.hs[sigmet.idx.hs] brindan los sets mas limpios de vias solamente de señalizacion y metabolicas.

#Cargar librerias
library(pathview)
library(gage)
library(gageData)

#Subir data
data("kegg.sets.hs")
data("sigmet.idx.hs")

#Nombrar las variables de la informacion subida
kegg.sets.hs<- kegg.sets.hs[sigmet.idx.hs]

#Imprimir el head de la tabla
head(kegg.sets.hs, 3)

#La funcion gage() requiero un verctor nombrado de fold change (veces de cambio), donde los nombres de los valores estan en Entrez gene IDs

foldchanges<- res$log2FoldChange
names(foldchanges)<- res$entrez
head(foldchanges)

#Ahora, vamos a correr el Análisis de vias.
#Para ver ayuda sobre la funcion escribir "?gage"
#Especialmente, tu quizá querras tratar de cambiar el valor de same.dir. Este valor determina si hay cambios en un gen establecido a solo una direccion
#(Todos los genes up o down regulated) o cambios hacia ambas direcciones simultaneamente (algunos genes en la vía mal regulados)
##Para los conjuntos de genes derivados experimentalmente, terminos GO, etc. La corregulacon es un caso comun (samer.dir=TRUE), por default;
##En KEGG, BioCarta pathways, la frecuencia de los genes no estan coreguladas, por lo tanto podria ser informativo dejas (same.dir=FALSE)
##Tambien (same.dir=TRUE) podría ser interesante para las vias

#En este caso se usaran (same.dir=TRUE), nos dara listas separadas para las vias que estan upregulated vs las que estan downregulated
#Veamos los primeros resultados para cada uno

#Obtener resultados
keggres<- gage(foldchanges, gsets=kegg.sets.hs, same.dir = TRUE)

#Ver ambas up(más grandes), abajo (menos) y estadísticas
lapply(keggres, head)

#Ahora vamos a procesar los resultados para extraer las 5 principales vias reguladas hacia arriba, y luego procesarlas mas para obetener los ID
#Se utilizan estos ID de ruta KEGG downstream para graficar

#Obtener rutas (vias)
keggrespathways= data.frame(id=row.names(keggres$greater), keggres$greater) %>% 
  tibble::as_tibble()%>% 
  filter(row_number()<=5) %>% 
  .$id %>%   
  as.character()  
keggrespathways

#Vias menos expresadas
keggrespathwaysless= data.frame(id=row.names(keggres$less), keggres$less) %>% 
  tibble::as_tibble()%>% 
  filter(row_number()<=5) %>% 
  .$id %>%   
  as.character()  
keggrespathwaysless

#Obtencion de los IDs
keggresids= substr(keggrespathways, start= 1, stop=8)
keggresids

#Finalmente, la fucion "pathview ()" en el paquete pathview hace graficos.
#Correr una funcion para recorrer y dibujar graficos para los 5 caminos principales que creamos con anterioridad

#Definir la funcion de graficar para aplicar despues
plot_pathway<- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

#Graficar multiples rutas (plots guardados en el disco y que regresan un objeto de lista desechable )
tmp<- sapply(keggresids, function(pid) pathview(gene.data=foldchanges,pathway.id=pid, species="hsa"))
