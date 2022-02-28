#Volcano Plot datos GSE188915(Comparison of organoids)
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
##Comenzar el Volcano Plot

res<- results(dds, contrast<- c("Muestra","Hormone_treated_Endometrial_Organoid", "Endometrial_Organoid"))
res<- res[order(res$pvalue),]

#Convertir el objeto S4 con clase de DESeqResults a un data frame para poder graficar
res_volcano<- data.frame(res)

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

res_volcano$symbol<- mapIds(org.Hs.eg.db,
                    keys = row.names(res_volcano),
                    column ="SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")

res_volcano$entrez<- mapIds(org.Hs.eg.db,
                    keys = row.names(res_volcano),
                    column ="ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals = "first")

res_volcano$name<- mapIds(org.Hs.eg.db,
                  keys = row.names(res_volcano),
                  column ="GENENAME",
                  keytype = "ENSEMBL",
                  multiVals = "first")
head(res_volcano, 10)

#Crear una nueva columna DF
res_volcano$category <- "no change"
str(res_volcano)
head(res_volcano)

#Modificar nuestra nueva columna usando valores de corte [] y $ para modificar data frame, aplicar valores de corte para encontrar UP, DOWN y significativo, & es de agregar

res_volcano$category[res_volcano$log2FoldChange > 1.5 & res_volcano$pvalue < 0.05] <- "Up regulated"
res_volcano$category[res_volcano$log2FoldChange < -1.5 & res_volcano$pvalue < 0.05] <- "Down regulated"
View(res_volcano)

#ggplot para hacer volcano plot

library(ggplot2)


volc= ggplot(data=res_volcano, aes(x=log2FoldChange,y= -log(pvalue),colour=category))+
  geom_point()+ ggtitle("ANALISIS DE EXPRESION DIFERENCIAL DE TEJIDO ENDOMETRIAL CON Y SIN TRATAMIENTO HORMONAL")

volc+geom_text(data=head(res_volcano, 10), aes(label=symbol)



