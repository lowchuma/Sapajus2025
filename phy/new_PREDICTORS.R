# Instalar pacotes necessários
install.packages("ncdf4")
library(raster)
library(rgdal)
library(writexl)
library(ncdf4)

# Definir o diretório de trabalho
setwd("C:/Users/louch/OneDrive/Pasta da Faculdade/LabMastozoo/Lourenço _ 𝘗𝘳𝘪𝘮𝘢𝘵𝘦𝘴/Data_Chuma/Phylogenetics/Sapajus_10000_all_spp")

# Carregar a planilha de localidades
cs_size_sapajus <- read.csv("67indv_data.csv")
species <- as.factor(cs_size_sapajus$sp)
summary(species)

# Carregar dados bioclimáticos do WorldClim
bioclim_dir <- "C:/Users/louch/OneDrive/Pasta da Faculdade/LabMastozoo/Lourenço _ 𝘗𝘳𝘪𝘮𝘢𝘵𝘦𝘴/Data_Chuma/Bioclim/wc2.1_2.5m_bio"  # Diretório onde os arquivos bioclim estão extraídos
bioclim_files <- list.files(bioclim_dir, pattern = ".tif$", full.names = TRUE)

# Carregar os 19 arquivos raster bioclimáticos
bioclim_stack <- stack(bioclim_files)
names(bioclim_stack) <- paste0("BIO", 1:19)
print(bioclim_stack)

# Função para carregar dados netCDF
load_nc_data <- function(file_path, var_name) {
  nc <- nc_open(file_path)
  data <- ncvar_get(nc, var_name)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  nc_close(nc)
  raster_data <- raster(data, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +datum=WGS84"))
  return(raster_data)
}

# Caminhos dos arquivos netCDF de NPP e AVGANNRH
npp_file_path <- "C:/Users/louch/OneDrive/Pasta da Faculdade/LabMastozoo/Lourenço _ 𝘗𝘳𝘪𝘮𝘢𝘵𝘦𝘴/Data_Chuma/Bioclim/NPP/npp"
avgannrh_file_path <- "C:/Users/louch/OneDrive/Pasta da Faculdade/LabMastozoo/Lourenço _ 𝘗𝘳𝘪𝘮𝘢𝘵𝘦𝘴/Data_Chuma/Bioclim/AVGANNRH/avgannrh"

# Carregar dados de NPP e AVGANNRH
npp_data <- list.files(npp_file_path, pattern = ".adf", full.names = TRUE)
avgannrh_data <- list.files(avgannrh_file_path, pattern = ".adf", full.names = TRUE)

# Carregar dados de NPP e AVGANNRH
npp_data <- raster(npp_file_path)
avgannrh_data <- raster(avgannrh_file_path)

# Função para extrair dados para cada espécie
extract_climate_data <- function(data, species_name) {
  species_data <- subset(data, sp == species_name)
  longitude <- species_data$long
  latitude <- species_data$lat
  coordinates <- data.frame(x = longitude, y = latitude)
  points <- SpatialPoints(coordinates, proj4string = CRS(projection(bioclim_stack)))
  
  for (i in 1:19) {
    var_name <- paste0("BIO", i)
    species_data[[var_name]] <- extract(bioclim_stack[[var_name]], points)
  }
  
  species_data$NPP <- extract(npp_data, points)
  species_data$AVGANNRH <- extract(avgannrh_data, points)
  
  return(species_data)
}

# Lista de espécies
species_list <- unique(cs_size_sapajus$sp)

# Extrair dados para todas as espécies
all_species_data <- do.call(rbind, lapply(species_list, function(sp) extract_climate_data(cs_size_sapajus, sp)))

# Verificar os dados combinados
head(all_species_data)
print(all_species_data$NPP)
all_species_data <- all_species_data[, -c(2-60)]
all_species_data$NPP[is.na(all_species_data$NPP)] <- 0

# Salvar os dados em um arquivo Excel
write_xlsx(all_species_data, "NPP.xlsx")




