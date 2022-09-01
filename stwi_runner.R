library(raster)

# Chargement des fonctions
rstudioapi::getSourceEditorContext()$path |>
  dirname() |>
  file.path("stwi_functions.R") |>
  source()


# Matrice obligatoire (modèle numérique de terrain)
wdir <- "Y:/Developpement/Programmation/JFB/Saga_TWI"
path_dem <- file.path(wdir, "dem.tif")
nodata <- -9999  # valeur réelle de NoData du DEM
cellsize <- 1  # en unité

# Matrice facultative (matrice de poids)
path_weight <- file.path(wdir, "weight.tif")


# Paramétrage de l'algorithme
SUCTION <- 10
AREA_TYPE <- 2    # 0 = "total catchment area" | 1 = "square root of catchment area" | 2 = "specific catchment area"
SLOPE_TYPE <- 1   # 0 = "local slope" | 1 = "catchment slope"
SLOPE_MIN <- 0    # en degrés
SLOPE_OFF <- 1    # en degrés
SLOPE_WEIGHT <- 1


# Production à partir de SAGA GIS
# https://saga-gis.sourceforge.io/saga_tool_doc/7.0.0/ta_hydrology_15.html
saga_exe <- "//ulysse/LIDAR/Developpement/Hydrographie/Outils/Utilitaires/saga-7.3.0_x64/saga_cmd.exe"
suffixe <- paste0("area", AREA_TYPE, "_slopeType", SLOPE_TYPE)
path_TWI_SDAT <- file.path(wdir, paste0("TWI_saga_", suffixe, ".sdat"))
path_AREA_SDAT <- file.path(wdir, paste0("AREA_saga_", suffixe, ".sdat"))
path_AREA_MOD_SDAT <- file.path(wdir, paste0("AREAMOD_saga_", suffixe, ".sdat"))
path_SLOPE_SDAT <- file.path(wdir, "SLOPECATCHMENT_saga.sdat")

cmd <- paste(
  saga_exe,
  "ta_hydrology 15",
  "-DEM", path_dem,
  "-WEIGHT", path_weight,
  "-TWI", path_TWI_SDAT,
  "-AREA", path_AREA_SDAT,
  "-AREA_MOD", path_AREA_MOD_SDAT,
  "-SLOPE", path_SLOPE_SDAT,
  "-SUCTION", SUCTION,
  "-AREA_TYPE", AREA_TYPE,
  "-SLOPE_TYPE", SLOPE_TYPE,
  "-SLOPE_MIN", SLOPE_MIN,
  "-SLOPE_OFF", SLOPE_OFF,
  "-SLOPE_WEIGHT", SLOPE_WEIGHT
)
system(cmd)



# Production à partir de R
rDEM <- raster(path_dem)
rWeight <- raster(path_weight)
outpath_twi <- file.path(wdir, paste0("TWI_r_area", AREA_TYPE, "_slopeType", SLOPE_TYPE, ".tif"))

rTWI <- Get_STWI(rDEM, SUCTION, AREA_TYPE, SLOPE_TYPE, SLOPE_WEIGHT, SLOPE_MIN, SLOPE_OFF, nodata, cellsize, rWeight)
raster::writeRaster(rTWI, outpath_twi)
