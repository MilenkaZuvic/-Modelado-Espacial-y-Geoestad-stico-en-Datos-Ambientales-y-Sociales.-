rm(list = ls())
gc()

# Datos de área - London Suicides -----------------------------
# LDNSuicides.shp (la geometría)
# LDNSuicides.shx (índice espacial)
# LDNSuicides.dbf (atributos)

## Datos de suicidios en 32 municipios de Londres
## LDNSuicides.shp: información espacial de los municipios (shapefile/archivo de forma)
## LondonSuicides.RData: datos/variables relevantes 
## el número de suicidios observados (y), 
## casos esperados (E), 
## los índices de privación (x1)
## fragmentación social (x2)


# Pregunta 2: Datos de área - London Suicides -----------------------------

# Cargar librerías necesarias
library(sf)         # shapefiles
library(spdep)      # análisis espacial y matrices de pesos
library(spatialreg) # Para modelos espaciales

# Cargar el shapefile
london_shape <- st_read("LDNSuicides.shp")
str(london_shape)
names(london_shape)

# Cargar los datos de suicidios
load("LondonSuicides.RData")
ls()

# Crear un data frame con las variables cargadas
LondonSuicides <- data.frame(y = y, E = E, x1 = x1, x2 = x2) # datos de los suicidios
summary(LondonSuicides) 

# Vericar numero de filas
nrow(LondonSuicides)
nrow(london_shape)

# Unir y Añadir las columnas de los datos LondonSuicides a london_shape (polígono espacial)
london_shape$y <- LondonSuicides$y
london_shape$E <- LondonSuicides$E
london_shape$x1 <- LondonSuicides$x1
london_shape$x2 <- LondonSuicides$x2

str(london_shape)



# Crear variable SMR ---------------------------------------------------------------------

## 1. Aplique los test de Moran y Geary para la variable SMR. 
## 2. Utilice la matriz W estandarizada por fila. 
## 3. Se justifica un modelo espacial? 

## La razon de mortalidad estandarizada se calcula como 
## RME = muertes observadas (y) /muertes esperadas (E)

london_shape$SMR <- london_shape$y / london_shape$E


# Crear la matriz de pesos W  -------------------------------------------

library(spdep)
adyacencia <- poly2nb(london_shape)

# Matriz (W) estandarizada
W <- nb2listw(adyacencia, style = "W", zero.policy = TRUE) ## Utilzar Modelo SAR

# Matriz binaria (B)
W_binaria <- nb2listw(adyacencia, style = "B", zero.policy = TRUE) ## Utilizar Modelo CAR



# Test de Moran y Geary -----------------------------------------------------------

# Test de Moran
test_moran <- moran.test(london_shape$SMR, W)

# Test de Geary
test_geary <- geary.test(london_shape$SMR, W)

# Resultados
test_moran
test_geary


#  Ajuste de modelos para SMR  -------------------------------

# Modelo 1: Lineal   ------------------------------------------------------

model_1 <- lm(SMR ~ x1 + x2 + I(x1*x2), data = london_shape)
summary(model_1)


# Modelo 2: SAR ----------------------------------------------------------------
model_2 <- errorsarlm(SMR ~ x1 + x2 + I(x1*x2), data = london_shape, listw = W)
summary(model_2)


# Modelo 3: CAR  ----------------------------------------------------------------
model_3 <- spautolm(SMR ~ x1 + x2 + I(x1*x2), data = london_shape, listw = W_binaria, family = "CAR")
summary(model_3)



# Comparacion AIC ---------------------------------------------------------

aic_model1 <- AIC(model_1)
aic_model2 <- AIC(model_2)
aic_model3 <- AIC(model_3)

aic_model1 #LINEAL
aic_model2 #SAR  
aic_model3 #CAR


# Modelo Seleccionado -----------------------------------------------------

# Modelo LM ---------------------------------------------------------------

model_car <- spautolm(SMR ~ x1 + x2, data = london_shape, listw = W_binaria, family = "CAR")
summary(model_car)

aic_model_car <- AIC(model_car)
aic_model_car


# Graficar Residuos ------------------------------------------------------

library(tmap)

## Model selecc: model_1
## Guardar el ajuste y los residuos

london_shape$SMR_est <- fitted(model_3) # modelo final seleccionado 
london_shape$residuals <- residuals(model_3)

mapa_srm_real <- tm_shape(london_shape) + tm_polygons("SMR", title = "SMR Real") + 
  tm_layout(title = "SMR Real en Municipios de Londres") # Graficar SMR real
tmap_save(mapa_srm_real, "SMR_mapa.png")

mapa_srm_estimado <- tm_shape(london_shape) + tm_polygons("SMR_est", title = "SMR Estimado") + 
  tm_layout(title = "SMR Estimado en Municipios de Londres") # Graficar SMR estimado 
tmap_save(mapa_srm_estimado, "SMR_mapa_estimado.png")

moran_resid <- moran.test(london_shape$residuals, W) 
geary_resid <- geary.test(london_shape$residuals, W) 

print(moran_resid)
print(geary_resid)

mapa_residuos_car <- tm_shape(london_shape) + tm_polygons("residuals", title = "Residuos del Modelo Final") + 
  tm_layout(title = "Mapa de Residuos Modelo CAR")
tmap_save(mapa_residuos_car, "SMR_mapa_residuos_car.png")

london_shape$residuals_ln <- residuals(model_1)
mapa_residuos_ln<- tm_shape(london_shape) + tm_polygons("residuals_ln", title = "Residuos del Modelo Final") + 
  tm_layout(title = "Mapa de Residuos Modelo Lineal")
tmap_save(mapa_residuos_ln, "SMR_mapa_residuos_ln.png")

mapa_combinado <- tmap_arrange(mapa_residuos_car, mapa_residuos_ln, ncol = 2)

# Guardar el mapa combinado
tmap_save(mapa_combinado, "mapa_combinado.png", width = 12, height = 5)


#   -----------------------------------------------------------------------


