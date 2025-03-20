# Taller 4
library(geoR)

# cargamos los datos
data(camg)

# dejamos solo el atributo de interes que es mg020 en formato de geodata
mg020 <- as.geodata(camg,data.col = 6)
class(mg020)

# graficamos para ver tendencias y posible boxcox
plot(mg020)
plot(mg020,scatter3d = TRUE)
plot(mg020,lowess = TRUE)

# no se ven tendencias respecto a los valores de magnesio referente a las coordenadas
# la distribucion del magnesio es simetrica, por ende no es necesario corregirla con boxcox

# usaremos de rango la distancia maxima divida en dos
hmax= max(dist(mg020$coords))

# variograma
vg1 <- variog(mg020,trend = "cte",max.dist = hmax/2)

# graficamos el variograma
plot(vg1,max.dist = hmax/2,pch=20)

# vamos a proponer 3 modelos para el variograma
fit1_visual <- eyefit(vg1)

fit1_visual


# luego del ajuste visual validamos recalculamos los 3 modelos
# proponer modelo teorico (maxima verosimilitud)

lik1 <- likfit(mg020,
               ini.cov.pars=fit1_visual[[1]]$cov.pars,
               fix.nugget = TRUE,
               nugget = fit1_visual[[1]]$nugget, 
               cov.model=fit1_visual[[1]]$cov.model,
               lik.method = "ML")


lik2 <- likfit(mg020,
               ini.cov.pars=fit1_visual[[2]]$cov.pars,
               fix.nugget = TRUE,
               nugget = fit1_visual[[2]]$nugget, 
               cov.model=fit1_visual[[2]]$cov.model,
               lik.method = "ML")

lik3 <- likfit(mg020,
               ini.cov.pars=fit1_visual[[3]]$cov.pars,
               fix.nugget = TRUE,
               nugget = fit1_visual[[3]]$nugget, 
               cov.model=fit1_visual[[3]]$cov.model,
               lik.method = "ML")


AIC(lik1)
AIC(lik2)
AIC(lik3)

lik1$cov.model
lik2$cov.model
lik3$cov.model

# grafico consolidando toda la informacion
plot(vg1,pch=20,ylim=c(0,85))
lines(lik1,col=1,lty=1)
lines(lik2,col=2,lty=2)
lines(lik3,col=3,lty=3)

# leyenda para el grafico
legend("topright", 
       legend = c("Exponencial", "Gaussiano", "Cubico"), 
       col = c(1, 2, 3), 
       lty = c(1, 2, 3))


# si la media no es constante, sino lineal

vg_lineal <- variog(mg020,trend = "1st",max.dist = hmax/2)
plot(vg_lineal,pch=20,max.dist = hmax/2)

# ahora hacemos el ajuste para estos casos
fit2_visual <- eyefit(vg_lineal)

# ajustes y graficas
lik_l_1 <- likfit(mg020,
               ini.cov.pars=fit2_visual[[1]]$cov.pars,
               fix.nugget = TRUE,
               nugget = fit2_visual[[1]]$nugget, 
               cov.model=fit2_visual[[1]]$cov.model,
               lik.method = "ML")


lik_l_2 <- likfit(mg020,
               ini.cov.pars=fit2_visual[[2]]$cov.pars,
               fix.nugget = TRUE,
               nugget = fit2_visual[[2]]$nugget, 
               cov.model=fit2_visual[[2]]$cov.model,
               lik.method = "ML")

lik_l_3 <- likfit(mg020,
               ini.cov.pars=fit2_visual[[3]]$cov.pars,
               fix.nugget = TRUE,
               nugget = fit2_visual[[3]]$nugget, 
               cov.model=fit2_visual[[3]]$cov.model,
               lik.method = "ML")


AIC(lik_l_1)
AIC(lik_l_2)
AIC(lik_l_3)

# grafico consolidando toda la informacion
plot(vg_lineal,pch=20,ylim=c(0,85))
lines(lik_l_1,col=1,lty=1)
lines(lik_l_2,col=2,lty=2)
lines(lik_l_3,col=3,lty=3)

# leyenda para el grafico
legend("topright", 
       legend = c("Exponencial", "Gaussiano", "Cubico"), 
       col = c(1, 2, 3), 
       lty = c(1, 2, 3))


# el modelo escogido corresponde al siguiente
lik1 # exponencial con media constante

# creamos la grilla para la prediccion
n = 50
grilla_mg <- expand.grid(east=seq(4957,5961,l=n),
                         north=seq(4829,5710,l=n))


plot(grilla_mg,pch=20,col="gray70")
points(mg020$coords, col=2, pch=20)

# kraging ordinario

okc <- krige.conv(mg020,
                  locations=grilla_mg,
                  krige=krige.control(obj.model=lik1,trend.d = "cte",trend.l = "cte"))


image(okc,val=okc$predict,main="Predicciones para el valor de magnesio")
contour(okc,add = TRUE,main="Predicciones para el valor de magnesio")

image(okc,val=okc$krige.var,main="Varianza de las predicciones")
contour(okc,val=okc$krige.var,filled = TRUE,main="Varianza de las predicciones")
