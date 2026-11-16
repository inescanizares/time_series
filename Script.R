# Trabajo Practico series de tiempo 
# Este archivo es para probar algunas cosas antes de hacer el informe 

# Carga de paquetes
library(readxl)
library(readr)
library(tidyverse)
library(feasts)
library(fpp3)
library(tsibble)
library(lubridate)
library(forecast)
library(rmarkdown)
library(knitr)
library(gridExtra)
library(kableExtra)
library(patchwork)
options(scipen = 999)

# Cargar los datos
datos <- read.delim(
  "Multiple_Cause_of_Death__1999-2020.xls",
  sep = "\t"
)

# Transformar los datos en una serie temporal mensual
datos <- datos |>
  mutate(Month = yearmonth(`Month.Code`)) |>
  select(Month, muertes = Deaths) |>
  as_tsibble(index = Month)

# Se seleccionó el período comprendido entre enero de 1999 y diciembre de 2018.
datos_filtrado <- datos |>
  filter_index("1999 Jan" ~ "2018 Dec")

graf_original <- ggplot(
  data = datos_filtrado,
  aes(x = Month, y = muertes)
) +
  geom_line(
    color = "cadetblue4",
    linewidth = 0.6
  ) +
  labs(
    x = "Fecha",
    y = "Cantidad de muertes"
  ) +
  theme_minimal() +
  scale_x_yearmonth(
    date_breaks = "24 months",
    date_labels = "%Y"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1
    )
  ) 

datos_trans <- datos_filtrado |>
  mutate(
    Original  = muertes,
    Logaritmo = log(muertes),
    BoxCox = box_cox(muertes, -1)
  )

datos_trans |>
  as_tibble() |>
  pivot_longer(c(Original, Logaritmo, BoxCox),
               names_to = "serie", values_to = "valor") |>
  mutate(serie = factor(serie, levels = c("Original", "Logaritmo", "BoxCox"))) |>
  ggplot(aes(x = Month, y = valor)) +
  geom_line(color = "cadetblue4") +
  facet_wrap(~ serie, ncol = 1, scales = "free_y") +
  labs(x = "Fecha", y = NULL) +
  theme_minimal()

datos_trans |>
  as_tibble() |>
  mutate(año = year(Month)) |>
  pivot_longer(c(Original, Logaritmo, BoxCox),
               names_to = "serie", values_to = "valor") |>
  mutate(serie = factor(serie, levels = c("Original", "Logaritmo", "BoxCox"))) |>
  group_by(serie, año) |>
  summarise(media = mean(valor), desvio = sd(valor), .groups = "drop") |>
  ggplot(aes(x = media, y = desvio)) +
  geom_point(color = "cadetblue4") +
  geom_smooth(method = "lm", se = FALSE, color = "gray40", linewidth = 0.5) +
  facet_wrap(~ serie, scales = "free") +
  scale_x_continuous(n.breaks = 3) +
  labs(x = "Media anual", y = "Desvío estándar anual") +
  theme_minimal()

# Seguimos con la serie sin transformar ya que los graficos de la serie dan iguales 
# ¿Cuántas diferencias estacionales?
datos_filtrado |> features(muertes, unitroot_nsdiffs)

# ¿Cuántas diferencias regulares, tras la estacional?
datos_filtrado |>
  mutate(dif12 = difference(muertes, 12)) |>
  features(dif12, unitroot_ndiffs)

grafico_dif <- function(d, ylab) {
  d <- d |> filter(!is.na(dif))   # saca los NA de la diferenciación (y el aviso)
  
  p_serie <- ggplot(d, aes(x = Month, y = dif)) +
    geom_line(color = "cadetblue4") +
    geom_point(color = "cadetblue4", size = 1) +
    labs(x = "Fecha", y = ylab) +
    theme_minimal()
  
  p_acf <- d |> ACF(dif, lag_max = 36) |> autoplot() +
    labs(x = "Rezago (meses)", y = "FAC") +
    theme_minimal()
  
  p_pacf <- d |> PACF(dif, lag_max = 36) |> autoplot() +
    labs(x = "Rezago (meses)", y = "FACP") +
    theme_minimal()
  
  p_serie / (p_acf | p_pacf)
}
#graficamos las serie con dos diferenciaciones 
datos_filtrado |>
  mutate(dif = difference(muertes, 1)) |>
  grafico_dif("Serie con diferencia regular")

datos_filtrado |>
  mutate(dif = difference(difference(muertes, 12),1)) |>
  grafico_dif("Serie con diferencia regula y estacional")

# $\text{SARIMA}(p,d,q)(P,D,Q)_{12}$: (1) $(2,0,0)(2,0,0)$, (2) $(0,1,1)(0,1,1)$

# Crear la serie temporal
# No solo tenerlo en "data frame" sino convertirlo en un objeto de tipo "serie temporal".

serie_acv <- ts(datos_filtrado$muertes,
                   frequency = 12,
                   start = c(1999, 1)) 

fit <- datos_filtrado |>
  model(
    modelo_1    = ARIMA(muertes ~ 1 + pdq(2,0,0) + PDQ(2,0,0)),
    modelo_2    = ARIMA(muertes ~ 0 + pdq(0,1,1) + PDQ(0,1,1)),
    modelo_3    = ARIMA(muertes ~ 0 + pdq(0,1,2) + PDQ(0,1,1)),
    modelo_auto = ARIMA(muertes ~ 0 + pdq(0,1,2) + PDQ(1,1,2))
  )

glance(fit) |>
  arrange(AICc) |>
  mutate(Modelo = recode(.model,
                         m1   = "ARIMA(2, 0, 0)(2, 0, 0)[12]",
                         m2   = "ARIMA(0,1,1)(0,1,1)[12]",
                         m3   = "ARIMA(0,1,2)(0,1,1)[12]",
                         auto = "Automático: ARIMA(0,1,2)(1,1,2)[12]")) |>
  select(Modelo, AICc, BIC) |>
  kable(digits = 2, align = "c")

inicio <- yearmonth("2000 Feb")   # primera observación con residuo válido

# Residuos de innovación estandarizados (escala Box-Cox)
res_std <- function(nombre) {
  augment(fit |> select(all_of(nombre))) |>
    filter(!is.na(.innov), Month >= inicio) |>
    mutate(est = .innov / sd(.innov))
}

# 1) Coeficientes con IC al 95%
tab_coef <- function(nombre) {
  fit |> select(all_of(nombre)) |> tidy() |>
    mutate(
      IC_inf = estimate - 1.96 * std.error,
      IC_sup = estimate + 1.96 * std.error,
      Significativo = ifelse(IC_inf > 0 | IC_sup < 0, "Sí", "No")
    ) |>
    select(Parámetro = term, Estimación = estimate, `Error estándar` = std.error,
           `IC 95% inferior` = IC_inf, `IC 95% superior` = IC_sup) |>
    kable(digits = 4, align = "c")
}

# 2) Test de Portmanteau (Ljung-Box)
tab_portmanteau <- function(nombre, lags = c(6, 12, 18, 24)) {
  r   <- augment(fit |> select(all_of(nombre))) |> filter(Month >= inicio)
  dof <- fit |> select(all_of(nombre)) |> tidy() |> nrow()
  map_df(lags, \(L) r |> features(.innov, ljung_box, lag = L, dof = dof) |>
           mutate(Rezago = L, gl = L - dof)) |>
    select(Rezago, Estadística = lb_stat, gl, `p-value` = lb_pvalue) |>
    kable(digits = c(0, 3, 0, 3), align = "c")
}

# 3) Residuos estandarizados + FAC + FACP
graf_residuos <- function(nombre) {
  r <- res_std(nombre)
  p1 <- ggplot(r, aes(Month, est)) +
    geom_line(color = "cadetblue4") +
    geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "gray40") +
    labs(x = "Fecha", y = "Residuos estandarizados") + theme_minimal()
  p2 <- r |> ACF(est, lag_max = 36) |> autoplot() +
    labs(x = "Rezago (meses)", y = "FAC") + theme_minimal()
  p3 <- r |> PACF(est, lag_max = 36) |> autoplot() +
    labs(x = "Rezago (meses)", y = "FACP") + theme_minimal()
  p1 / (p2 | p3)
}

# 4) Histograma + QQ-plot
graf_normalidad <- function(nombre) {
  r <- res_std(nombre)
  p1 <- ggplot(r, aes(est)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "cadetblue4", color = "white") +
    stat_function(fun = dnorm, args = list(mean = mean(r$est), sd = sd(r$est)),
                  color = "violetred", linewidth = 1) +
    labs(x = "Residuos estandarizados", y = "Densidad") + theme_minimal()
  p2 <- ggplot(r, aes(sample = est)) +
    stat_qq(color = "cadetblue4", size = 1.6) +
    stat_qq_line(color = "gray40") +
    labs(x = "Cuantiles teóricos", y = "Cuantiles de los residuos") + theme_minimal()
  p1 | p2
}

tab_coef("modelo_3")
tab_portmanteau("modelo_3")
graf_residuos("modelo_3")
graf_normalidad("modelo_3")
shapiro.test(res_std("modelo_3")$est)

tab_coef("modelo_auto")
tab_portmanteau("modelo_auto")
graf_residuos("modelo_auto")
graf_normalidad("modelo_auto")
shapiro.test(res_std("modelo_auto")$est)


tab_coef("modelo_2")
tab_portmanteau("modelo_2")
graf_residuos("modelo_2")
graf_normalidad("modelo_2")
shapiro.test(res_std("modelo_2")$est)


# AHORA PRUEBO LA SERIE CON TRANSFORMAICON LOGARITMICA 

fit <- datos_filtrado |>
  model(
    modelo_1  = ARIMA(log(muertes) ~ 1 + pdq(2,0,0) + PDQ(2,0,0)),
    modelo_2  = ARIMA(log(muertes) ~ 0 + pdq(0,1,1) + PDQ(0,1,1)),
    modelo_2b = ARIMA(log(muertes) ~ 0 + pdq(0,1,2) + PDQ(0,1,1))
  )

inicio_modelo <- list(
  modelo_1  = yearmonth("1999 Jan"),
  modelo_2  = yearmonth("2000 Feb"),
  modelo_2b = yearmonth("2000 Feb")
)

glance(fit) |>
  arrange(AICc) |>
  mutate(Modelo = recode(.model,
                         m1   = "ARIMA(2, 0, 0)(2, 0, 0)[12]",
                         m2   = "ARIMA(0,1,1)(0,1,1)[12]",
                         auto = "Automático: ARIMA(0,1,2)(0,1,1)[12]")) |>
  select(Modelo, AICc, BIC) |>
  kable(digits = 2, align = "c")

tab_coef("modelo_2b")
tab_portmanteau("modelo_2b")
graf_residuos("modelo_2b")
graf_normalidad("modelo_2b")
shapiro.test(res_std("modelo_2b")$est)


tab_coef("modelo_2")
tab_portmanteau("modelo_2")



# Los 12 datos reservados (2019)
datos_test <- datos |>
  filter_index("2019 Jan" ~ "2019 Dec")

# Pronóstico del modelo seleccionado (asegúrate de que m2 sea tu modelo SARIMA(0,1,2)(0,1,1))
fc <- fit |>
  select(modelo_3) |>
  forecast(h = nrow(datos_test))

# Generación de la tabla de pronósticos con el intervalo del 95%
tabla_pronosticos <- fc |>
  hilo(level = 95) |>
  unpack_hilo(`95%`) |>
  as_tibble() |>
  left_join(
    datos_test |> as_tibble() |> rename(real = muertes),
    by = "Month"
  ) |>
  transmute(
    Mes = format(as.Date(Month), "%B"),
    Valor_real = real,
    Pronostico = round(.mean),
    IC_95 = paste0("(", round(`95%_lower`), "; ", round(`95%_upper`), ")")
  )

# Visualización en tabla con kable
kable(tabla_pronosticos,
      col.names = c("Mes", "Valor real", "Pronóstico", "IC del 95%"),
      align = "c",
      format.args = list(
        big.mark = ".",
        decimal.mark = ","
      ))



