# packages
library(terra)         # for SpatRaster and c()
library(virtualspecies)
library(imageRy)
library(geodata)
library(ggplot2)

# climate (WorldClim v2.1 at 0.5 arc-min)
tmin <- worldclim_country("Lux", "tmin", path = tempdir(), res = 0.5, version = "2.1")
tmax <- worldclim_country("Lux", "tmax", path = tempdir(), res = 0.5, version = "2.1")
prec <- worldclim_country("Lux", "prec", path = tempdir(), res = 0.5, version = "2.1")

# July layers
tmin_jul <- tmin[[7]]
tmax_jul <- tmax[[7]]
prec_jul <- prec[[7]]

# environmental stack
env_stack_base <- c(tmin_jul, tmax_jul, prec_jul)
names(env_stack_base) <- c("tmin", "tmax", "prec")
summary(env_stack_base)

# 1) Species response curves consistent with the area's summary (tmin~12, tmax~22, prec~75):
#    - dec: prefers cooler & drier -> worsens with warming/wetting
#    - inc: prefers warmer & wetter -> improves with warming/wetting
#    - stab: broad tolerance around current values -> roughly stable
curve_list <- list(
  dec  = formatFunctions( # decreasing: slightly cooler/drier than 2000
    tmin = c(fun="dnorm", mean=11.5, sd=1.2),
    tmax = c(fun="dnorm", mean=21.0, sd=1.5),
    prec = c(fun="dnorm", mean=70,   sd=6)
  ),
  inc  = formatFunctions( # increasing: slightly warmer/wetter than 2000
    tmin = c(fun="dnorm", mean=13.5, sd=1.2),
    tmax = c(fun="dnorm", mean=24.0, sd=1.5),
    prec = c(fun="dnorm", mean=82,   sd=6)
  ),
  stab = formatFunctions( # stable: centered and tolerant
    tmin = c(fun="dnorm", mean=12.0, sd=3.5),
    tmax = c(fun="dnorm", mean=22.0, sd=3.5),
    prec = c(fun="dnorm", mean=75,   sd=12)
  )
)

# 2) Parametric climate change: default ~ +0.2 °C/year and +0.5 precip units/year
get_env_for_year <- function(year, base_stack,
                             dt_tmin = 0.2,  # °C/year
                             dt_tmax = 0.2,  # °C/year
                             dt_prec = 0.5   # precip units/year (same units as 'prec')
) {
  d <- year - 2000
  r <- if (inherits(base_stack, "SpatRaster")) base_stack else terra::rast(base_stack)
  r[["tmin"]] <- r[["tmin"]] + dt_tmin * d
  r[["tmax"]] <- r[["tmax"]] + dt_tmax * d
  r[["prec"]] <- r[["prec"]] + dt_prec * d
  r
}

# 3) Simulation: returns a stack per species (layers = years), values in [0,1]
simulate_three_species <- function(env_stack_base,
                                   years = 2000:2020,
                                   curves = curve_list,
                                   seed = 8888,
                                   dt_tmin = 0.2, dt_tmax = 0.2, dt_prec = 0.5,
                                   rescale = TRUE) {
  set.seed(seed)
  base_spat <- if (inherits(env_stack_base, "SpatRaster")) env_stack_base else terra::rast(env_stack_base)
  
  get_suit <- function(vs_obj) if (is.list(vs_obj) && "suitab.raster" %in% names(vs_obj)) vs_obj$suitab.raster else vs_obj
  
  out <- vector("list", length(curves))
  names(out) <- names(curves)
  
  for (nm in names(curves)) {
    stk <- NULL
    for (yr in years) {
      env_y <- get_env_for_year(yr, base_spat, dt_tmin, dt_tmax, dt_prec)
      vs    <- virtualspecies::generateSpFromFun(env_y, parameters = curves[[nm]], rescale = rescale)
      s     <- get_suit(vs)
      s     <- if (inherits(s, "SpatRaster")) s else terra::rast(s)
      names(s) <- paste0("y", yr)
      stk <- if (is.null(stk)) s else c(stk, s)
    }
    out[[nm]] <- stk
  }
  out
}

species_suitability <- simulate_three_species(
  env_stack_base,
  years   = 2000:2005,                 # or 2000:2020
  curves  = curve_list,
  dt_tmin = 0.30, dt_tmax = 0.30,      # ~ +1.5 °C in 5 years
  dt_prec = 1.0,                       # ~ +5 in 5 years
  rescale = TRUE                       # important, avoids collapse to ~0
)

plot(species_suitability$dec)   # SpatRaster: layers y2000..y2020 (declining species)

im.ridgeline(species_suitability$inc,  scale = 2, palette = "viridis")
im.ridgeline(species_suitability$dec,  scale = 2, palette = "viridis")
im.ridgeline(species_suitability$stab, scale = 2, palette = "magma")

im.ridgeline(species_suitability$stab, scale = 2, palette = "viridis") +
  labs(
    y    = "year",
    x    = "values",
    fill = "suitability"
  )

# ---- Map plotting ----
# rename layers removing the "y" prefix (FIX: use 'stab' names rather than 'dec')
names(species_suitability$stab) <- gsub("y", "", names(species_suitability$stab))

plot(
  species_suitability$stab,
  col  = hcl.colors(20, "viridis"),
  axes = FALSE,
  mar  = c(3, 3, 2, 6),
  plg  = list(title = "Suitability")
)

# =========================
# Saver: write stacks to disk
# =========================
# Save the whole object (quickest way to reload everything later)
saveRDS(species_suitability, "species_suitability.Rds")

# Save ONE multi-band GeoTIFF per species (bands = years)
out_dir <- "ss_out"
dir.create(out_dir, showWarnings = FALSE)

for (sp in names(species_suitability)) {
  terra::writeRaster(
    species_suitability[[sp]],
    filename = file.path(out_dir, paste0(sp, "_stack.tif")),
    overwrite = TRUE
  )
}
