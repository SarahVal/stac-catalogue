

#' @name create_projection
#' @param lon string, name of the longitude column
#' @param lat string, name of the latitude column
#' @param proj_from character, initial projection of the xy coordinates
#' @param proj_to character, target projection
#' @param new_lon character, name of the new longitude column
#' @param new_lat character, name of the new latitude column
#' @return a dataframe with two columns in the proj_to projection
#' @import dplyr
#' @export
create_projection <- function(obs,
                              lon,
                              lat,
                              proj_from,
                              proj_to,
                              new_lon = NULL,
                              new_lat = NULL) {
  if (is.null(new_lon)) {
    new_lon <- lon
  }
  
  if (is.null(new_lat)) {
    new_lat <- lat
  }
  
  new.coords <- project_coords(obs, lon, lat, proj_from, proj_to)
  new.coords.df <- data.frame(new.coords) %>%
    setNames(c(new_lon, new_lat))
  
  suppressWarnings(obs <- obs %>%
                     dplyr::select(-one_of(c(new_lon, new_lat))) %>% dplyr::bind_cols(new.coords.df))
  
  return(obs)
}

#' @name project_coords
#' @param xy data frame, containing the coordinates to reproject
#' @param lon string, name of the longitude column
#' @param lat string, name of the latitude column
#' @param proj_from character, initial projection of the xy coordinates
#' @param proj_to character, target projection
#' @import sp dplyr
#' @return spatial points in the proj_to projection
#' @export
project_coords <-
  function(xy,
           lon = "lon",
           lat = "lat",
           proj_from,
           proj_to = NULL) {
    xy <- dplyr::select(xy, dplyr::all_of(c(lon, lat)))
    sp::coordinates(xy) <- c(lon, lat)
    sp::proj4string(xy) <- sp::CRS(proj_from)
    
    if (!is.null(proj_to)) {
      xy <- sp::spTransform(xy, sp::CRS(proj_to))
    }
    xy
  }

#' @name points_to_bbox
#' @param xy data frame, containing the coordinates to reproject
#' @param buffer integer, buffer to add around the observations
#' @param proj_from character, initial projection of the xy coordinates
#' @param proj_to character, target projection
#' @return a box extent
#' @export
points_to_bbox <-
  function(xy,
           buffer = 0,
           proj_from = NULL,
           proj_to = NULL) {
    if (!inherits(xy, "SpatialPoints")) {
      sp::coordinates(xy) <- colnames(xy)
      proj4string(xy) <- sp::CRS(proj_from)
    }
    
    bbox <-
      sf::st_buffer(sf::st_as_sfc(sf::st_bbox(xy)), dist = buffer)
    
    if (!is.null(proj_to)) {
      bbox <- bbox %>%
        sf::st_transform(crs = sp::CRS(proj_to))
      proj_from <- proj_to
    }
    
    bbox %>% sf::st_bbox(crs = proj_from)
  }

#' @export
shp_to_bbox <- function(shp,
                        proj_from = NULL,
                        proj_to = NULL) {
  if (is.na(sf::st_crs(shp)) && is.null(proj_from)) {
    stop("proj.fom is null and shapefile has no crs.")
  }
  
  if (is.na(sf::st_crs(shp))) {
    crs(shp) <- proj_from
    shp <- shp %>% sf::st_set_crs(proj_from)
  }
  
  if (!is.null(proj_to)) {
    shp <- shp %>%
      sf::st_transform(crs = sp::CRS(proj_to))
  }
  
  
  bbox <- sf::st_bbox(shp, crs = proj)
  
  bbox
}


#' Create a proxy data cube for current climate,
#' which loads data from a given image collection according to a data cube view based
#' on a specific box coordinates or using a set of observations
#'
#' @name load_cube
#'
#' @param stac_path, a character, base url of a STAC web service.
#' @param limit, an integer defining the maximum number of results to return.
#' @param collections, a character vector of collection IDs to include
#' subsetLayers, a vector, containing the name of layers to select. If NULL, all layers in dir.pred selected by default.
#' @param use.obs, a boolean. If TRUE, the provided observations will be sued as a basis for calculating the extent and bbox.
#' @param obs, a data.frame containg the observations (used if use.obs is T)
#' @param srs.obs, string, observations spatial reference system. Can be a proj4 definition, WKT, or in the form "EPSG:XXXX".
#' @param lon, a string, column from obs containing longitude
#' @param lat, a string, column from obs containing latitude
#' @param buffer.box, an integer, buffer to apply around the obs to calculate extent and bbox
#' @param bbox, a numeric vector of size 4 or 6. Coordinates of the bounding box (if use.obs is FALSE). Details in rstac::stac_search documentation.
#' @param layers, a string vector, names of bands to be used,. By default (NULL), all bands with "eo:bands" attributes will be used.
#' @param srs.cube, string, target spatial reference system. Can be a proj4 definition, WKT, or in the form "EPSG:XXXX".
#' @param t0, ISO8601 datetime string, start date.
#' @param t1, ISO8601 datetime string, end date.
#' @param left, a float. Left coordinate of the extent. Used if use.obs = F
#' @param right, a float. Right coordinate of the extent. Used if use.obs = F
#' @param top, a float. Top coordinate of the extent. Used if use.obs = F
#' @param bottom, a float. Bottom coordinate of the extent. Used if use.obs = F
#' @param spatial.res, a float, size of pixels in longitude and latitude directions, in the unit of srs.cube spatial reference system.
#' @param temporal.res, size of pixels in time-direction, expressed as ISO8601 period string (only 1 number and unit is allowed) such as "P16D"
#' @param aggregation, a character, aggregation method as string, defining how to deal with pixels containing data from multiple images, can be "min", "max", "mean", "median", or "first"
#' @param resampling, a character, resampling method used in gdalwarp when images are read, can be "near", "bilinear", "bicubic" or others as supported by gdalwarp (see https://gdal.org/programs/gdalwarp.html)
#' @return a raster stack of variables not intercorrelated
#' @import gdalcubes dplyr sp sf rstac
#' @return a proxy raster data cube
#' @export
load_cube <-
  function(stac_path = "https://io.biodiversite-quebec.ca/stac/",
           limit = NULL,
           collections = c("chelsa-clim"),
           bbox = NULL,
           layers = NULL,
           variable = NULL,
           mask = NULL,
           srs.cube = "EPSG:6623",
           t0 = NULL,
           t1 = NULL,
           spatial.res = NULL,
           temporal.res = "P1Y",
           aggregation = "mean",
           resampling = "near") {
    
    s <- rstac::stac(stac_path)
   
    if(!inherits(bbox, "bbox")) stop("The bbox is not a bbox object.")
    
      left <- bbox$xmin
      right <- bbox$xmax
      bottom <- bbox$ymin
      top <- bbox$ymax
      
      if (left > right) {
        stop("left and right seem reversed")
      }
      if (bottom > top) {
        stop("bottom and top seem reversed")
      }
      
      # Create the bbxo (WGS84 projection)
      bbox.wgs84 <- bbox %>% sf::st_bbox(crs = srs.cube) %>% 
        sf::st_as_sfc() %>% sf::st_transform(crs = "EPSG:4326") %>% 
        sf::st_bbox()
   
    
    if (!is.null(t0)) {
      datetime <- format(lubridate::as_datetime(t0), "%Y-%m-%dT%H:%M:%SZ")
    } else {
      it_obj_tmp <- s %>%
        rstac::stac_search(bbox = bbox.wgs84,
                           collections = collections,
                           limit = limit) %>%
        rstac::get_request()
      datetime <- it_obj_tmp$features[[1]]$properties$datetime
      t0 <- datetime
      t1 <- datetime
    }
    if (!is.null(t1) && t1 != t0) {
      datetime <- paste(datetime,
                        format(lubridate::as_datetime(t1),
                               "%Y-%m-%dT%H:%M:%SZ"),
                        sep = "/")
    }
    RCurl::url.exists(stac_path)
    it_obj <- s %>%
      rstac::stac_search(bbox = bbox.wgs84,
                         collections = collections,
                         datetime = datetime,
                         limit = limit) %>%
      rstac::get_request()
    if (is.null(spatial.res)) {
      name1 <- unlist(lapply(it_obj$features, function(x) {
        names(x$assets)
      }))[1]
      spatial.res <-
        it_obj$features[[1]]$assets[[name1]]$`raster:bands`[[1]]$spatial_resolution
    }
    if (is.null(layers)) {
      layers <- unlist(lapply(it_obj$features, function(x) {
        names(x$assets)
      }))
    }
    if (!is.null(variable)) {
      st <- gdalcubes::stac_image_collection(
        it_obj$features,
        asset_names = layers,
        property_filter = function(x) {
          x[["variable"]] %in% variable
        }
      )
    } else {
      st <- gdalcubes::stac_image_collection(it_obj$features,
                                             asset_names = layers)
    }
    v <- gdalcubes::cube_view(
      srs = srs.cube,
      extent = list(
        t0 = t0,
        t1 = t1,
        left = left,
        right = right,
        top = top,
        bottom = bottom
      ),
      dx = spatial.res,
      dy = spatial.res,
      dt = temporal.res,
      aggregation = aggregation,
      resampling = resampling
    )
    gdalcubes::gdalcubes_options(parallel = T)
    
    cube <- gdalcubes::raster_cube(st, v, mask)
    return(cube)
  }


#' Create a proxy data cube for future climate,
#' which loads data from a given image collection according to a data cube view based
#' on a specific box coordinates or using a set of observations
#'
#' @name load_cube_projection
#'
#' @param stac_path, a character, base url of a STAC web service.
#' @param limit, an integer defining the maximum number of results to return.
#' @param collections, a character vector of collection IDs to include
#' subsetLayers, a vector, containing the name of layers to select. If NULL, all layers in dir.pred selected by default.
#' @param use.obs, a boolean. If TRUE, the provided observations will be sued as a basis for calculating the extent and bbox.
#' @param obs, a data.frame containg the observations (used if use.obs is T)
#' @param srs.obs, string, observations spatial reference system. Can be a proj4 definition, WKT, or in the form "EPSG:XXXX".
#' @param lon, a string, column from obs containing longitude
#' @param lat, a string, column from obs containing latitude
#' @param buffer.box, an integer, buffer to apply around the obs to calculate extent and bbox
#' @param bbox, a numeric vector of size 4 or 6. Coordinates of the bounding box (if use.obs is FALSE). Details in rstac::stac_search documentation.
#' @param layers, a string vector, names of bands to be used,. By default (NULL), all bands with "eo:bands" attributes will be used.
#' @param srs.cube, string, target spatial reference system. Can be a proj4 definition, WKT, or in the form "EPSG:XXXX".
#' @param time.span, a string, time interval of the projection model.
#' @param rcp, a string, climatic scenario
#' @param left, a float. Left coordinate of the extent. Used if use.obs = F
#' @param right, a float. Right coordinate of the extent. Used if use.obs = F
#' @param top, a float. Top coordinate of the extent. Used if use.obs = F
#' @param bottom, a float. Bottom coordinate of the extent. Used if use.obs = F
#' @param spatial.res, a float, size of pixels in longitude and latitude directions, in the unit of srs.cube spatial reference system.
#' @param temporal.res, size of pixels in time-direction, expressed as ISO8601 period string (only 1 number and unit is allowed) such as "P16D"
#' @param aggregation, a character, aggregation method as string, defining how to deal with pixels containing data from multiple images, can be "min", "max", "mean", "median", or "first"
#' @param resampling, a character, resampling method used in gdalwarp when images are read, can be "near", "bilinear", "bicubic" or others as supported by gdalwarp (see https://gdal.org/programs/gdalwarp.html)
#' @return a raster stack of variables not intercorrelated
#' @import gdalcubes dplyr sp sf rstac
#' @return a proxy raster data cube
#' @export
load_cube_projection <- function(stac_path =
                                   "https://io.biodiversite-quebec.ca/stac/",
                                 limit = NULL,
                                 collections = c("chelsa-clim-proj"),
                                 bbox = NULL,
                                 layers = NULL,
                                 variable = NULL,
                                 srs.cube = "EPSG:6623",
                                 time.span = "2041-2070",
                                 rcp = "ssp585",
                                 left = -2009488,
                                 right = 1401061,
                                 bottom = -715776,
                                 top = 2597757,
                                 spatial.res = 2000,
                                 temporal.res = "P1Y",
                                 aggregation = "mean",
                                 resampling = "near") {
  # t0 param
  if (time.span == "2011-2040") {
    t0 <- "2011-01-01"
  }
  
  if (time.span == "2041-2070") {
    t0 <- "2041-01-01"
  }
  
  if (time.span == "2071-2100") {
    t0 <- "2071-01-01"
  }
  
  datetime <-
    format(lubridate::as_datetime(t0), "%Y-%m-%dT%H:%M:%SZ")
  s <- rstac::stac(stac_path)
  
    left <- bbox$xmin
    right <- bbox$xmax
    bottom <- bbox$ymin
    top <- bbox$ymax
    
    if (left > right) {
      stop("left and right seem reversed")
    }
    if (bottom > top) {
      stop("bottom and top seem reversed")
    }
    
    # Create the bbxo (WGS84 projection)

    bbox.wgs84 <- bbox %>% sf::st_bbox(crs = srs.cube) %>% 
      sf::st_as_sfc() %>% sf::st_transform(crs = "EPSG:4326") %>% 
      sf::st_bbox()
 
  
  it_obj <- s %>%
    rstac::stac_search(bbox = bbox.wgs84,
                       collections = collections,
                       datetime = datetime,
                       limit = limit) %>%
    rstac::get_request() # bbox in decimal lon/lat
  
  # If no layers is selected, get all the layers by default
  if (is.null(layers)) {
    layers <- unlist(lapply(it_obj$features, function(x) {
      names(x$assets)
    }))
  }
  
  #
  # Creates an image collection
  if (!is.null(variable)) {
    st <- gdalcubes::stac_image_collection(
      it_obj$features,
      asset_names = layers,
      property_filter = function(x) {
        x[["variable"]] %in% variable & x[["rcp"]] == rcp
      }
    )
  } else {
    st <- gdalcubes::stac_image_collection(
      it_obj$features,
      asset_names = layers,
      property_filter = function(x) {
        x[["rcp"]] == rcp
      }
    )
  }
  
  
  # if layers = NULL, load all the layers
  v <- gdalcubes::cube_view(
    srs = srs.cube,
    extent = list(
      t0 = t0,
      t1 = t0,
      left = left,
      right = right,
      top = top,
      bottom = bottom
    ),
    dx = spatial.res,
    dy = spatial.res,
    dt = temporal.res,
    aggregation = aggregation,
    resampling = resampling
  )
  gdalcubes::gdalcubes_options(threads = 4)
  cube <- gdalcubes::raster_cube(st, v)
  return(cube)
}



#' @export
cube_to_raster <- function(cube, format = "raster") {
  # Transform to a star object
  cube.xy <- cube %>%
    stars::st_as_stars()
  
  # If not, names are concatenated with temp file names
  
  # We remove the temporal dimension
  cube.xy <- cube.xy %>% abind::adrop(c(F, F, T))
  
  # Conversion to a spatial object
  
  if (format == "raster") {
    # Raster format
    cube.xy <- raster::stack(as(cube.xy, "Spatial"))
  } else {
    # Terra format
    cube.xy <- terra::rast(cube.xy)
  }
  names(cube.xy) <- names(cube)
  
  cube.xy
}

#' @export
extract_gdal_cube <-
  function(cube,
           n_sample = 5000,
           simplify = T,
           xy = F) {
    x_min <- gdalcubes::dimensions(cube)$x$low
    x_max <- gdalcubes::dimensions(cube)$x$high
    y_min <- gdalcubes::dimensions(cube)$y$low
    y_max <- gdalcubes::dimensions(cube)$y$high
    
    # ensure that the number of samplis is smaller than the number of pixels
    n_sample <-
      min(
        n_sample,
        gdalcubes::dimensions(cube)$x$count * gdalcubes::dimensions(cube)$y$count
      )
    
    x <- runif(n_sample, x_min, x_max)
    y <- runif(n_sample, y_min, y_max)
    
    df <-
      sf::st_as_sf(data.frame(x = x, y = y),
                   coords = c("x", "y"),
                   crs = srs(cube))
    value_points <- gdalcubes::extract_geom(cube, df)
    
    if (xy) {
      value_points <- bind_cols(value_points, data.frame(x = x, y = y))
    }
    if (simplify) {
      value_points <- value_points %>% dplyr::select(-FID,-time)
    }
    value_points
  }



#' @export
get_info_collection <- function(stac_path =
                                  "https://io.biodiversite-quebec.ca/stac/",
                                limit = 5000,
                                collections = c("chelsa-clim"),
                                bbox = NULL) {
  # Creating RSTACQuery  query
  s <- rstac::stac(stac_path)
  
  if (is.null(bbox)) {
    bbox <- c(
      xmin = -180,
      xmax = 180,
      ymax = 90,
      ymin = -90
    )  %>% 
      sf::st_bbox(crs = "EPSG:4326")
  }
  
  
  # CreateRSTACQuery object with the subclass search containing all search field parameters
  it_obj <- s %>% # think changing it for %>%
    rstac::stac_search(bbox = bbox,
                       collections = collections,
                       limit = limit) %>% rstac::get_request()
  
  layers <- unlist(lapply(it_obj$features, function(x) {
    names(x$assets)
  }))
  temporal_extent <- unlist(lapply(it_obj$features, function(x) {
    x$properties$datetime
  }))
  variable <- unique(unlist(lapply(it_obj$features, function(x) {
    x$properties$variable
  })))
  t0 <- min(temporal_extent)
  t1 <- max(temporal_extent)
  spatial_res <- unique(unlist(lapply(it_obj$features,
                                      function(x) {
                                        lapply(x$assets,
                                               function(x) {
                                                 lapply(x$`raster:bands`,
                                                        function(x) {
                                                          x$spatial_resolution
                                                        })
                                               })
                                      }), use.names = F))
  
  return(
    list(
      "layers" = layers,
      "variables" = variable,
      "t0" = t0,
      "t1" = t1,
      "spatial_resolution" = spatial_res
    )
  )
}

#' Create a proxy data cube for current climate,
#' which loads data from a given image collection according to a data cube view based
#' on a specific box coordinates or using a set of observations
#'
#' @name load_prop_values
#'
#' @param stac_path, a character, base url of a STAC web service.
#' @param limit, an integer defining the maximum number of results to return.
#' @param collections, a character vector of collection IDs to include
#' subsetLayers, a vector, containing the name of layers to select. If NULL, all layers in dir.pred selected by default.
#' @param use.obs, a boolean. If TRUE, the provided observations will be sued as a basis for calculating the extent and bbox.
#' @param obs, a data.frame containg the observations (used if use.obs is T)
#' @param srs.obs, string, observations spatial reference system. Can be a proj4 definition, WKT, or in the form "EPSG:XXXX".
#' @param lon, a string, column from obs containing longitude
#' @param lat, a string, column from obs containing latitude
#' @param buffer.box, an integer, buffer to apply around the obs to calculate extent and bbox
#' @param bbox, a numeric vector of size 4 or 6. Coordinates of the bounding box (if use.obs is FALSE). Details in rstac::stac_search documentation.
#' @param layers, a string vector, names of bands to be used,. By default (NULL), all bands with "eo:bands" attributes will be used.
#' @param srs.cube, string, target spatial reference system. Can be a proj4 definition, WKT, or in the form "EPSG:XXXX".
#' @param t0, ISO8601 datetime string, start date.
#' @param t1, ISO8601 datetime string, end date.
#' @param left, a float. Left coordinate of the extent. Used if use.obs = F
#' @param right, a float. Right coordinate of the extent. Used if use.obs = F
#' @param top, a float. Top coordinate of the extent. Used if use.obs = F
#' @param bottom, a float. Bottom coordinate of the extent. Used if use.obs = F
#' @param spatial.res, a float, size of pixels in longitude and latitude directions, in the unit of srs.cube spatial reference system.
#' @param temporal.res, size of pixels in time-direction, expressed as ISO8601 period string (only 1 number and unit is allowed) such as "P16D"
#' @param prop, a boolean. If TRUE, the proportion of each classes from selected_values is calculated. If FALSE, only the values in selected_values are retrieved.
#' @param prop.res, an integer, resolution to calculate the proportion of land cover classes (if prop is TRUE)
#' @param selected_values, a vector, classes values to select.
#' @return a raster stack of variables not intercorrelated
#' @import gdalcubes dplyr sp sf rstac
#' @return a proxy raster data cube
#' @export
load_prop_values <-
  function(stac_path = "https://io.biodiversite-quebec.ca/stac/",
           collections = c("esacci-lc"),
           srs.cube = "EPSG:6623",
           t0 = "2000-01-01",
           t1 = "2001-12-31",
           spatial.res = 250,
           bbox = NULL,
           limit = 5000,
           prop = F,
           prop.res = 1000,
           select_values = NULL,
           layers = NULL,
           variable = NULL,
           temporal.res = "P1Y") {
    s <- rstac::stac(stac_path)

    left <- bbox$xmin
    right <- bbox$xmax
    bottom <- bbox$ymin
    top <- bbox$ymax
    
    
    if (left > right) {
      stop("left and right seem reversed")
    }
    if (bottom > top) {
      stop("bottom and top seem reversed")
    }
    
    # Create the bbox in WGS84 projection)
    bbox.wgs84 <- bbox %>% sf::st_bbox(crs = srs.cube) %>% 
      sf::st_as_sfc() %>% sf::st_transform(crs = "EPSG:4326") %>% 
      sf::st_bbox()


    if (!is.null(t0)) {
      datetime <- format(lubridate::as_datetime(t0), "%Y-%m-%dT%H:%M:%SZ")
    } else {
      it_obj_tmp <- s %>%
        rstac::stac_search(bbox = bbox.wgs84,
                           collections = collections,
                           limit = limit) %>%
        rstac::get_request()
      datetime <- it_obj_tmp$features[[1]]$properties$datetime
      t0 <- datetime
      t1 <- datetime
    }
    if (!is.null(t1) && t1 != t0) {
      datetime <- paste(datetime,
                        format(lubridate::as_datetime(t1),
                               "%Y-%m-%dT%H:%M:%SZ"),
                        sep = "/")
    }
    RCurl::url.exists(stac_path)
    it_obj <- s %>%
      rstac::stac_search(bbox = bbox.wgs84,
                         collections = collections,
                         datetime = datetime,
                         limit = limit) %>%
      rstac::get_request()
    if (is.null(spatial.res)) {
      name1 <- unlist(lapply(it_obj$features, function(x) {
        names(x$assets)
      }))[1]
      spatial.res <-
        it_obj$features[[1]]$assets[[name1]]$`raster:bands`[[1]]$spatial_resolution
    }
    if (is.null(layers)) {
      layers <- unlist(lapply(it_obj$features, function(x) {
        names(x$assets)
      }))
    }
    
    v <- gdalcubes::cube_view(
      srs = srs.cube,
      extent = list(
        t0 = t0,
        t1 = t1,
        left = left,
        right = right,
        top = top,
        bottom = bottom
      ),
      dx = spatial.res,
      dy = spatial.res,
      dt = temporal.res,
      aggregation = "mode",
      resampling = "near"
    )
    gdalcubes::gdalcubes_options(parallel = T)
    
    cube_class_rstack <- raster::stack()
    
    for (i in 1:length(unique(select_values))) {
      for (j in 1:length(layers)) {
        if (!is.null(variable)) {
          st <- gdalcubes::stac_image_collection(
            it_obj$features,
            asset_names = layers[j],
            property_filter = function(x) {
              x[["variable"]] %in% variable
            }
          )
        } else {
          st <- gdalcubes::stac_image_collection(it_obj$features,
                                                 asset_names = layers[j])
        }
        
        
        cube_class <-
          raster_cube(st,
                      v,
                      mask = image_mask(
                        band = layers[j],
                        values = select_values[i],
                        invert = T
                      ))
        
        if (prop) {
          cube_class <- cube_class %>%
            aggregate_space(dx = prop.res,
                            dy = prop.res,
                            method = "count") %>%
            stars::st_as_stars() %>%
            as("Raster")
          cube_class <- cube_class / ((prop.res / spatial.res) ^ 2)
        } else {
          cube_class <- cube_class %>%
            stars::st_as_stars() %>%
            as("Raster")
        }
        
        cube_class <- cube_class[[j]]
        cube_class
        names(cube_class) <-
          paste0("y", substring(layers[j], 11, 15), "_class", select_values[i])
        cube_class_rstack <-
          raster::stack(cube_class_rstack, cube_class)
      }
    }
    
    return(cube_class_rstack)
  }
