#require(dplyr)
require(raster)



## --
# Fonction pour récupérer une valeur maximale locale
Get_Local_Maximum <- function(m_grid, row, col, nodata)
{
  dcol <- c(0, 1, 1, 1, 0,-1,-1,-1)
  drow <- c(1, 1, 0,-1,-1,-1, 0, 1)
  ext <- dim(m_grid)
  z <- m_grid[row,col]

  for (ii in 1:8)
  {
    coln <- col + dcol[ii]
    rown <- row + drow[ii]
    
    if ( (rown <= ext[1] & rown >= 1) & (coln <= ext[2] & coln >= 1) )
    {
      val <- m_grid[rown, coln]
      if (val != nodata & val > z) z <- val
    }
  }
  return(z)
}


## --
# Calcul de pente
# Selon Zevenbergen & Thorne (1987)
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/esp.3290120107
# Tiré de CSG_Grid::Get_Gradient dans grid_operation.cpp (L1011)
Get_Gradient <- function(m_grid, row, col, nodata, cellsize = 1)
{
  dcol <- c(0, 1, 0,-1)
  drow <- c(1, 0,-1, 0)
  ext <- dim(m_grid)
  z <- m_grid[row, col]
  dz <- rep(0, 4)

  for (ii in 1:4)
  {
    colTo <- col + dcol[ii]
    rowTo <- row + drow[ii]
    colFrom <- col - dcol[ii]
    rowFrom <- row - drow[ii]
    
    # En fait toutes les conversion de is_InGrid devraient être comme ça
    # ce sera moins compliqué dans WhiteboxTools car on peut extraire une valeur
    # en dehors de l'extent et ça nous retournera la valeur de nodata plutôt qu'une
    # erreur comme dans R avec une matrice
    is_InGrid_To <- ifelse(
      (rowTo <= ext[1] & rowTo >= 1) & (colTo <= ext[2] & colTo >= 1),
      ifelse(m_grid[rowTo, colTo] != nodata, TRUE, FALSE),
      FALSE
      )
    is_InGrid_From <- ifelse(
      (rowFrom <= ext[1] & rowFrom >= 1) & (colFrom <= ext[2] & colFrom >= 1),
      ifelse(m_grid[rowFrom, colFrom] != nodata, TRUE, FALSE),
      FALSE
    )
    
    if (is_InGrid_To)
    {
      dz[ii] <- m_grid[rowTo, colTo] - z
    } else if (is_InGrid_From)
    {
      dz[ii] <- z - m_grid[rowFrom, colFrom]
    }
  }

  G	<- (dz[1] - dz[3]) / (2 * cellsize)
  H	<- (dz[2] - dz[4]) / (2 * cellsize)

  Slope	<- atan(sqrt(G*G + H*H))
  Aspect <- atan(-H/-G)
  return(c("Slope" = Slope, "Aspect" = Aspect))
}


## --
# Calcul des superficies (accumulation)
Get_Area <- function(m_pDEM, m_pWeight, suction, slope_weight, nodata, cellsize = 1)
{
  MFD_Converge <- 1.1
  ext <- dim(m_pDEM)
  nb_cells <- ext[1]*ext[2]

  # Initialisation des rasters (sous forme de matrices)
  m_Suction <- m_pDEM
  m_pArea <- m_pDEM
  m_pSlope <- m_pDEM

  m_Suction[m_pDEM != nodata] <- 0
  m_pArea[m_pDEM != nodata] <- 0
  m_pSlope[m_pDEM != nodata] <- 0

  # Variables pour analyse de voisinage
  dcol <- c(0, 1, 1, 1, 0,-1,-1,-1)
  drow <- c(1, 1, 0,-1,-1,-1, 0, 1)
  length_dcolrow <- sqrt((dcol*cellsize)^2 + (drow*cellsize)^2)

  # Création de l'index permettant de classer en ordre décroissant des valeurs d'élèvation
  tbl_index <- Create_Index(m_pDEM, nodata)
  
  
  for (n in 1:nrow(tbl_index))
  {
    z <- tbl_index[n, "value"]
    row <- tbl_index[n, "row"]
    col <- tbl_index[n, "col"]
    
    # Calcul de la pente initiale et de la suction
    Slope <- Get_Gradient(m_pDEM, row, col, nodata)[["Slope"]]
    t_param <- suction^(slope_weight * Slope)
    m_Suction[row,col] <- (1/t_param)^exp(t_param)

    # Ajustement de l'accumulation
    Area <- m_pArea[row,col] + m_pWeight[row,col]
    m_pArea[row,col] <- Area
    
    # Ajustement de la pente du catchment en fonction de l'accumulation
    Slope	<- m_pSlope[row,col] + Slope
    m_pSlope[row,col] <- Slope / Area


    dz <- rep(0, 8)

    for (ii in 1:8)
    {
      coln <- col + dcol[ii]
      rown <- row + drow[ii]
      if ( (rown <= ext[1] & rown >= 1) & (coln <= ext[2] & coln >= 1) )
      {
        val <- m_pDEM[rown, coln]
        if (val != nodata)
        {
          d <- z - val
          if (d > 0) { dz[ii] <- atan(d / length_dcolrow[ii]) ^ MFD_Converge }
        }
      }
    }
    dzSum <- sum(dz)


    if (dzSum > 0)
    {
      for (ii in which(dz > 0))
      {
        coln <- col + dcol[ii]
        rown <- row + drow[ii]
        if ( (rown <= ext[1] & rown >= 1) & (coln <= ext[2] & coln >= 1) )
        {
          m_pArea[rown,coln] <- m_pArea[rown,coln] + Area * dz[ii] / dzSum
          m_pSlope[rown,coln] <- m_pSlope[rown,coln] + Slope * dz[ii] / dzSum
        }
      }
    }
  }
  
  bool_nodata <- m_pArea != nodata
  m_pArea[bool_nodata] <- m_pArea[bool_nodata] * cellsize^2
  return(list(m_pArea=m_pArea, m_pSlope=m_pSlope, m_Suction=m_Suction))
}


## --
# Établissement de l'ordre de traitement des cellules
Create_Index <- function(m_grid, nodata)
{
  ext <- dim(m_grid)
  tbl_index <- dplyr::tibble(value = as.vector(m_grid))
  tbl_index["row"] <- rep(1:ext[1], ext[2])
  tbl_index["col"] <- rep(1:ext[2], each = ext[1])
  tbl_index <- tbl_index |>
    dplyr::filter(value != nodata) |>
    dplyr::arrange(dplyr::desc(value)) |>
    as.data.frame()
  
  return(tbl_index)
}


## --
# Modification de l'accumulation
Get_Modified <- function(m_pArea, m_pSlope, m_Suction, nodata)
{
  m_pAmod <- m_pArea
  Area <- m_pArea
  ext <- dim(m_pArea)

  nChanges <- 1
  Iteration <- 0

  while(nChanges > 0)
  {
    Iteration <- Iteration + 1
    nChanges <- 0

    for (row in 1:ext[1])
    {
      for (col in 1:ext[2])
      {
        if (m_Suction[row,col] != nodata)
        {
          z <- m_Suction[row,col] * Get_Local_Maximum(Area, row, col, nodata)
          if (z > Area[row,col])
          {
            nChanges <- nChanges + 1
            Area[row,col] <- z
          }
        }
      }
    }

    if (nChanges > 0)
    {
      nChanges <- 0

      for (row in 1:ext[1])
      {
        for (col in 1:ext[2])
        {
          if (Area[row,col] != m_pAmod[row,col])
          {
            nChanges <- nChanges + 1
            m_pAmod[row,col] <- Area[row,col]
          }
        }
      }
    }

    cat("\npass", Iteration, "(", nChanges, "> 0)")
  }

  cat("\npost-processing...")
  for (row in 1:ext[1])
  {
    for (col in 1:ext[2])
    {
      if (Area[row,col] != nodata)
      {
        bModify <- FALSE
        n <- 0
        z <- 0

        for (drow in -1:1)
        {
          rown <- row + drow
          for (dcol in -1:1)
          {
            coln <- col + dcol
            
            is_InGrid_To <- ifelse(
              (rown <= ext[1] & rown >= 1) & (coln <= ext[2] & coln >= 1),
              ifelse(m_pArea[rown, coln] != nodata, TRUE, FALSE),
              FALSE
            )
            
            if (is_InGrid_To)
            {
              if (Area[rown,coln] > m_pArea[rown,coln]) bModify <- TRUE
              n <- n + 1
              z <- z + Area[rown,coln]
            }
          }
        }
        m_pAmod[row,col] <- ifelse(bModify, z / n, Area[row,col])
      } else {
        m_pAmod[row,col] <- nodata
      }
    }
  }
  
  return(m_pAmod)
}


## --
# Calcul du TWI
Get_TWI <- function(m_pAmod, m_pSlope, area_type, slope_type, slope_min, slope_offset, cellsize)
{
  ext <- dim(m_pAmod)
  m_pTWI <- m_pAmod
  slope_min <- slope_min * pi / 180
  slope_offset <- slope_offset * pi / 180

  for (row in 1:ext[1])
  {
    for (col in 1:ext[2])
    {
      if (m_pAmod[row,col] == nodata | m_pSlope[row,col] == nodata)
      {
        m_pTWI[row,col] <- nodata
      } else {
        if (slope_type == "catchment slope")
        {
          Slope	<- m_pSlope[row,col]
        } else {
          Slope <- Get_Gradient(m_pDEM, row, col, nodata, cellsize)[["Slope"]]
        }

        Slope	<- tan(max(slope_min, Slope + slope_offset))
        Area <- m_pAmod[row,col]

        if (area_type == "square root of catchment area")
        {
          Area <- sqrt(Area) 
        } else if (area_type == "specific catchment area") {
          Area <- Area / cellsize
        }

        m_pTWI[row,col] <- log(Area / Slope)
      }
    }
  }

  return(m_pTWI)
}


## --
# Fonction principale
Get_STWI <- function(rDEM, suction, area_type, slope_type, slope_weight, slope_min, slope_offset, nodata, cellsize, rWeight = NULL)
{
  # Conversion des objets raster en matrices pour les traitements
  m_pDEM <- as.matrix(rDEM)
  m_pDEM[is.na(m_pDEM)] <- nodata
  
  if (is.null(rWeight))
  {
    m_pWeight <- m_pDEM
    m_pWeight[m_pDEM != nodata] <- 1
  }
  else
  {
    m_pWeight <- as.matrix(rWeight)
    m_pWeight[is.na(m_pWeight)] <- nodata
  }
  
  
  # Vérifie que le NoData de m_pWeight est au même endroit que celui de m_pDEM
  if ( !all( (as.vector(m_pDEM) == nodata) == (as.vector(m_pWeight) == nodata) ) )
    stop("Le NoData de rWeight doit être aux mêmes endroits que celui de rDEM", call. = FALSE)
  
  
  # Conversion des types d'accumulation et de pente
  area_type <- c("total catchment area", "square root of catchment area", "specific catchment area")[area_type + 1]
  slope_type <- c("local slope", "catchment slope")[slope_type + 1]
  
 
  # Traitement lui-même
  ls_area <- Get_Area(m_pDEM, m_pWeight, suction, slope_weight, nodata, cellsize)
	m_pAmod <- Get_Modified(ls_area[["m_pArea"]], ls_area[["m_pSlope"]], ls_area[["m_Suction"]], nodata)
	m_pTWI <- Get_TWI(m_pAmod, ls_area[["m_pSlope"]], area_type, slope_type, slope_min, slope_offset, cellsize)
	
	
	# Conversion de la matrice de TWI en objet raster
	rTWI <- rDEM
	rTWI[] <- as.vector(t(m_pTWI))
	rTWI[rTWI == nodata] <- NA

  return(rTWI)
}
