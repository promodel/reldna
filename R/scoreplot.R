scoreplot <- function#Function to make scoreplots
###
###
(resPLS, ##<< mvr-object
comps = c(1, 2),  ##<< numeric vector of length equals two, numbers of components for scoreplot
groups, ##<< numeric vector or list, containing indices of observation to show in different colors (optional)
col, ##<< if groups are specified, it is possible to choose colors for them (optional)
ellipse = TRUE, ##<< logical, if TRUE: ellipse that contains 95% of the observed data is drawn
labels = TRUE, ##<< logical, if TRUE: observations are labeled 
pch = 20, ##<< integer, specifies symbol type in plot 
save = FALSE,  ##<< logical, to save plot or not?
filetype, ##<< if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.
filename='PLS_scores', ##<< character string, if save=TRUE specifies a name of the file
show = TRUE ##<< logical, to show plot or not?
) {
		makegroups <- TRUE
    gtest <- try(length(groups), TRUE)
    if (class(gtest) == "try-error") {
        makegroups <- FALSE
    }
    ctest <- try(length(col), TRUE)
    if (class(ctest) == "try-error") {
        col <- palette()
    }
 sp <- function() {
            if (ellipse == TRUE) {
                dataEllipse <- car::dataEllipse  # maybe some other way?
                dataEllipse(resPLS$scores[, comps[1]], resPLS$scores[, comps[2]], 
                  levels = 0.95, lty = 1, fill = FALSE, pch = pch, main = "Scores", 
                  xlab = paste("Comp", comps[1], sep = " "), ylab = paste("Comp", 
                    comps[2], sep = " "), col = "black")
                abline(v = 0, h = 0)
                if (makegroups == TRUE) {
                  grclass <- is(groups)
                  if (grclass[1] == "list") {
                    for (i in 1:length(groups)) {
                      points(resPLS$scores[groups[[i]], comps[1]], resPLS$scores[groups[[i]], 
                        comps[2]], pch = pch, col = col[i])
                    }
                  } else {
                    points(resPLS$scores[groups, comps[1]], resPLS$scores[groups, 
                      comps[2]], pch = pch, col = col[1])
                    points(resPLS$scores[-groups, comps[1]], resPLS$scores[-groups, 
                      comps[2]], pch = pch, col = col[2])
                  }
                }
                if (labels == TRUE) {
                  text(resPLS$scores[, comps[1]], resPLS$scores[, comps[2]], pos = 3, 
                    names(resPLS$scores[, comps[1]]))
                }
            }
            if (ellipse == FALSE) {
                if (exists("groups") == TRUE) {
                  grclass <- is(groups)
                  if (grclass[1] == "list") {
                    for (i in 1:length(groups)) {
                      points(resPLS$scores[groups[[i]], comps[1]], resPLS$scores[groups[[i]], 
                        comps[2]], pch = pch, col = col[i])
                    }
                  } else {
                    plot(resPLS$scores[, comps[1]], resPLS$scores[, comps[2]], type = "p", 
                      main = "Scores", xlab = paste("Comp", comps[1], sep = " "), 
                      ylab = paste("Comp", comps[2], sep = " "), pch = pch, col = col[1])
                    points(resPLS$scores[groups, comps[1]], resPLS$scores[groups, 
                      comps[2]], pch = pch, col = col[1])
                    points(resPLS$scores[-groups, comps[1]], resPLS$scores[-groups, 
                      comps[2]], pch = pch, col = col[2])
                  }
                }
                if (labels == TRUE) {
                  text(resPLS$scores[, comps[1]], resPLS$scores[, comps[2]], pos = 3, 
                    names(resPLS$scores[, comps[1]]))
                }
            }
        }
    if (show == TRUE) {
       sp()
    }
    
    # save pic
    if (save == TRUE) {
       filesave(sp(), filetype, filename)
    }
}  
