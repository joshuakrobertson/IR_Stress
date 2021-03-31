## COLI Thermal Responses to Stress Exposure

#------------------------------------------------------------------------------#
## Loading in necessary packages, fonts, and plotting themes.
#------------------------------------------------------------------------------#

library('easypackages')
libraries('curl','emmeans','grid','gridExtra','itsadug',
    'MASS','mgcv','showtext','tidyverse','viridis')

font_add_google(name = "Noto Sans", family = "Noto Sans")
showtext_auto()

my.theme = theme(
  panel.grid.minor = element_blank(),
  axis.title = element_text(size = 18, family = "Noto Sans"),
  axis.text = element_text(size = 14, colour = "black", family = "Noto Sans"),
  axis.title.y = element_text(vjust = 0.5), panel.grid.major =
    element_line(colour = "grey75"), legend.title = element_text(
    size = 14,
    colour = "black", family = "Noto Sans"
  ), legend.text = element_text(
    size = 14,
    colour = "black", family = "Noto Sans"
  )
)

#------------------------------------------------------------------------------#
# Adding additional functions
#------------------------------------------------------------------------------#

Dotplot = function(data, y, fill = "gold2", colour = "black", pch = 21, size = 2, error_sd = 4){
    require('ggplot2')

    if(!exists("data")){
        return("No data file provided")
    }
    if(!is.data.frame(data)){
        return("data argument must be a dataframe object")
    }

    if(!(y %in% colnames(data))){
        return("y argument cannot be found in the supplied data frame")
    }

    if(!is.integer(error_sd) & !is.numeric(error_sd)){
        return("error_sd must be numeric or integerial")
    }

    if(!is.integer(pch) & !is.numeric(pch)){
        return("pch must be an integer")
    }

    if(!is.character(fill)){
        return("fill must be a character string")
    }

    if(!is.character(colour)){
        return("colour must be a character string")
    }

    dplot = data %>% ggplot(aes(x = 1:nrow(.), y = get(y))) + 
        geom_point(size = size, pch = pch, fill = fill, colour = colour) + 
        annotate(geom = "rect", xmin = 1, xmax = nrow(data), 
        ymin = mean(data[,y], na.rm = T) - 
            error_sd*sd(data[,y], na.rm = T),
        ymax = mean(data[,y], na.rm = T) + 
            error_sd*sd(data[,y], na.rm = T),
        fill = "grey20", alpha = 0.5) + 
    theme_bw() + xlab("Row Number") + ylab(y)    

    return(dplot)
}

ACF_Plot = function(model, residual_type = "deviance", colour = "black", fill = "dodgerblue2"){
    
    require('ggplot2')

    if (residual_type == "deviance" | residual_type == "response" | residual_type == "pearson" | 
        residual_type == "scaled.pearson" | residual_type == "working") { 
    
    ACF_Dat = data.frame("Lag" = c(acf(residuals(model, type = residual_type), plot = F)$lag),
        "acf" = c(acf(residuals(model, type = residual_type), plot = F)$acf))

    acf_plot = ACF_Dat %>% ggplot(aes(x = Lag, y =  acf)) + 
        geom_col(fill = fill, colour = colour) + 
        theme_bw()

    return(acf_plot)   

    } else if (residual_type == "normalized" | residual_type == "normalised"){
    
        if (class(model)[1] == "gam" | class(model)[1] == "bam"){
           require('itsadug')
    
           ACF_Dat = data.frame("Lag" = c(acf(resid_gam(model), plot = F)$lag),
            "acf" = c(acf(resid_gam(model), plot = F)$acf))

           acf_plot = ACF_Dat %>% ggplot(aes(x = Lag, y =  acf)) + 
            geom_col(fill = fill, colour = colour) + 
           theme_bw()

           return(acf_plot)    
 
        } else if (class(model)[1] == "lme"){
 
           if (residual_type == "normalised"){
               residual_type = "normalized"
           } else {
           }

           ACF_Dat = data.frame("Lag" = c(acf(residuals(model, type = residual_type), plot = F)$lag),
           "acf" = c(acf(residuals(model, type = residual_type), plot = F)$acf))

           acf_plot = ACF_Dat %>% ggplot(aes(x = Lag, y =  acf)) + 
           geom_col(fill = fill, colour = colour) + 
           theme_bw()

           return(acf_plot)  

        }   else {
            
            return("Normalised/Normalized residuals cannot be calculated for this model type")
        }
    
    } else {
        
    return("Residual type is not supported")

    }
}

HeatT = function(Ts,Ta,atm = 101.325, dim = 2*sqrt(SA/pi), SA, Emis = 0.95){
  
  if (is.numeric(dim) == FALSE){
    return ("dim must be a numeric value equal to the horizontal length of the object, in m.")
  }
  if (is.numeric(SA) == FALSE){
    return ("SA must be numeric value equal to the surface area of the object in m2.")
  }

  require('bigleaf')
  require('dplyr')
  require('Thermimage')
  require('doParallel')
  registerDoParallel()

  if (length(Ta) != length(Ts)){
    if (length(Ta) != 1){
      return("Ta must be of length 1, or equal to the length of Ts")
    } else if (length(Ta) == 1){
      Ta = c(rep(Ta, length(Ts)))
    }
  }

   TD = data.frame(Ta = Ta, Ts = Ts, ID = seq(from = 1, to = length(Ts), by = 1))

   if (identical(c(which(is.na(Ts))), integer(0)) == TRUE){
     NAs_TR = TD  
   } else {
    NAs_TR = TD[c(which(is.na(Ts))),]
   }

   CTD = na.omit(TD)

    Kt = c()
    Vis = c()
    Alpha = c()
    Gr = c()
    Nu = c()
    Hc = c()

    Kt = unlist(foreach(i=1:length(CTD$Ta)) %dopar% {
        airtconductivity(CTD$Ta[i])
    })

    Vis = unlist(foreach(i=1:length(CTD$Ta)) %dopar% {
        airtconductivity(CTD$Ta[i])
    })
    
    Alpha = unlist(foreach(i=1:length(CTD$Ta)) %dopar% {
        airtconductivity(CTD$Ta[i])
    })
    
    Gr = unlist(foreach(i=1:length(CTD$Ta)) %dopar% {
        ((Alpha[i]) * (9.81) * (dim)^3 *
        (CTD$Ts[i] - CTD$Ta[i])) / (Vis[i])^2
    })

    Nu = unlist(foreach(i=1:length(CTD$Ts)) %dopar% {
        sign(Gr[i]) * 0.50 * (abs(Gr[i]))^0.25
      })

    Hc = unlist(foreach(i=1:length(CTD$Ts)) %dopar% {
        Nu[i] * (Kt[i] / dim)
    })

    Qconv = unlist(foreach(i=1:length(CTD$Ts)) %dopar% {
        SA * (Hc[i]) * (CTD$Ts[i] - CTD$Ta[i])
    })

    Qrad = unlist(foreach(i=1:length(CTD$Ts)) %dopar% {
        SA * (5.67 * 10^-8) * Emis^2*
        ((CTD$Ts[i] + 273.15)^4 - (CTD$Ta[i] + 273.15)^4)
    })

    Qtot = unlist(foreach(i=1:length(CTD$Ts)) %dopar% {
        Qconv[i] + Qrad[i]
    })

    CTD$Qtot = Qtot
    NAs_TR$Qtot = NA

    if (identical(c(which(is.na(Ts))), integer(0)) == TRUE){
     to_return = CTD  
   } else {
    to_return = rbind(CTD, NAs_TR)
    to_return = dplyr::arrange(to_return, ID)
    }
    return(to_return$Qtot)
}

Group_Diag = function(model, data, group, view, label = NA, cutoff = 2.5, fill = "grey70", method = "gam", bs = "cr", k = 3){

    require('ggplot2')
    require('grid')
    require('gridExtra')
    require('dplyr')

    if(class(model)[1] == "gamm") {

    Response = suppressWarnings(gsub("[[:space:]]", "", summary(model$gam)$formula)[[2]])

    Pred_Dat = as.data.frame(data %>% 
      mutate("Pred" = predict(model$gam, newdata = data, type = "response"), 
             "Res" = residuals(model$lme, type = "normalized")))

    # Predict curves 

    Pred_Dat$Res_Wide = Pred_Dat$Pred
    Pred_Dat$Res_Wide[c(which(abs(Pred_Dat$Res) <= cutoff))] = NA

    if (is.na(label)){

    if (method == "gam"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, group = Pred_Dat[,group])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +        
        theme_bw()
    } else if (method == "lm"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +            
        theme_bw() 
    } else if (method == "loess"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +              
        theme_bw() 
    } else {
        return("method must be either gam, loess, or lm")
    }

    # Raw curves

    Pred_Dat$Res_Wide = Pred_Dat[,Response]
    Pred_Dat$Res_Wide[c(which(abs(Pred_Dat$Res) <= cutoff))] = NA
    
    if (method == "gam"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +        
        theme_bw()
    } else if (method == "lm"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() 
    } else if (method == "loess"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() 
    } else {
        return("method must be either gam, loess, or lm")
    }

    Joined_Plots = grid.arrange(Pred_Plot, Res_Plot, nrow = 1, top = "Group-wise Outlier Assessment")
    return(Joined_Plots)

    } else if (!is.na(label)){

    if (method == "gam"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +        
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "lm"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +            
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "loess"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +              
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else {
        return("method must be either gam, loess, or lm")
    }

    # Raw curves

    Pred_Dat$Res_Wide = Pred_Dat[,Response]
    Pred_Dat$Res_Wide[c(which(abs(Pred_Dat$Res) <= cutoff))] = NA
    
    if (method == "gam"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +        
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "lm"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "loess"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else {
        return("method must be either gam, loess, or lm")
    }

    Joined_Plots = grid.arrange(Pred_Plot, Res_Plot, nrow = 1, top = "Group-wise Outlier Assessment")
    return(Joined_Plots)
    }

    } else if (class(model)[1] == "gam" | class(model)[1] == "bam"){

    require('itsadug')

    Response = suppressWarnings(gsub("[[:space:]]", "", summary(model)$formula)[[2]])

    Pred_Dat = as.data.frame(data %>% 
      mutate("Pred" = predict(model, newdata = data, type = "response"), 
             "Res" = resid_gam(model, incl_na = TRUE)))

    # Predict curves 

    Pred_Dat$Res_Wide = Pred_Dat$Pred
    Pred_Dat$Res_Wide[c(which(abs(Pred_Dat$Res) <= cutoff))] = NA

    if (is.na(label)){

    if (method == "gam"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, group = Pred_Dat[,group])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +        
        theme_bw()
    } else if (method == "lm"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +            
        theme_bw() 
    } else if (method == "loess"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +              
        theme_bw() 
    } else {
        return("method must be either gam, loess, or lm")
    }

    # Raw curves

    Pred_Dat$Res_Wide = Pred_Dat[,Response]
    Pred_Dat$Res_Wide[c(which(abs(Pred_Dat$Res) <= cutoff))] = NA
    
    if (method == "gam"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +        
        theme_bw()
    } else if (method == "lm"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() 
    } else if (method == "loess"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() 
    } else {
        return("method must be either gam, loess, or lm")
    }

    Joined_Plots = grid.arrange(Pred_Plot, Res_Plot, nrow = 1, top = "Group-wise Outlier Assessment")
    return(Joined_Plots)

    } else if (!is.na(label)){

    if (method == "gam"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +        
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "lm"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +            
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "loess"){
        Pred_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat$Pred, Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Predicted Values") + xlab(view) +              
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else {
        return("method must be either gam, loess, or lm")
    }

    # Raw curves

    Pred_Dat$Res_Wide = Pred_Dat[,Response]
    Pred_Dat$Res_Wide[c(which(abs(Pred_Dat$Res) <= cutoff))] = NA
    
    if (method == "gam"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, formula = y~s(x, bs = bs, k = k), 
        size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +        
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "lm"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else if (method == "loess"){
        Res_Plot = ggplot(Pred_Dat, aes(x = Pred_Dat[,view], y = Pred_Dat[,Response], group = Pred_Dat[,group],
            linetype = Pred_Dat[,label])) + 
        geom_smooth(method = method, size = 0.5, colour = "grey50") + 
        geom_point(aes(y = Res_Wide), size = 3, pch = 21, colour = "black", 
        fill = fill) + ylab("Raw Values") + xlab(view) +           
        theme_bw() + guides(linetype=guide_legend(title="Label"))
    } else {
        return("method must be either gam, loess, or lm")
    }

    Joined_Plots = grid.arrange(Pred_Plot, Res_Plot, nrow = 1, top = "Group-wise Outlier Assessment")
    return(Joined_Plots)
    }     
    }    
}

Coef_Dist = function(model, incl_int = FALSE, incl_re = FALSE, facet = TRUE){

    require('ggplot2')
    require('grid')
    require('gridExtra')

    if (class(model)[1] == "gam" | class(model)[1] == "bam"){

    Coefs = data.frame("Predictors" = names(model$coefficients),
     "Coefficients" = as.numeric(model$coefficients),
     "E.D.F." = as.numeric(model$edf),
     "S.E.M" = as.numeric(sqrt(diag(vcov(model, unconditional = TRUE))))
    )

    Preds = unlist(str_split(as.character(summary(model)$formula)[3], pattern = "\\+"))
    Pred_Type = c() 

    for (i in 1:length(Preds)){

        if (identical(grep("s\\(", Preds[i]), integer(0)) == TRUE){
            Pred_Type[i] = "Parametric" 
        } else if (identical(grep("s\\(", Preds[i]), integer(0)) == FALSE) {
            if (identical((grep("\"re", unlist(str_split(Preds[i], pattern = "\\,")))), integer(0)) == TRUE & 
            identical((grep("\"fs", unlist(str_split(Preds[i], pattern = "\\,")))), integer(0)) == TRUE){
            Pred_Type[i] = "Fixed" 
        } else {
            Pred_Type[i] = "Random" 
        } 
      }
    }

    Pred_Type = c("Parametric", Pred_Type)
    Coefs$Pred_Whole = gsub("\\.[[:digit:]]*.", "", Coefs$Predictors)

    Pairing = data.frame("Pred_Whole" = unique(Coefs$Pred_Whole), "Type" = Pred_Type)
    Coefs = left_join(Coefs, Pairing, by = "Pred_Whole")

    if (incl_re == FALSE){

    Coefs = Coefs %>% filter(Type != "Random")
    
    } else if (incl_re == TRUE){
      message("Random effects included")  
    } else {
        return("incl_re must be TRUE or FALSE")
    }

    if (incl_int == FALSE){

       Coefs = Coefs[-1,]

    } else if (incl_int == TRUE){
      message("Intercept included")  
    } else {
        return("incl_int must be TRUE or FALSE")
    }

    Coef_Plots = vector('list', nrow(Coefs))
    Coefs$LCL = Coefs$Coefficients - 2.576*Coefs$S.E.M
    Coefs$UCL = Coefs$Coefficients + 2.576*Coefs$S.E.M

    Range = c(min(Coefs$LCL), max(Coefs$UCL))

        for (i in 1:nrow(Coefs)){

        Coef_Plots[[i]] = ggplot(data = data.frame(x = Range), aes(x)) +    
            stat_function(fun = dnorm, n = 100, args = list(mean = Coefs$Coefficients[i], 
            sd = Coefs$S.E.M[i]), geom = "area", alpha = 0.4, fill = viridis::viridis(n = 20)[i]) +   
            ylab(Coefs$Predictors[i]) + scale_y_continuous(breaks = NULL) + theme_bw() + 
            theme(axis.title.y = element_text(angle = 0, vjust=0.5)) + xlab("Coefficient") + 
            xlim(Range)   
        }

       if (facet == FALSE){
          return(Coef_Plots)
       } else if (facet == TRUE){
          All_Plots = do.call("grid.arrange", c(Coef_Plots, ncol=1))
          return(All_Plots)
       } else {
         return("facet must be TRUE or FALSE")
       }

    } else {
        return("Model class currently not supported")
    }
}

post_simGam = function(model, data, view, group_var, sim = 100, greyscale = TRUE, seed = 12){

    require('dplyr')
    require('ggplot2')
    require('MASS')
    require('viridis')

    rows_tk = seq(1, nrow(data))

    Mod_Terms = gsub("\\.[[:digit:]].*", "", gsub("\\)", "", 
        gsub("s\\(", "", names(model$coefficients))))

    if (identical(grep(":", Mod_Terms), integer(0)) == TRUE){
         Mod_Terms = Mod_Terms[-1]    
    } else {
       throw = grep(":", Mod_Terms)
       Mod_Terms = Mod_Terms[-c(1,throw)] 
    } 

    if (identical(grep(".L", Mod_Terms), integer(0)) == TRUE & 
        identical(grep(".U", Mod_Terms), integer(0)) == TRUE ){
       Mod_Terms = Mod_Terms
    } else {
       Mod_Terms = gsub(".U", "", gsub(".L", "", Mod_Terms))
    }

    Mod_Terms = unique(Mod_Terms)

    # Creating prediction data

    pred_dat = data[,Mod_Terms]

    lp_mat = predict(model, newdata = pred_dat, type = "lpmatrix")
    coefficients = model$coefficients
    var_covar = vcov(model)

    set.seed(seed)
    
    sim_gam = mvrnorm(sim, mu = coefficients, Sigma = var_covar)
    keep_gam = grep(view, colnames(lp_mat))
    mod_sims_fit = lp_mat[, keep_gam] %*% t(sim_gam[, keep_gam])

    ylims = range(mod_sims_fit)

    # Colouring by group

    Long_Pred = vector('list', ncol(mod_sims_fit))

        for (i in 1:length(Long_Pred)){
            Long_Pred[[i]] = pred_dat %>% mutate("Pred" = mod_sims_fit[,i], "Sim_Number" = i)
        }

    Long_Pred = bind_rows(Long_Pred)
    Long_Pred$Group = paste(Long_Pred[,group_var], Long_Pred$Sim_Number, sep = "_")

    if (greyscale == TRUE){

        Sim_Plot = ggplot(Long_Pred, aes(x = Long_Pred[,view], y = Pred, group = Group, colour = Long_Pred[,group_var])) + 
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1, se = FALSE, alpha = 0.3) + 
        theme_bw() + scale_colour_manual(values = c("black", "grey50")) + 
        ylab("Centered Response") + xlab(view) + guides(colour=guide_legend(title=group_var))

    } else if (greyscale == FALSE){

        Sim_Plot = ggplot(Long_Pred, aes(x = Long_Pred[,view], y = Pred, group = Group, colour = Long_Pred[,group_var])) + 
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3), size = 1, se = FALSE, alpha = 0.3) + 
        theme_bw() + scale_colour_viridis_d() + 
        ylab("Centered Response") + xlab(view) + guides(colour=guide_legend(title=group_var))

    }

    return(Sim_Plot)
} 

#------------------------------------------------------------------------------#
# Loading, assessing, and cleaning data
#------------------------------------------------------------------------------#

setwd("/home/joshk/git_repositories/COLI_Thermal/Data")
Therm_Dat = read.csv("COLI_Surface_Temperature_Data.csv")

# Binding in sex data 

Sex = read.csv("COLI_Sex.csv")
Sex = Sex %>% rename("ID" = TZ_ID)
Therm_Dat = left_join(Therm_Dat, Sex, by = c("ID"))

# Cutting data to include only measurements captured within 420 seconds 
# before handling, and 420 seconds after handling.

Therm_Dat = filter(Therm_Dat, Slice <=840)

# Removing measurments with poor image quality.

Therm_Dat$Eye.Temp[c(which(Therm_Dat$Eye.Quality != "Good"))] = NA
Therm_Dat$Bill.Temp[c(which(Therm_Dat$Bill.Quality != "Good"))] = NA

# Counting total frames used and total individuals used

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>%
    summarise(count = n())

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>%
    group_by(Eye.Quality) %>% summarise(count = n())

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>%
    group_by(Bill.Quality) %>% summarise(count = n())

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>% 
    group_by(ID, Sex) %>% summarise(n())

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>%
    group_by(Treatment) %>% summarise("Count" = n())

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>% 
    group_by(Treatment, Sex, ID) %>% summarise(n()) %>% 
    group_by(Treatment, Sex) %>% summarise(Count = n())

Therm_Dat %>% filter(Eye.Quality == "Good" | Bill.Quality == "Good") %>% 
    group_by(Treatment, Sex, ID) %>% summarise("Count" = n()) %>% 
    group_by(Treatment) %>% summarise("Mean" = mean(Count), "SD" = sd(Count))

# Assessing distribution of data points using Cleveland Dot-plots

# Eyes

Dotplot(data = Therm_Dat, y = "Eye.Temp", fill = "dodgerblue3")

# A few eye temperature values well below 4 standard deviations. Sounting and replacing these values with NAs.

length(Therm_Dat$Eye.Temp[c(which(Therm_Dat$Eye.Temp <= 30))])
Therm_Dat$Eye.Temp[c(which(Therm_Dat$Eye.Temp <= 30))] = NA

#Bill

Dotplot(data = Therm_Dat, y = "Bill.Temp", fill = "dodgerblue3")

# Three extreme outliers. Counting and replacing values with NAs.

length(Therm_Dat$Bill.Temp[c(which(Therm_Dat$Bill.Temp <= 20 | Therm_Dat$Bill.Temp >= 60))])
Therm_Dat$Bill.Temp[c(which(Therm_Dat$Bill.Temp <= 20 | Therm_Dat$Bill.Temp >= 60))] = NA

#------------------------------------------------------------------------------#
# Calculating heat-transfer across images. Surface area of bill is estimated according to
# measurements drawn from: Johnston, R. F. (1992). Evolution in the rock dove: skeletal 
# morphology. The Auk, 109(3), 530-542 and formulas adapted from: Greenberg, R., Cadena, 
# V., Danner, R. M., & Tattersall, G. (2012). Heat loss may explain bill size differences 
# between birds occupying different habitats. PloS one, 7(7).

# NOTE that heat transfer analyses are not included in the manuscript accompanying this 
# R file.
#------------------------------------------------------------------------------#

# Males
## Bill length = 23.7 mm
## Bill width = 3.3 mm

# Females
## Bill length = 23.7 mm
## Bill width = 3.1 mm

Bill_Width = 3.2
Bill_Length = 23.7
Bill_Depth = Bill_Width

# Surface area formula 

Bill_SA = (((Bill_Width+Bill_Depth)/4)*(Bill_Length/3)*pi) + 
   (pi*sqrt(2*((Bill_Width/2)^2 + (Bill_Depth/2)^2) - 
   0.5*(Bill_Width-Bill_Depth)^2))*((Bill_Length/3)*2)*2

# Note that two thirds of the bill length is treated as a cylinder, and 
# one third of the bill is treated a cone. Width and depth are assumed
# to be equal.

# Eye diameter estimated as 1.3 cm from: 
# Donovan, W. J. (1978). Structure and function of the pigeon visual system. 
# Physiological Psychology, 6(4), 403-437.

Eye_Diameter = 13 
Eye_SA = (pi*(Eye_Diameter/2)^2)*2

# Note that eye surface area is multiplied by two to account for
# heat at both eyes

# Converting measurements to metres and calculating heat transfer.

Bill_SA = Bill_SA*(1*10^-6)
Eye_SA = Eye_SA*(1*10^-6)
Bill_Length = Bill_Length*(1*10^-3)
Eye_Diameter = Eye_Diameter*(1*10^-3)

Therm_Dat$Eye_HT = HeatT(Ts = Therm_Dat$Eye.Temp, Ta = 18, dim = Eye_Diameter, SA = Eye_SA)
Therm_Dat$Bill_HT = HeatT(Ts = Therm_Dat$Bill.Temp, Ta = 18, dim = Bill_Length, SA = Bill_SA)
Therm_Dat$qTot = Therm_Dat$Bill_HT + Therm_Dat$Eye_HT

# Plotting time series'

{
    EyeTemp_TS = ggplot(Therm_Dat, aes(x = Slice, y = Eye.Temp, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Eye Temperature (°C)") + 
    theme(legend.position = "none")   

    BillTemp_TS = ggplot(Therm_Dat, aes(x = Slice, y = Bill.Temp, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Bill Temperature (°C)") + 
    theme(legend.position = "none") 

    EyeHT_TS = ggplot(Therm_Dat, aes(x = Slice, y = Eye_HT, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Eye Heat Transfer (W)") + 
    theme(legend.position = "none")     

    BillHT_TS = ggplot(Therm_Dat, aes(x = Slice, y = Bill_HT, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Bill Heat Transfer (W)") + 
    theme(legend.position = "none") 

    TotalHT_TS = ggplot(Therm_Dat, aes(x = Slice, y = qTot, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Total Heat Transfer (W)")
}

lay = rbind(c(1,2,5,5),
             c(3,4,5,5))    

Raw_Plots = grid.arrange(EyeTemp_TS, BillTemp_TS, EyeHT_TS, BillHT_TS, TotalHT_TS, 
    layout_matrix = lay, top = "Raw Time Series Trends")

Raw_Plots

# Averaging measurements across every individual at each 5 second interval

{
    EyeTemp_TS = ggplot(Therm_Dat, aes(x = Slice, y = Eye.Temp, 
    fill = Treatment, linetype = Treatment)) + 
    stat_summary_bin(fun = "mean", binwidth = 5, geom = "point", size = 2, 
        pch = 21, colour = "black", alpha = 0.5) + 
    stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, geom = "ribbon",
        alpha = 0.5) +     
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Eye Temperature (°C)") + 
    theme(legend.position = "none") + 
    geom_vline(xintercept = 420, size = 1, colour = "grey30", linetype = "longdash") +
    annotate("text", x=-Inf, y = Inf, label = "A", vjust = 1.1, hjust = 0, size = 8)

    BillTemp_TS = ggplot(Therm_Dat, aes(x = Slice, y = Bill.Temp, 
    fill = Treatment, linetype = Treatment)) + 
    stat_summary_bin(fun = "mean", binwidth = 5, geom = "point", size = 2, 
        pch = 21, colour = "black", alpha = 0.5) + 
    stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, geom = "ribbon",
        alpha = 0.5) +     
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Bill Temperature (°C)") + 
    theme(legend.position = "none") + 
    geom_vline(xintercept = 420, size = 1, colour = "grey30", linetype = "longdash") +
    annotate("text", x=-Inf, y = Inf, label = "B", vjust = 1.1, hjust = 0, size = 8)    

    EyeHT_TS = ggplot(Therm_Dat, aes(x = Slice, y = Eye_HT, 
    fill = Treatment, linetype = Treatment)) + 
    stat_summary_bin(fun = "mean", binwidth = 5, geom = "point", size = 2, 
        pch = 21, colour = "black", alpha = 0.5) + 
    stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, geom = "ribbon",
        alpha = 0.5) +     
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Eye Heat Transfer (W)") + 
    theme(legend.position = "none") + 
    geom_vline(xintercept = 420, size = 1, colour = "grey30", linetype = "longdash") +
    annotate("text", x=-Inf, y = Inf, label = "C", vjust = 1.1, hjust = 0, size = 8)        

    BillHT_TS = ggplot(Therm_Dat, aes(x = Slice, y = Bill_HT, 
    fill = Treatment, linetype = Treatment)) + 
    stat_summary_bin(fun = "mean", binwidth = 5, geom = "point", size = 2, 
        pch = 21, colour = "black", alpha = 0.5) + 
    stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, geom = "ribbon",
        alpha = 0.5) +     
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Bill Heat Transfer (W)") + 
    theme(legend.position = "none") + 
    geom_vline(xintercept = 420, size = 1, colour = "grey30", linetype = "longdash") +
    annotate("text", x=-Inf, y = Inf, label = "D", vjust = 1.1, hjust = 0, size = 8)    

    TotalHT_TS = ggplot(Therm_Dat, aes(x = Slice, y = qTot, 
    fill = Treatment, linetype = Treatment)) + 
    stat_summary_bin(fun = "mean", binwidth = 5, geom = "point", size = 2, 
        pch = 21, colour = "black", alpha = 0.5) + 
    stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, geom = "ribbon",
        alpha = 0.5) +     
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3, bs = "cr"),
        colour = "black", size = 1, se = FALSE) + 
    scale_fill_manual(values = c("#00CC99", "#6633CC")) +
    theme_bw() + xlab("Frame") + ylab("Total Heat Transfer (W)") + 
    geom_vline(xintercept = 420, size = 1, colour = "grey30", linetype = "longdash") +
    annotate("text", x=-Inf, y = Inf, label = "E", vjust = 1.1, hjust = 0, size = 8)   
}

lay = rbind(c(1,2,5,5),
             c(3,4,5,5))    

Av_Plots = grid.arrange(EyeTemp_TS, BillTemp_TS, EyeHT_TS, BillHT_TS, TotalHT_TS, 
    layout_matrix = lay)

Av_Plots

# Note that trends in heat transfer rates should appear identical to those in surface temperature
# given that environmental conditions were held constant across birds.

#------------------------------------------------------------------------------#
## Modelling eye temperature responses
#------------------------------------------------------------------------------#

# Assessing distribution of possible response variables

h1 = ggplot(Therm_Dat, aes(x = Eye.Temp)) + 
    geom_histogram(binwidth = 0.5, colour = "black", fill = "dodgerblue3")+
      geom_density(aes(y=0.5 * ..count..), 
      colour = "black", fill = "dodgerblue3", alpha = 0.5) + 
      xlab("Eye Temperature (°C)")

h2 = ggplot(Therm_Dat, aes(x = Bill.Temp)) + 
    geom_histogram(binwidth = 1, colour = "black", fill = "magenta")+
      geom_density(aes(y=1 * ..count..), 
      colour = "black", fill = "magenta", alpha = 0.5) + 
      xlab("Bill Temperautre (°C)")    

h3 = ggplot(Therm_Dat, aes(x = qTot)) + 
    geom_histogram(binwidth = 0.003, colour = "black", fill = "grey70")+
      geom_density(aes(y=0.003 * ..count..), 
      colour = "black", fill = "black", alpha = 0.5) + 
      xlab("Total Heat Transfer (W)")

lay = rbind(c(1,3,3),
            c(2,3,3))    

grid.arrange(h1, h2, h3, 
    layout_matrix = lay, top = "Response variable distributions")

# Irregular. Eye and bill temperature distributions appear bimodal, perhaps due to
# differences by treatment. Proceeding with models.

Therm_Dat$OTreatment = factor(Therm_Dat$Treatment, ordered = T)
Therm_Dat$ID = factor(Therm_Dat$ID)

EyeTemp_Mod = gam(Eye.Temp ~ OTreatment +
    s(Slice, bs = "cr", k = 3) + 
    s(Slice, by = OTreatment, bs = "cr", k = 3, m = 1) + 
    s(ID, bs = "re"),  
    method = "REML", na.action = na.omit, data = Therm_Dat)

# Assessing temporal autocorrelation

ACF_Plot(model = EyeTemp_Mod, residual_type = "normalised")
acf(residuals(EyeTemp_Mod), plot = FALSE)$acf[2]

# Very high autocorrelation. Adjusting using an AR1 structure in bam.
# Fist binning on ten second intervals to speed model convergence.

bin_Fun = function(x) floor(x/10)*10
Therm_Dat$Time_Bin = bin_Fun(Therm_Dat$Slice)

Therm_Bin = as.data.frame(Therm_Dat %>% group_by(ID, New.ID, Time_Bin, Treatment, Sex) %>%
    summarise(Max.Eye = mean(Eye.Temp, na.rm = T),
    Max.Bill = mean(Bill.Temp, na.rm = T),
    Eye.HT = mean(Eye_HT, na.rm = T),
    Bill.HT = mean(Bill_HT, na.rm = T),
    m.qTot = mean(qTot, na.rm = T)))

Therm_Bin$oTreatment = factor(Therm_Bin$Treatment, ordered = T)

# Assessing sample sizes 

Therm_Bin %>% filter(!is.na(Max.Eye)) %>% 
group_by(Treatment, ID, Sex) %>% summarise(n()) %>%
group_by(Treatment, Sex) %>% summarise(n())

Therm_Bin %>% filter(!is.na(Max.Eye)) %>% group_by(Treatment) %>% summarise(n())
Therm_Bin %>% filter(!is.na(Max.Eye)) %>% group_by(Treatment, New.ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())
Therm_Bin %>% filter(!is.na(Max.Eye)) %>% group_by(ID, Treatment) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())

EyeTemp_Mod = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, data = Therm_Bin)
head(predict(EyeTemp_Mod, data = Therm_Bin, type = "lpmatrix"))

Therm_Bin$oTreatment = factor(Therm_Bin$oTreatment, ordered = FALSE)
AR_Eye = start_event(Therm_Bin, column="Time_Bin", event=c("New.ID"))
AR_Eye$oTreatment = factor(AR_Eye$oTreatment, ordered = TRUE)
Rho = acf(resid_gam(EyeTemp_Mod), plot = F)$acf[2]
Rho

EyeTemp_Mod_AR = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Eye$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = AR_Eye)

# Assessing model residuals

{
p1 = ggplot(data = data.frame("Res" = resid_gam(EyeTemp_Mod_AR)),
  aes(x = 1:length(resid_gam(EyeTemp_Mod_AR)), y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalised Residuals")

p2 = ggplot(data = data.frame("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE), 
  "Fit" = predict(EyeTemp_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)))),
  aes(x = Fit, y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalised Residuals")

p3 = ggplot(data.frame("Res" = resid_gam(EyeTemp_Mod_AR)), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(resid_gam(EyeTemp_Mod_AR)), 
    sd = sd(resid_gam(EyeTemp_Mod_AR)))) +
  theme_bw() + xlab("Normalised Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(resid_gam(EyeTemp_Mod_AR)) - 
    3*sd(resid_gam(EyeTemp_Mod_AR))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(resid_gam(EyeTemp_Mod_AR)) + 
    3*sd(resid_gam(EyeTemp_Mod_AR))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE), 
    "Fit" = predict(EyeTemp_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

p5 = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)) %>% 
  mutate("Pred" = predict(EyeTemp_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)))) %>%
    ggplot(aes(x = Pred, y = Max.Eye)) + 
    geom_point(size = 2, colour = "black", fill = "orchid", pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = FALSE) + 
    my.theme + theme_bw()

p6 = ACF_Plot(model = EyeTemp_Mod_AR, residual_type = "normalised")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, top = "Model Residuals")
}    

# Notable tails in qqplot, and predictive capacity of model appears quite low. 
# Plotting by individual to assess oddities. 

Group_Diag(model = EyeTemp_Mod_AR, 
    data = as.data.frame(Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)), group = "ID", view = "Time_Bin")

Group_Diag(model = EyeTemp_Mod_AR, 
    data = as.data.frame(Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)), group = "ID", label = "oTreatment", view = "Time_Bin")

# Raw responses and extreme residuals appear evenly dispersed. Checking whether
# extreme residuals fall within only one treatment, one bird, or one time bin.

p1 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = oTreatment, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p2 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Res, fill = oTreatment)) + 
    geom_density(alpha = 0.5, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p3 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Time_Bin, y = Res, fill = oTreatment, linetype = oTreatment)) + 
    stat_summary_bin(fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.5, 
    binwidth = 50) + stat_summary_bin(fun = "mean", geom = "point", size = 2, pch = 21,
        colour = "black", binwidth = 10) + 
    geom_smooth(method = "lm", size = 1, colour = "black") +    
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p4 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = ID, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4, position = position_dodge(width = 0.5)) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black",
        position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Residual Distributions")

# Handled means appear particularly deviant above 800 seconds, and TZ19 appears extremely 
# variable when compared with other individuals. Checking whether heteroskedasticity is supported by a levene test

Var_Test = Therm_Bin %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.) %>% mutate("Res" = resid_gam(EyeTemp_Mod_AR, incl_na = TRUE))
car::leveneTest(Res ~ ID, Var_Test)

ggplot(subset(Var_Test, ID != "TZ19"), aes(x = ID, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4, position = position_dodge(width = 0.5)) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black",
        position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw() # Much better.

# Removing TZ19 and reassessing model.

Therm_Bin_Rev = Therm_Bin %>% filter(ID != "TZ19") %>% mutate("ID" = factor(ID))

EyeTemp_Mod_Rev = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, data = Therm_Bin_Rev)

Therm_Bin_Rev$oTreatment = factor(Therm_Bin_Rev$oTreatment, ordered = FALSE)
AR_Eye_Rev = start_event(Therm_Bin_Rev, column="Time_Bin", event=c("New.ID"))
Therm_Bin_Rev$oTreatment = factor(Therm_Bin_Rev$oTreatment, ordered = TRUE)
AR_Eye_Rev$oTreatment = factor(AR_Eye_Rev$oTreatment, ordered = TRUE)

Rho_Rev = acf(resid_gam(EyeTemp_Mod_Rev), plot = F)$acf[2]
Rho_Rev

# Rho remains high. 

EyeTemp_Mod_Rev_AR = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

# Reassessing residuals.

{
p1 = ggplot(data = data.frame("Res" = resid_gam(EyeTemp_Mod_Rev_AR)),
  aes(x = 1:length(resid_gam(EyeTemp_Mod_Rev_AR)), y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalised Residuals")

p2 = ggplot(data = data.frame("Res" = resid_gam(EyeTemp_Mod_Rev_AR, incl_na = TRUE), 
  "Fit" = predict(EyeTemp_Mod_Rev_AR, newdata = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)))),
  aes(x = Fit, y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalised Residuals")

p3 = ggplot(data.frame("Res" = resid_gam(EyeTemp_Mod_Rev_AR)), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(resid_gam(EyeTemp_Mod_Rev_AR)), 
    sd = sd(resid_gam(EyeTemp_Mod_Rev_AR)))) +
  theme_bw() + xlab("Normalised Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(resid_gam(EyeTemp_Mod_Rev_AR)) - 
    3*sd(resid_gam(EyeTemp_Mod_Rev_AR))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(resid_gam(EyeTemp_Mod_Rev_AR)) + 
    3*sd(resid_gam(EyeTemp_Mod_Rev_AR))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = resid_gam(EyeTemp_Mod_Rev_AR, incl_na = TRUE), 
    "Fit" = predict(EyeTemp_Mod_Rev_AR, newdata = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

p5 = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)) %>% 
  mutate("Pred" = predict(EyeTemp_Mod_Rev_AR, newdata = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, oTreatment, Time_Bin, ID) %>% drop_na(.)))) %>%
    ggplot(aes(x = Pred, y = Max.Eye)) + 
    geom_point(size = 2, colour = "black", fill = "orchid", pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = FALSE) + 
    my.theme + theme_bw()

p6 = ACF_Plot(model = EyeTemp_Mod_Rev_AR, residual_type = "normalised")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, top = "Model Residuals")
}    

# Little improvement. Assessing residuals across predictors again.

p1 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_Rev_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = oTreatment, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p2 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_Rev_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Res, fill = oTreatment)) + 
    geom_density(alpha = 0.5, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p3 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_Rev_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Time_Bin, y = Res, fill = oTreatment, linetype = oTreatment)) + 
    geom_point(size = 2, pch = 21, colour = "black") + 
    geom_smooth(method = "lm", size = 1, colour = "black") +    
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p4 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(EyeTemp_Mod_Rev_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = ID, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4, position = position_dodge(width = 0.5)) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black",
        position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Residual Distributions")

# Assessing whether time-splines carry enough curvature.

gam.check(EyeTemp_Mod_Rev_AR, type = "pearson")

# Both responses are largely linear. Testing a glm approach using 
# nlme.

require('nlme')

lin_Eyemod = lme(Max.Eye ~ Treatment*Time_Bin, random = ~1|ID,
    correlation = corAR1(form = ~1|ID),
    na.action = na.omit, 
    data = Therm_Bin_Rev)

acf(residuals(lin_Eyemod, type = "normalized"))

# Autocorrelation appears clear. Plotting residuals.

{
p1 = ggplot(data = data.frame("Res" = residuals(lin_Eyemod, type = "normalized")),
  aes(x = 1:length(residuals(lin_Eyemod, type = "normalized")), y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalised Residuals")

p2 = ggplot(data = data.frame("Res" = residuals(lin_Eyemod, type = "normalized"), 
  "Fit" = predict(lin_Eyemod, newdata = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, Treatment, Time_Bin, ID) %>% drop_na(.)))),
  aes(x = Fit, y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalised Residuals")

p3 = ggplot(data.frame("Res" = residuals(lin_Eyemod, type = "normalized")), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(residuals(lin_Eyemod, type = "normalized")), 
    sd = sd(residuals(lin_Eyemod, type = "normalized")))) +
  theme_bw() + xlab("Normalised Residuals") + ylab("Count") +
  geom_vline(xintercept = mean(residuals(lin_Eyemod, type = "normalized") - 
    3*sd(residuals(lin_Eyemod, type = "normalized"))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(residuals(lin_Eyemod, type = "normalized")) + 
    3*sd(residuals(lin_Eyemod, type = "normalized"))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = residuals(lin_Eyemod, type = "normalized"), 
    "Fit" = predict(lin_Eyemod, newdata = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, Treatment, Time_Bin, ID) %>% drop_na(.)))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

p5 = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, Treatment, Time_Bin, ID) %>% drop_na(.)) %>% 
  mutate("Pred" = predict(lin_Eyemod, newdata = as.data.frame(Therm_Bin_Rev %>% 
  dplyr::select(Max.Eye, Treatment, Time_Bin, ID) %>% drop_na(.)))) %>%
    ggplot(aes(x = Pred, y = Max.Eye)) + 
    geom_point(size = 2, colour = "black", fill = "orchid", pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = FALSE) + 
    my.theme + theme_bw()

p6 = ACF_Plot(model = lin_Eyemod, residual_type = "normalised")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, top = "Model Residuals")
}    

# By predictors.

p1 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, Treatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = residuals(lin_Eyemod, type = "normalized")) %>% 
        ggplot(aes(x = Treatment, y = Res, fill = Treatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p2 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, Treatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = residuals(lin_Eyemod, type = "normalized")) %>% 
        ggplot(aes(x = Res, fill = Treatment)) + 
    geom_density(alpha = 0.5, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p3 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, Treatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = residuals(lin_Eyemod, type = "normalized")) %>% 
        ggplot(aes(x = Time_Bin, y = Res, fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black") + 
    geom_smooth(method = "lm", size = 1, colour = "black") +    
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p4 = as.data.frame(Therm_Bin_Rev %>% dplyr::select(Max.Eye, Time_Bin, Treatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = residuals(lin_Eyemod, type = "normalized")) %>% 
        ggplot(aes(x = ID, y = Res, fill = Treatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4, position = position_dodge(width = 0.5)) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black",
        position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Residual Distributions")

# Again, few points deviating; probably holding little bearing on model. Assessing Cook's distance.

require('predictmeans')

mod_inf = CookD(model = lin_Eyemod, group = NULL)
Therm_Bin_Rev$CooksD = mod_inf

ggplot(Therm_Bin_Rev, aes(x = Time_Bin, y = CooksD, fill = Treatment, linetype = Treatment)) + 
    geom_point(pch = 21, size = 2, colour = "black") + 
    geom_smooth(method = "lm", size = 1, colour = "black") + 
    theme_bw()

# U-shape, suggesting that a non-linear model is probably more appropriate. Testing gam with penalisation removed. 

EyeTemp_Mod_Rev_AR = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

EyeTemp_Mod_Rev_AR_NP = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3, fx = TRUE) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1, fx = TRUE) + 
    s(ID, bs = "re"),
    AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

summary(EyeTemp_Mod_Rev_AR)
summary(EyeTemp_Mod_Rev_AR_NP)

# Significance levels do not change across models. Assessing influence of non-penalised model.

mod_fits = fitted(EyeTemp_Mod_Rev_AR_NP, type = "response")

leave_out = function(x) {
  updated_data = AR_Eye_Rev[-x, ]
  model = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3, fx = TRUE) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1, fx = TRUE) + 
    s(ID, bs = "re"),
    AR.start = updated_data$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = updated_data)

  sum((mod_fits[-x] - fitted(model))^2)
}

influence = sapply(1:nrow(AR_Eye_Rev), leave_out)
AR_Eye_Rev$Influence = influence

Crit_I = AR_Eye_Rev %>% filter(!is.na(Max.Eye) & start.event == "FALSE") %>%
    summarise(Crit = 3*mean(Influence, na.rm = T))

ggplot(as.data.frame(AR_Eye_Rev %>% filter(!is.na(Max.Eye) & start.event == "FALSE")), 
    aes(x = Time_Bin, y = Influence, fill = oTreatment, linetype = oTreatment)) + 
    geom_point(pch = 21, size = 2, colour = "black") + 
    geom_smooth(method = "lm", size = 1, colour = "black") + 
    theme_bw() + geom_hline(yintercept = Crit_I$Crit, linetype = "longdash",
        size = 1, colour = "mediumvioletred") + 
    scale_fill_manual(values = c("black", "grey70"))

# Indentifying those with high influence

AR_Eye_Rev %>% filter(!is.na(Max.Eye) & start.event == "FALSE") %>%
    filter(Influence > Crit_I$Crit)

# Interesting. Largely individual TZ24. Plotting raw reponses per individual.

ggplot(AR_Eye_Rev, aes(x = Time_Bin, y = Max.Eye, colour = oTreatment)) + 
    geom_path(size = 1, aes(group = New.ID)) + facet_wrap(~ID) + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw()

# No apparent reason to remove TZ24 in either treatment. Proceeding with model as is.
# Assessing deviance explained.

Null_Mod = bam(Max.Eye ~ 1 + 
  s(ID, bs = "re"),
  AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

EyeTemp_Mod_Rev_AR = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

(sum(as.numeric(residuals(Null_Mod, type = "deviance"))^2) - 
sum(as.numeric(residuals(EyeTemp_Mod_Rev_AR, type = "deviance"))^2))/
(sum(as.numeric(residuals(Null_Mod, type = "deviance"))^2))

# <10% deviance explained. Calculating R-squared. 

summary(lm(AR_Eye_Rev$Max.Eye~predict(EyeTemp_Mod_Rev_AR, type = "response", newdata = AR_Eye_Rev)))

# ~ 58.6. Acceptable model, but presumably eye temperature is not well predicted by time or treatment.
# Comparing coefficients by group.

Coef_Dist(model = EyeTemp_Mod_Rev_AR, incl_int = FALSE)

# Slight negative trend in second time interval in controls, but non-significant.
# Lastly, assessing uncertainty around coefficients.

post_simGam(EyeTemp_Mod_Rev_AR, AR_Eye_Rev, view = "Time_Bin", group_var = "oTreatment", sim = 50)

# Considerable uncertainty and overlap between groups. Plotting marginal means.

Eye_MMDat = emmip(EyeTemp_Mod_Rev_AR, oTreatment ~ Time_Bin, data = AR_Eye_Rev,
    cov.reduce = FALSE, nesting = NULL, nesting.order = FALSE,
    CIs = TRUE, plot = F)

ggplot(Eye_MMDat, aes(x = Time_Bin, y = yvar, 
    linetype = oTreatment, fill = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cr"), size = 1, se = F) +
    stat_summary_bin(data = AR_Eye_Rev, aes(x = Time_Bin, y = Max.Eye),
        geom = "point", fun = "mean", binwidth = 10, pch = 21, size = 2) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    theme_bw() + xlab("Frame Number") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "black", size = 1, linetype = "dashed")

# One considerable outlier. Adjusting y-limits to ignore this point for now

Eye_plot = ggplot(Eye_MMDat, aes(x = Time_Bin, y = yvar, 
    linetype = oTreatment, fill = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cr"), size = 1, se = F) +
    stat_summary_bin(data = AR_Eye_Rev, aes(x = Time_Bin, y = Max.Eye),
        geom = "point", fun = "mean", binwidth = 10, pch = 21, size = 2) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "black", size = 1, linetype = "dashed") + 
    ylim(c(33.5, 36)) + my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210"))

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Eye Temperature Response.jpeg", 
    Eye_plot, dpi = 800, width = 7.67, height = 6.67, units = "in")

# Computing simultaneous confidence intervals around coefficient estimates. This method
# follows Ruppert et al. (2003), as transcribed by G. Simpson 
# (url = https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/)
# with slight alterations to expand on random variables.

rmvn_custom = function(n, mu, sig) {
    L = mroot(sig)
    m = ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}

# Extracting model covariance matrix and generating predictions for an expanded grid
# of predictors. Note that the covariance matrix has been estimated as conditional
# upon smoothing parameters. An expanded grid is used here to average uncertainty
# across influence of random effect. Importantly, variance around random intercepts
# has been set to zero for simulation purposes; intercepts are therefor fixed to 
# a convergence at zero. 

var_cov = vcov(EyeTemp_Mod_Rev_AR)
diag(var_cov)[7:15] = 0

# Expansion

pred_dat = with(Therm_Bin_Rev, expand.grid("ID" = unique(ID), 
    "Time_Bin" = seq(0, 840, by = 10), "oTreatment" = unique(oTreatment)))

pred_vals = predict(EyeTemp_Mod_Rev_AR, newdata = pred_dat, se.fit = TRUE)
Cg = predict(EyeTemp_Mod_Rev_AR, pred_dat, type = "lpmatrix")

fitted_se = pred_vals$se.fit

# Simulating maximal absolute deviation from model estimates 
# (using standard errors).

set.seed(20)
N_sim = 10000

# Setting mean of simulated multivariate normal distribution to zero, 
# with variance estimated from our model.

BUdiff = rmvn_custom(N_sim, mu = rep(0, nrow(var_cov)), sig = var_cov)

# Calcuting deviation from predicted model across simulations, then calculating
# absolute deviation from true model using predicted standard errors.

simDeviation = Cg %*% t(BUdiff)

# Absolute deviation is estimated by simulated deviation/fitted standard error values

absDeviation = abs(sweep(simDeviation, 1, fitted_se, FUN = "/"))

# Maximal absolute deviation is pulled from each simulation-model comparison, then 
# critical value estimated from maximal deviations (95% confidence is assumed here).
# Note that because maximum is truly supremum, max values are calculated for each
# categorical subset.

MA_Cnt = apply(absDeviation[c(which(pred_dat$oTreatment == "Control")),], 2L, max)
MA_Hand = apply(absDeviation[c(which(pred_dat$oTreatment == "Handled")),], 2L, max)

# Assessing deviation of maximal values by treatment

par(mfrow=c(1,2))
hist(MA_Cnt) 
hist(MA_Hand) 
par(mfrow=c(1,1))

# Deviance values are right-skewed, so using type 8 quantiles.

CV_Cnt = quantile(MA_Cnt, prob = 0.95, type = 8)
CV_Hand = quantile(MA_Hand, prob = 0.95, type = 8)

# Creating data frame to plot confidence intervals, then averaging at time intervals
# by treatment.

pred_df = transform(cbind(data.frame(pred_vals), pred_dat),
                  uprS = fit + (CV_Cnt * fitted_se),
                  lwrS = fit - (CV_Cnt * fitted_se))

pred_df$uprS[c(which(pred_dat$oTreatment == "Handled"))] = 
    pred_vals$fit[c(which(pred_dat$oTreatment == "Handled"))] + 
    (CV_Hand*fitted_se[c(which(pred_dat$oTreatment == "Handled"))])

pred_df$lwrS[c(which(pred_dat$oTreatment == "Handled"))] = 
    pred_vals$fit[c(which(pred_dat$oTreatment == "Handled"))] - 
    (CV_Hand*fitted_se[c(which(pred_dat$oTreatment == "Handled"))])

pred_df_grouped = pred_df %>% group_by(Time_Bin, oTreatment) %>%
    summarise(fit = mean(fit), se.fit = mean(se.fit), 
    UCL = mean(uprS), LCL = mean(lwrS))

# Plotting 

ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw()

# Simulating response again to overlay.

Eye_Coef = EyeTemp_Mod_Rev_AR$coefficients

set.seed(20)

# Now simulating multivariate normal distribution around coefficients

var_cov = vcov(EyeTemp_Mod_Rev_AR)
sim_gam = rmvn_custom(1000, mu = Eye_Coef, sig = var_cov)
mod_sims_fit = Cg %*% t(sim_gam)
ylims = range(mod_sims_fit)

# Setting up data for shading per group

Pred_Long_Form = vector('list', ncol(mod_sims_fit))

for (i in 1:length(Pred_Long_Form)){
    Pred_Long_Form[[i]] = pred_dat %>% mutate("Pred" = mod_sims_fit[,i], "Sim_Number" = i)
}

Pred_Long_Form = bind_rows(Pred_Long_Form)
Pred_Long_Form$Group = paste(Pred_Long_Form[,"oTreatment"], Pred_Long_Form$Sim_Number, sep = "_")

# Binding with confidence intervals in plot

ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    geom_path(data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin, Group) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = Group),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("mediumvioletred", "mediumvioletred"), name = "Treatment") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw()

# CIs appear suitable. Plotting with mean trend-lines and raw data points

ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.4) +
    stat_summary(fun = "mean", geom = "smooth", data = Pred_Long_Form,
        aes(x = Time_Bin, y = Pred, linetype = oTreatment, colour = oTreatment),
        size = 1) +  
    geom_point(data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye),
        pch = 21, size = 2) +  
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw() + xlab("Frame Number") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, size = 1, linetype = "longdash", colour = "black") + 
    my.theme + facet_wrap(~ID)

# Significant variation in points around trends.

# Comparing output with itsadug

ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    scale_linetype_manual(values = c("solid", "dashed"), name = "Treatment") +
    theme_bw()

dev.new() 

plot_smooth(EyeTemp_Mod_Rev_AR, view = "Time_Bin", plot_all = "oTreatment", rm.ranef = FALSE, 
    sim.ci = TRUE)

get_predictions(EyeTemp_Mod_Rev_AR, cond = list("Time_Bin" = seq(0, 840, 10), "oTreatment" = c("Control", "Handled")), 
    rm.ranef = FALSE, sim.ci = TRUE)
    
# Itsadug intervals are similar, but set for the values of one individual alone.
# Saving final plots with custom estimated simultaneous confidence intervals.

Eye_plot_correct_CI = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Eye Correct CIs.jpeg", 
    Eye_plot_correct_CI, dpi = 800, width = 7.76, height = 6.67, units = "in")

Therm_Dat$oTreatment = Therm_Dat$OTreatment
Therm_Dat$oTreatment = Therm_Dat$Treatment
P_Dat = as.data.frame(Therm_Dat %>% filter(ID != "TZ19"))

Eye_plot_with_points = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye),
        binwidth = 10, pch = 21, size = 2, alpha = 0.6) + 
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) +     
    annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 38.75) +  
    #stat_summary_bin(data = as.data.frame(Therm_Dat %>% filter(ID != "TZ19")), 
    #    binwidth = 0.5, geom = "point",
    #    fun = "mean", aes(x = Slice, y = Eye.Temp, fill = oTreatment),
    #    pch = 21) +     
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Eye Correct CIs with points per second - Colour Matched and Legend on Plot.jpeg", Eye_plot_with_points, dpi = 800, width = 7, height = 6.67, units = "in")

# And again, with confidence intervals around points

Eye_plot_with_points_and_ses = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) +     
    annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 38.75) +  
    #stat_summary_bin(data = as.data.frame(Therm_Dat %>% filter(ID != "TZ19")), 
    #    binwidth = 0.5, geom = "point",
    #    fun = "mean", aes(x = Slice, y = Eye.Temp, fill = oTreatment),
    #    pch = 21) +     
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_plot_with_points_and_ses

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Eye Correct CIs with points per second - Colour Matched and Legend on Plot - CIs added.jpeg", Eye_plot_with_points_and_ses, dpi = 800, width = 7, height = 6.67, units = "in")

# Lastly, with grey rectangles representing duration of handling

Eye_plot_with_points_and_ses = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) +     
    annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 790, y = 38.75) +  
    #stat_summary_bin(data = as.data.frame(Therm_Dat %>% filter(ID != "TZ19")), 
    #    binwidth = 0.5, geom = "point",
    #    fun = "mean", aes(x = Slice, y = Eye.Temp, fill = oTreatment),
    #    pch = 21) +     
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Eye Region Temperature (°C)") + 
    annotate("rect", xmin = 420, xmax = 840, ymin = -Inf, ymax = Inf, colour = "black", fill = "grey40", alpha = 0.4) + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + 
    theme(legend.position = c(0.21,0.1), legend.background = NULL)    

Eye_plot_with_points_and_ses

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Eye Correct CIs with points per second - Colour Matched and Legend on Plot - CIs and Grey Bands added.jpeg", Eye_plot_with_points_and_ses, dpi = 800, width = 7, height = 6.67, units = "in")

# Summarising model

summary(EyeTemp_Mod_Rev_AR)

# Pulling out standard errors and mean coefficients, along with marginal means per time group

Coefs = EyeTemp_Mod_Rev_AR$coefficients
Coef.SE = as.data.frame(sqrt(diag(vcov(EyeTemp_Mod_Rev_AR, unconditional = TRUE))))
Coef.SE$Row = c(1:nrow(Coef.SE))
Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Coef.SE))))

data.frame("Coefficient" = Coef_Names, "Estimate" = Coefs, "SE" = Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

Eye_Grid = ref_grid(EyeTemp_Mod_Rev_AR,  
    at = list(Time_Bin = seq(420, 840, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

Eye_Grid_Before = ref_grid(EyeTemp_Mod_Rev_AR,  
    at = list(Time_Bin = seq(0, 419, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

emmeans(Eye_Grid, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

emmeans(Eye_Grid_Before, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

# Assessing baseline temperatures (here, within the first 30 seconds and minute of observation)

Eye_Grid_Base_30 = ref_grid(EyeTemp_Mod_Rev_AR,  
    at = list(Time_Bin = seq(0, 30, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

Eye_Grid_Base_60 = ref_grid(EyeTemp_Mod_Rev_AR,  
    at = list(Time_Bin = seq(0, 60, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

emmeans(Eye_Grid_Base_30, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

emmeans(Eye_Grid_Base_60, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

# For curiosity sake, calculating first order derivatives of eye responses across time

pred_dat = with(Therm_Bin_Rev, expand.grid("ID" = unique(ID), 
    "Time_Bin" = seq(0, 840, by = 10), "oTreatment" = unique(oTreatment)))

X0 = predict(EyeTemp_Mod_Rev_AR, pred_dat, type = "lpmatrix")

increment = 1*10^-6
pred_dat_adj = pred_dat %>% mutate("Time_Bin" = Time_Bin + increment)
coefficients = EyeTemp_Mod_Rev_AR$coefficients

var_covar = vcov(EyeTemp_Mod_Rev_AR)
diag(var_covar)[7:15] = 0
sim_gam = rmvn_custom(10000, mu = coefficients, sig = var_covar)

X1 = predict(EyeTemp_Mod_Rev_AR, pred_dat_adj, type = "lpmatrix")

Pred_X0 = X0 %*% t(BUdiff)
Pred_X1 = X1 %*% t(BUdiff)

Diff_Mat = Pred_X0*0

for (i in 1:ncol(Diff_Mat)){
    Diff_Mat[,i] = Pred_X1[,i] - Pred_X0[,i]
}

Diff_Mat_Corrected = sweep(Diff_Mat, 1, (60*10^6), FUN = "*")

Long_derivs_Aft = vector('list', ncol(Diff_Mat_Corrected))
for (i in 1:length(Long_derivs_Aft)){
     Long_derivs_Aft[[i]] = pred_dat %>% mutate("Pred" = Diff_Mat_Corrected[,i], "Sim_Number" = i)
}

Long_derivs_Aft = bind_rows(Long_derivs_Aft)

Long_derivs_Aft_Av = Long_derivs_Aft %>% group_by(Time_Bin, oTreatment) %>%
    summarise("Deriv" = mean(Pred), LCL = quantile(Pred, probs = c(0.025, 0.975), type = 8)[1],
    UCL = quantile(Pred, probs = c(0.025, 0.975), type = 8)[2])

Long_derivs_Aft_Av$Increasing = Long_derivs_Aft_Av$Deriv
Long_derivs_Aft_Av$Decreasing = Long_derivs_Aft_Av$Deriv
Long_derivs_Aft_Av$Increasing[c(which(Long_derivs_Aft_Av$Deriv < 0))] = NA
Long_derivs_Aft_Av$Decreasing[c(which(Long_derivs_Aft_Av$Deriv >= 0))] = NA

Eye_Derivs = ggplot(Long_derivs_Aft_Av, aes(x = Time_Bin, y = Deriv, 
    group = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) +
    geom_line(aes(y = Decreasing), size = 2.0, alpha = 0.5) + 
    geom_line(size = 0.5) + 
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature Slope (°C/min)") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    geom_vline(xintercept = 420, linetype = "longdash", colour = "black") +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-420", "-210", "0", "210", "420"))

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Eye Derivatives.jpeg", 
    Eye_Derivs, dpi = 800, width = 7.76, height = 6.67, units = "in")

# Finally, quantifying whether inclusion of a random linear slope of time per 
# individual improves log-likelihood of the model

EyeTemp_Base = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

summary(EyeTemp_Base)
AR_Eye_Rev = AR_Eye_Rev %>% 
    mutate(New.ID = factor(New.ID))

EyeTemp_RandomSlope = bam(Max.Eye ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re") + 
    s(Time_Bin, ID, oTreatment, bs = "re"),
    AR.start = AR_Eye_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_Eye_Rev)

test_stat = -2 * (logLik(EyeTemp_Base) - logLik(EyeTemp_RandomSlope))
pchisq(as.numeric(test_stat), df = (attr(logLik(EyeTemp_RandomSlope), "df") - 
    attr(logLik(EyeTemp_Base), "df")), lower.tail = FALSE)

# Plotting individual responses to treatments 

Exp_Eye = expand.grid("Time_Bin" = seq(1, 840, by = 1),
    "ID" = unique(na.omit(AR_Eye_Rev$ID)),
    "oTreatment" = unique(na.omit(AR_Eye_Rev$oTreatment))) %>% 
    mutate("Group" = paste(oTreatment, ID, sep = "_"))

# Removing combinations that did not occur

Compare = AR_Eye_Rev %>% 
    mutate("Group" = paste(oTreatment, ID, sep = "_"))

Catch = unique(Exp_Eye$Group)[c(which(!(unique(Exp_Eye$Group) %in% unique(Compare$Group))))]
Exp_Eye = Exp_Eye[-c(which(Exp_Eye$Group %in% Catch)),]

# Predicting

Response_Predictions = predict(EyeTemp_RandomSlope, newdata = Exp_Eye, exclude = "s(ID)", se.fit = TRUE)

Eye_Ind = Exp_Eye %>%
    mutate("Fit" = Response_Predictions$fit,
        "LCL" = Response_Predictions$fit - 1.96*Response_Predictions$se.fit,
        "UCL" = Response_Predictions$fit + 1.96*Response_Predictions$se.fit) %>% 
    mutate("Groups" = paste(oTreatment, ID, sep = "_")) %>% 
    mutate(Groups = factor(Groups))

Eye_Individual_Slopes = ggplot(Eye_Ind, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, group = Groups), alpha = 0.15) +
    #stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 10, pch = 21, size = 2, alpha = 0.6) + 
    geom_line(aes(x = Time_Bin, y = Fit, group = Groups), size = 1, alpha = 0.7,
        colour = "black") +     
    annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 38.75) +    
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_Individual_Slopes

extractLegend <- function(xplot) {
  grobs <- ggplot_gtable(ggplot_build(xplot))
  g_title <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[g_title]]
}

Legend_Pull = ggplot(Eye_Ind, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, group = Groups), alpha = 0.75) +
    #stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 10, pch = 21, size = 2, alpha = 0.6) + 
    geom_line(aes(x = Time_Bin, y = Fit, group = Groups), size = 1, alpha = 0.4,
        colour = "black") +     
    annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 38.75) +    
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

NewLeg = extractLegend(Legend_Pull)
replace_Grob = ggplot_gtable(ggplot_build(Eye_Individual_Slopes))
rep_new = which(sapply(replace_Grob$grobs, function(x) x$name) == "guide-box")
replace_Grob$grobs[[rep_new]] = NewLeg

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/Individual_Responses_At_Eye.jpeg", replace_Grob, height = 7, width = 8, dpi = 800
)

###### Running model with change from baseline by individual

TBin_Treat = split(Therm_Bin, f = Therm_Bin$Treatment)
for (i in 1:2){
    TBin_Treat[[i]] = split(TBin_Treat[[i]], f = TBin_Treat[[i]]$New.ID)
}  

for (j in 1:2){
    for (i in 1:length(TBin_Treat[[j]])){
        Base_Eye = TBin_Treat[[j]][[i]]$Max.Eye[which(TBin_Treat[[j]][[i]]$Time_Bin == min(TBin_Treat[[j]][[i]]$Time_Bin, na.rm = T))]
        TBin_Treat[[j]][[i]]$Eye_Diff = TBin_Treat[[j]][[i]]$Max.Eye - Base_Eye
    }
}    

for (i in 1:2){
    TBin_Treat[[i]] = bind_rows(TBin_Treat[[i]])
}

Eye_Diff = bind_rows(TBin_Treat)

# Modeling 

Eye_Diff_Mod = bam(Eye_Diff ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = Eye_Diff)

Eye_Diff$oTreatment = factor(Eye_Diff$oTreatment, ordered = FALSE)
Eye_Diff_AR = start_event(Eye_Diff, column="Time_Bin", event=c("New.ID"))
Eye_Diff_AR$oTreatment = factor(Eye_Diff_AR$oTreatment, ordered = TRUE)
Eye_Diff$oTreatment = factor(Eye_Diff$oTreatment, ordered = TRUE)

Eye_Diff_Rho = acf(resid_gam(Eye_Diff_Mod), plot = F)$acf[2]
Eye_Diff_Rho

Eye_Diff_Mod_AR = bam(Eye_Diff ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = Eye_Diff_AR$start.event, rho = Eye_Diff_Rho,
    method = "REML", na.action = na.omit, 
    data = Eye_Diff_AR)

emmip(Eye_Diff_Mod_AR, oTreatment ~ Time_Bin, CIs = TRUE, 
    at = list(Time_Bin = seq(1, 840, 1)),
    data = Eye_Diff_AR)

summary(Eye_Diff_Mod_AR)

# Interesting. Again, effect of treatment on trend not quite significant, although nearing so. Deviance explained is considerably lower as well (oddly). Plotting.

EM_Data_Diff = emmip(Eye_Diff_Mod_AR, oTreatment ~ Time_Bin, CIs = TRUE, 
    at = list(Time_Bin = seq(1, 840, 1)),
    data = Eye_Diff_AR, plot = FALSE)

ggplot(EM_Data_Diff, aes(x = Time_Bin, y = yvar, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.55, colour = NA) + 
    geom_line(size = 1, colour = "black") + 
    stat_summary_bin(data = Eye_Diff, aes(x = Time_Bin, y = Eye_Diff), geom = "errorbar", size = 1, width = 1, colour = "black",fun.data = "mean_cl_boot", binwidth = 50, position = position_dodge(width = 25)) + 
    stat_summary_bin(data = Eye_Diff, aes(x = Time_Bin, y = Eye_Diff), geom = "point", size = 3, pch = 21, colour = "black",fun = "mean", binwidth = 50, position = position_dodge(width = 25), alpha = 0.8) + theme_bw() +    
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + xlab("Time (s)") + ylab("Difference in Eye Temperature from Baseline (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme +  
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

# Results look nice. Perhaps unaveraged data is considerably noisy?

# Lastly, comparing results with those of model including TZ19 just in case

summary(EyeTemp_Mod_AR)
emmip(EyeTemp_Mod_AR, oTreatment ~ Time_Bin, data = Therm_Bin,
    CIs = TRUE, at = list(Time_Bin = seq(1, 840, by = 1)))
Coefs = EyeTemp_Mod_AR$coefficients
Coef.SE = as.data.frame(sqrt(diag(vcov(EyeTemp_Mod_AR, unconditional = TRUE))))
Coef.SE$Row = c(1:nrow(Coef.SE))
Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Coef.SE))))

data.frame("Coefficient" = Coef_Names, "Estimate" = Coefs, "SE" = Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

# Yes, results appear strikingly similar with inclusion of this individual.

#------------------------------------------------------------------------------#
## Modeling bill temperature responses ##
#------------------------------------------------------------------------------#

# Assessing sample sizes 

Therm_Bin %>% filter(!is.na(Max.Bill)) %>% group_by(Treatment) %>% summarise(n())
Therm_Bin %>% filter(!is.na(Max.Bill)) %>% group_by(Treatment, New.ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())
Therm_Bin %>% filter(!is.na(Max.Bill)) %>% group_by(Treatment, ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())

# Quickly assessing rho without time averages 

BillTemp_Mod_Long = gam(Bill.Temp ~ OTreatment +
    s(Slice, bs = "cr", k = 3) + 
    s(Slice, by = OTreatment, bs = "cr", k = 3, m = 1) + 
    s(ID, bs = "re"),  
    method = "REML", na.action = na.omit, data = Therm_Dat)

ACF_Plot(model = BillTemp_Mod_Long, residual_type = "normalised")
acf(residuals(BillTemp_Mod_Long), plot = FALSE)$acf[2]

# Continuing with averaged models

BillTemp_Mod = gam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),  
    method = "REML", na.action = na.omit, data = Therm_Bin)

ACF_Plot(model = BillTemp_Mod, residual_type = "normalised")

# Horrendous autocorrelation. Correcting with an AR1 structure in bam.

BillTemp_Mod = bam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),  
    method = "REML", na.action = na.omit, data = Therm_Bin)

Therm_Bin$oTreatment = factor(Therm_Bin$oTreatment, ordered = FALSE)
AR_Bill = start_event(Therm_Bin, column="Time_Bin", event=c("New.ID"))
Therm_Bin$oTreatment = factor(Therm_Bin$oTreatment, ordered = TRUE)
AR_Bill$oTreatment = factor(AR_Bill$oTreatment, ordered = TRUE)
Rho = acf(resid_gam(BillTemp_Mod), plot = F)$acf[2]
Rho

BillTemp_Mod_AR = bam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Bill$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = AR_Bill
)

ACF_Plot(model = BillTemp_Mod_AR, residual_type = "normalised")

# Much cleaner. Again, plotting model residuals. 

{
p1 = ggplot(data = data.frame("Res" = resid_gam(BillTemp_Mod_AR)),
  aes(x = 1:length(resid_gam(BillTemp_Mod_AR)), y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalised Residuals")

p2 = ggplot(data = data.frame("Res" = resid_gam(BillTemp_Mod_AR, incl_na = TRUE), 
  "Fit" = predict(BillTemp_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Bill, oTreatment, Time_Bin, ID) %>% drop_na(.)))),
  aes(x = Fit, y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalised Residuals")

p3 = ggplot(data.frame("Res" = resid_gam(BillTemp_Mod_AR)), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(resid_gam(BillTemp_Mod_AR)), 
    sd = sd(resid_gam(BillTemp_Mod_AR)))) +
  theme_bw() + xlab("Normalised Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(resid_gam(BillTemp_Mod_AR)) - 
    3*sd(resid_gam(BillTemp_Mod_AR))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(resid_gam(BillTemp_Mod_AR)) + 
    3*sd(resid_gam(BillTemp_Mod_AR))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = resid_gam(BillTemp_Mod_AR, incl_na = TRUE)), 
    "Fit" = predict(BillTemp_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Bill, oTreatment, Time_Bin, ID) %>% drop_na(.))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

p5 = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Bill, oTreatment, Time_Bin, ID) %>% drop_na(.)) %>% 
  mutate("Pred" = predict(BillTemp_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(Max.Bill, oTreatment, Time_Bin, ID) %>% drop_na(.)))) %>%
    ggplot(aes(x = Pred, y = Max.Bill)) + 
    geom_point(size = 2, colour = "black", fill = "orchid", pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = FALSE) + 
    my.theme + theme_bw()

p6 = ACF_Plot(model = BillTemp_Mod_AR, residual_type = "normalised")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, top = "Model Residuals")
}    

# Again, a few considerible outliers, however, cleaner than eye analyses. 
# Plotting individual slopes to identify.

Group_Diag(model = BillTemp_Mod_AR, cutoff = 5,
    data = as.data.frame(Therm_Bin %>% dplyr::select(Max.Bill, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)), group = "ID", view = "Time_Bin")

Group_Diag(model = BillTemp_Mod_AR, cutoff = 5,
    data = as.data.frame(Therm_Bin %>% dplyr::select(Max.Bill, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)), label = "oTreatment", group = "ID", view = "Time_Bin")

# Most outliers appears to fall at the beginning of a time-series.

# Assessing coefficients and residuals by predictors.

Coef_Dist(model = BillTemp_Mod_AR)

p1 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Bill, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(BillTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = oTreatment, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw() + theme(legend.position = "none")

p2 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Bill, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(BillTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Res, fill = oTreatment)) + 
    geom_density(alpha = 0.5, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw() + theme(legend.position = "none")

p3 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Bill, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(BillTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Time_Bin, y = Res, fill = oTreatment, 
        colour = oTreatment, linetype = oTreatment)) + 
    geom_point(pch = 21, size = 2, colour = "black") + 
    stat_summary_bin(fun.data = "mean_cl_boot", geom = "ribbon", binwidth = 20, alpha = 0.5) + 
    geom_smooth(method = "loess", size = 1) + 
    scale_fill_manual(values = c("black", "grey70")) +    
    scale_colour_manual(values = c("black", "grey70")) +                      
    theme_bw() + theme(legend.position = "none")

p4 = as.data.frame(Therm_Bin %>% dplyr::select(Max.Bill, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(BillTemp_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = ID, y = Res, fill = oTreatment)) + 
    geom_boxplot(alpha = 0.5, colour = "black", position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

lay = rbind(c(1,2,4,4,4),
            c(3,3,4,4,4),
            c(3,3,4,4,4))    

grid.arrange(p1, p2, p3, p4, layout_matrix = lay, top = "Residual Distributions")

# Some tailing in residuals, and mean is falling slightly lower among handled birds, but not considerably so.
# TZ19 appears far less variable here than in eye temperature modeling. No distict reason for concern. 

# Assessing influence, just in case.

BillTemp_Mod_AR_NP = bam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3, fx = TRUE) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1, fx = TRUE) + 
    s(ID, bs = "re"),
    AR.start = AR_Bill$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = AR_Bill
)

mod_fits = fitted(BillTemp_Mod_AR_NP, type = "response")

leave_out_AR = function(x) {

  if(AR_Bill[-x, "start.event"] == TRUE){
      return(NA)
  } else {
  updated_data = AR_Bill[-x, ]

  model = BillTemp_Mod_AR_NP = bam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3, fx = TRUE) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1, fx = TRUE) + 
    s(ID, bs = "re"),
    AR.start = updated_data$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = updated_data)

  return(sum((mod_fits[-x] - fitted(model))^2))
  }
}

influence = sapply(1:nrow(AR_Bill), leave_out)
AR_Bill$Influence = influence

hist(AR_Bill$Influence)

# Indentifying those with high influence again

Crit_I = 3*mean(AR_Bill$Influence, na.rm = TRUE)
Crit_I_Weak = mean(AR_Bill$Influence, na.rm = TRUE) + 
    3*sd(AR_Bill$Influence, na.rm = TRUE)

AR_Bill %>% filter(!is.na(Max.Bill)) %>%
    filter(Influence > Crit_I)    

AR_Bill %>% filter(!is.na(Max.Bill)) %>%
    filter(Influence > Crit_I_Weak)    
    
# Non with particularly high influence. Plotting raw reponses per individual.

ggplot(AR_Bill, aes(x = Time_Bin, y = Influence, fill = oTreatment, linetype = oTreatment)) + 
    geom_point(size = 2, pch = 21, colour = "black") + 
    geom_smooth(method = "lm", colour = "black", size = 1) + 
    geom_hline(yintercept = Crit_I) + 
    theme_bw() 

ggplot(AR_Bill, aes(x = Time_Bin, y = Influence, fill = oTreatment, linetype = oTreatment)) + 
    geom_point(size = 2, pch = 21, colour = "black") + 
    geom_smooth(method = "lm", colour = "black", size = 1) + 
    geom_hline(yintercept = Crit_I_Weak) + 
    theme_bw() 

# Three observations at weak influence cut off, with little biological purpose for their 
# exclusion. Proceeding with model as is.

# Plotting marginal means 

Bill_MMDat = emmip(BillTemp_Mod_AR, oTreatment ~ Time_Bin, cov.reduce = F, 
    data = Therm_Bin, nesting = NULL, nesting.order = FALSE, CIs = T, plot = F)

ggplot(Bill_MMDat, aes(x = Time_Bin, y = yvar, 
    linetype = oTreatment, fill = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cr"), size = 1, se = F) +
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    theme_bw() + xlab("Frame Number") + ylab("Bill Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "black", size = 1, linetype = "dashed") 

# Assessings uncertainty with posterior simulation

post_simGam(BillTemp_Mod_AR, Therm_Bin, view = "Time_Bin", group_var = "oTreatment", sim = 500)

# Distinct difference in slopes per group. Replotting marginal means with averaged points across
# birds.

Bill_Plot = ggplot(Bill_MMDat, aes(x = Time_Bin, y = yvar, 
    linetype = oTreatment, fill = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cr"), size = 1, se = F) +
    stat_summary_bin(data = Therm_Bin, aes(x = Time_Bin, y = Max.Bill),
        geom = "point", fun = "mean", binwidth = 10, pch = 21, size = 2) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "black", size = 1, linetype = "dashed") + 
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-420", "-210", "0", "210", "420"))

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Bill Temperature Response.jpeg", 
    Bill_Plot, dpi = 800, width = 7.67, height = 6.67, units = "in")

# Marginal means appear reasonably predictive, however, fall of handled curve appears
# to preceed that displayed by averaged points. 

# Plotting first order derivatives per group to assess change in response direction over times.

pred_dat = with(Therm_Bin, expand.grid("ID" = unique(ID), 
    "Time_Bin" = seq(0, 840, by = 10), "oTreatment" = unique(oTreatment)))

X0 = predict(BillTemp_Mod_AR, pred_dat, type = "lpmatrix")
increment = 1*10^-6
pred_dat_adj = pred_dat %>% mutate("Time_Bin" = Time_Bin + increment)
coefficients = BillTemp_Mod_AR$coefficients

var_covar = vcov(BillTemp_Mod_AR)
diag(var_covar)[7:16] = 0
sim_gam = rmvn_custom(10000, mu = coefficients, sig = var_covar)

X1 = predict(BillTemp_Mod_AR, pred_dat_adj, type = "lpmatrix")

Pred_X0 = X0 %*% t(sim_gam)
Pred_X1 = X1 %*% t(sim_gam)

Diff_Mat = Pred_X0*0

for (i in 1:ncol(Diff_Mat)){
    Diff_Mat[,i] = Pred_X1[,i] - Pred_X0[,i]
}

Diff_Mat_Corrected = sweep(Diff_Mat, 1, (60*10^6), FUN = "*")

Long_derivs_Aft = vector('list', ncol(Diff_Mat_Corrected))
for (i in 1:length(Long_derivs_Aft)){
     Long_derivs_Aft[[i]] = pred_dat %>% mutate("Pred" = Diff_Mat_Corrected[,i], "Sim_Number" = i)
}

Long_derivs_Aft = bind_rows(Long_derivs_Aft)

Long_derivs_Aft_Av = Long_derivs_Aft %>% group_by(Time_Bin, oTreatment) %>%
    summarise("Deriv" = mean(Pred), LCL = quantile(Pred, probs = c(0.025, 0.975), type = 8)[1],
    UCL = quantile(Pred, probs = c(0.025, 0.975), type = 8)[2])

Long_derivs_Aft_Av$Increasing = Long_derivs_Aft_Av$Deriv
Long_derivs_Aft_Av$Decreasing = Long_derivs_Aft_Av$Deriv
Long_derivs_Aft_Av$Increasing[c(which(Long_derivs_Aft_Av$Deriv < 0))] = NA
Long_derivs_Aft_Av$Decreasing[c(which(Long_derivs_Aft_Av$Deriv >= 0))] = NA

Bill_Derivs = ggplot(Long_derivs_Aft_Av, aes(x = Time_Bin, y = Deriv, 
    group = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) +
    geom_line(aes(y = Decreasing), size = 2.0, alpha = 0.5) + 
    geom_line(size = 0.5) + 
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature Slope (°C/min)") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    geom_vline(xintercept = 420, linetype = "longdash", colour = "black") +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-420", "-210", "0", "210", "420"))

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Bill Derivatives.jpeg", 
    Bill_Derivs, dpi = 800, width = 7.76, height = 6.67, units = "in")

# Preceeding fall is likely a consequence of model penalisation. Assessing trends 
# and mean data points. To do so, first computing simultaneous confidence intervals
# around coefficient estimates; refer to descriptions of each step provided for 
# plotting of eye temperature trends and 
# confidence intervals

var_cov = vcov(BillTemp_Mod_AR)
diag(var_cov)[7:16] = 0

# Expansion

pred_dat = with(Therm_Bin, expand.grid("ID" = unique(ID), 
    "Time_Bin" = seq(0, 840, by = 10), "oTreatment" = unique(oTreatment)))

pred_vals = predict(BillTemp_Mod_AR, newdata = pred_dat, se.fit = TRUE)
Cg = predict(BillTemp_Mod_AR, pred_dat, type = "lpmatrix")

fitted_se = pred_vals$se.fit

set.seed(20)
N_sim = 10000

BUdiff = rmvn_custom(N_sim, mu = rep(0, nrow(var_cov)), sig = var_cov)
simDeviation = Cg %*% t(BUdiff)
absDeviation = abs(sweep(simDeviation, 1, fitted_se, FUN = "/"))

# Max deviation distributions

MA_Cnt = apply(absDeviation[c(which(pred_dat$oTreatment == "Control")),], 2L, max)
MA_Hand = apply(absDeviation[c(which(pred_dat$oTreatment == "Handled")),], 2L, max)

par(mfrow=c(1,2))
hist(MA_Cnt) 
hist(MA_Hand) 
par(mfrow=c(1,1))

# Deviance values are right-skewed, so using type 8 quantiles again.

CV_Cnt = quantile(MA_Cnt, prob = 0.95, type = 8)
CV_Hand = quantile(MA_Hand, prob = 0.95, type = 8)

pred_df = transform(cbind(data.frame(pred_vals), pred_dat),
                  uprS = fit + (CV_Cnt * fitted_se),
                  lwrS = fit - (CV_Cnt * fitted_se))

pred_df$uprS[c(which(pred_dat$oTreatment == "Handled"))] = 
    pred_vals$fit[c(which(pred_dat$oTreatment == "Handled"))] + 
    (CV_Hand*fitted_se[c(which(pred_dat$oTreatment == "Handled"))])

pred_df$lwrS[c(which(pred_dat$oTreatment == "Handled"))] = 
    pred_vals$fit[c(which(pred_dat$oTreatment == "Handled"))] - 
    (CV_Hand*fitted_se[c(which(pred_dat$oTreatment == "Handled"))])

pred_df_grouped = pred_df %>% group_by(Time_Bin, oTreatment) %>%
    summarise(fit = mean(fit), se.fit = mean(se.fit), 
    UCL = mean(uprS), LCL = mean(lwrS))

# Plotting with simulated prediction curves 

Bill_Coef = BillTemp_Mod_AR$coefficients
set.seed(20)

var_cov = vcov(BillTemp_Mod_AR)
sim_gam = rmvn_custom(1000, mu = Bill_Coef, sig = var_cov)
mod_sims_fit = Cg %*% t(sim_gam)
ylims = range(mod_sims_fit)

Pred_Long_Form = vector('list', ncol(mod_sims_fit))

for (i in 1:length(Pred_Long_Form)){
    Pred_Long_Form[[i]] = pred_dat %>% mutate("Pred" = mod_sims_fit[,i], "Sim_Number" = i)
}

Pred_Long_Form = bind_rows(Pred_Long_Form)
Pred_Long_Form$Group = paste(Pred_Long_Form[,"oTreatment"], Pred_Long_Form$Sim_Number, sep = "_")

ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    geom_path(data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin, Group) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = Group),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("mediumvioletred", "mediumvioletred"), name = "Treatment") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw()

# Again, CIs appear suitable. Continuing with final plots.

Bill_plot_correct_CI = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210"))
 
showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Bill Correct CIs.jpeg", 
    Bill_plot_correct_CI, dpi = 800, width = 7.76, height = 6.67, units = "in")

Therm_Dat$oTreatment = Therm_Dat$OTreatment
Plot_Dat = Therm_Bin %>% 
    group_by(ID, oTreatment, Time_Bin) %>% 
    summarise("Max.Bill" = mean(Max.Bill, na.rm = T))

Bill_plot_with_points = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, alpha = pred_df_grouped$oTreatment), alpha = 0.55) +
    stat_summary_bin(geom = "point", fun = "mean", data = Plot_Dat, aes(x = Time_Bin, y = Max.Bill),
        binwidth = 10, pch = 21, size = 2, alpha = 0.6) +   
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1) +     
    annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 35.75) +  
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    scale_alpha_discrete(range = c(0.8, 0.4),name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    
 
showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Bill Correct CIs with points fine - Colour Matched with legend on panel.jpeg", Bill_plot_with_points, dpi = 800, width = 7, height = 6.67, units = "in")

# And again with confidence intervals around points 

Bill_plot_with_points_and_cis = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, alpha = pred_df_grouped$oTreatment), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Therm_Bin, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1) +     
    annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 35.75) +  
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    scale_alpha_discrete(range = c(0.8, 0.4),name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Bill_plot_with_points_and_cis

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Bill Correct CIs with points fine - Colour Matched with legend on panel - CIs included.jpeg", Bill_plot_with_points_and_cis, dpi = 800, width = 7, height = 6.67, units = "in")

# Lastly with grey box indicating duration of handling.

Bill_plot_with_points_and_cis = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, alpha = pred_df_grouped$oTreatment), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Therm_Bin, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1) +     
    annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 790, y = 35.5) +  
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    scale_alpha_discrete(range = c(0.8, 0.4),name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    annotate("rect", xmin = 420, xmax = 840, ymin = -Inf, ymax = Inf, colour = "black", fill = "grey40", alpha = 0.4) + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    theme(legend.position = c(0.21,0.1), legend.background = NULL)    

Bill_plot_with_points_and_cis

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - Bill Correct CIs with points fine - Colour Matched with legend on panel - CIs and Grey Box included.jpeg", Bill_plot_with_points_and_cis, dpi = 800, width = 7, height = 6.67, units = "in")

# Summarising model

summary(BillTemp_Mod_AR)

# Again, pulling out standard errors and mean coefficients, then comparing marginal means

Coefs = BillTemp_Mod_AR$coefficients
Coef.SE = as.data.frame(sqrt(diag(vcov(BillTemp_Mod_AR, unconditional = TRUE))))
Coef.SE$Row = c(1:nrow(Coef.SE))
Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Coef.SE))))

data.frame("Coefficient" = Coef_Names, "Estimate" = Coefs, "SE" = Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

Bill_Grid = ref_grid(BillTemp_Mod_AR,  
    at = list(Time_Bin = seq(420, 840, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin)

emmeans(Bill_Grid, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

# At beginning

Bill_Grid_Beginning = ref_grid(BillTemp_Mod_AR,  
    at = list(Time_Bin = seq(0, 419, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

emmeans(Bill_Grid_Beginning, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

Bill_Grid_Base_30 = ref_grid(BillTemp_Mod_AR,  
    at = list(Time_Bin = seq(0, 30, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

Bill_Grid_Base_60 = ref_grid(BillTemp_Mod_AR,  
    at = list(Time_Bin = seq(0, 60, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

emmeans(BillTemp_Mod_AR, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

emmeans(BillTemp_Mod_AR, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)
    
# Finally, testing for evidence of individualized responses to handling

BillTemp_Base = bam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_Bill$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = AR_Bill)

BillTemp_RandomSlope = bam(Max.Bill ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re") + 
    s(Time_Bin, ID, oTreatment, bs = "re"),
    AR.start = AR_Bill$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = AR_Bill)

test_stat = -2 * (logLik(BillTemp_Base) - logLik(BillTemp_RandomSlope))
pchisq(as.numeric(test_stat), df = (attr(logLik(BillTemp_RandomSlope), "df") - 
    attr(logLik(BillTemp_Base), "df")), lower.tail = FALSE)

# Plotting individual responses to treatments and adjusting for noise around baseline temperatures.

Exp_Bill = expand.grid("Time_Bin" = seq(1, 840, by = 1),
    "ID" = unique(na.omit(AR_Bill$ID)),
    "oTreatment" = unique(na.omit(AR_Bill$oTreatment))) %>% 
    mutate("Group" = paste(oTreatment, ID, sep = "_"))

# Again, removing combinations that did not occur

Compare = AR_Bill %>% 
    mutate("Group" = paste(oTreatment, ID, sep = "_"))

Catch = unique(Exp_Bill$Group)[c(which(!(unique(Exp_Bill$Group) %in% unique(Compare$Group))))]
Exp_Bill = Exp_Bill[-c(which(Exp_Bill$Group %in% Catch)),]

# Predicting

Bill_Response_Predictions = predict(BillTemp_RandomSlope, newdata = Exp_Bill, exclude = "s(ID)", se.fit = TRUE)

Bill_Ind = Exp_Bill %>%
    mutate("Fit" = Bill_Response_Predictions$fit,
        "LCL" = Bill_Response_Predictions$fit - 1.96*Bill_Response_Predictions$se.fit,
        "UCL" = Bill_Response_Predictions$fit + 1.96*Bill_Response_Predictions$se.fit) %>% 
    mutate("Groups" = paste(oTreatment, ID, sep = "_")) %>% 
    mutate(Groups = factor(Groups))

Bill_Individual_Slopes = ggplot(Bill_Ind, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, group = Groups), alpha = 0.15) +
    #stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 10, pch = 21, size = 2, alpha = 0.6) + 
    geom_line(aes(x = Time_Bin, y = Fit, group = Groups), size = 1, alpha = 0.7,
        colour = "black") +     
    annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 37) +  
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    scale_alpha_discrete(range = c(0.8, 0.4),name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)      

Bill_Individual_Slopes

extractLegend <- function(xplot) {
  grobs <- ggplot_gtable(ggplot_build(xplot))
  g_title <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[g_title]]
}

Legend_Pull = ggplot(Bill_Ind, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, group = Groups), alpha = 0.85) +
    #stat_summary_bin(geom = "point", fun = "mean", data = Therm_Bin_Rev, aes(x = Time_Bin, y = Max.Eye), binwidth = 10, pch = 21, size = 2, alpha = 0.6) + 
    geom_line(aes(x = Time_Bin, y = Fit, group = Groups), size = 1, alpha = 0.7,
        colour = "black") +     
    annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 37) +  
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    scale_alpha_discrete(range = c(0.8, 0.4),name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)   

NewLeg = extractLegend(Legend_Pull)
replace_Grob = ggplot_gtable(ggplot_build(Bill_Individual_Slopes))
rep_new = which(sapply(replace_Grob$grobs, function(x) x$name) == "guide-box")
replace_Grob$grobs[[rep_new]] = NewLeg

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/Individual_Responses_At_Bill.jpeg", replace_Grob, height = 7, width = 8, dpi = 800
)

# And assess trends with differences from baseline among individuals as response variable

TBin_Treat = split(Therm_Bin, f = Therm_Bin$Treatment)
for (i in 1:2){
    TBin_Treat[[i]] = split(TBin_Treat[[i]], f = TBin_Treat[[i]]$New.ID)
}  

for (j in 1:2){
    for (i in 1:length(TBin_Treat[[j]])){
        Base_Bill = TBin_Treat[[j]][[i]]$Max.Bill[which(TBin_Treat[[j]][[i]]$Time_Bin == min(TBin_Treat[[j]][[i]]$Time_Bin, na.rm = T))]
        TBin_Treat[[j]][[i]]$Bill_Diff = TBin_Treat[[j]][[i]]$Max.Bill - Base_Bill
    }
}    

for (i in 1:2){
    TBin_Treat[[i]] = bind_rows(TBin_Treat[[i]])
}

Bill_Diff = bind_rows(TBin_Treat)

# Modeling 

Bill_Diff$ID = factor(Bill_Diff$ID)

Bill_Diff_Mod = bam(Bill_Diff ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = Bill_Diff)

Bill_Diff$oTreatment = factor(Bill_Diff$oTreatment, ordered = FALSE)
Bill_Diff_AR = start_event(Bill_Diff, column="Time_Bin", event=c("New.ID"))
Bill_Diff_AR$oTreatment = factor(Bill_Diff$oTreatment, ordered = TRUE)
Bill_Diff$oTreatment = factor(Bill_Diff$oTreatment, ordered = TRUE)

Bill_Diff_Rho = acf(resid_gam(Bill_Diff_Mod), plot = F)$acf[2]
Bill_Diff_Rho

Bill_Diff_Mod_AR = bam(Bill_Diff ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = Bill_Diff_AR$start.event, rho = Bill_Diff_Rho,
    method = "REML", na.action = na.omit, 
    data = Bill_Diff_AR)

emmip(Bill_Diff_Mod_AR, oTreatment ~ Time_Bin, CIs = TRUE, 
    at = list(Time_Bin = seq(1, 840, 1)),
    data = Bill_Diff_AR)

summary(auto-Bill_Diff_Mod_AR)

# Odd. It appears that our residual autocorrelation is sapping up variance. Summarising original model and plotting to see.

summary(auto-Bill_Diff_Mod)

EM_Data_Diff = emmip(Bill_Diff_Mod_AR, oTreatment ~ Time_Bin, CIs = TRUE, 
    at = list(Time_Bin = seq(1, 840, 1)),
    data = Bill_Diff_AR, plot = FALSE)

ggplot(EM_Data_Diff, aes(x = Time_Bin, y = yvar, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.55, colour = NA) + 
    geom_line(size = 1, colour = "black") + 
    stat_summary_bin(data = Bill_Diff, aes(x = Time_Bin, y = Bill_Diff), geom = "errorbar", size = 1, width = 1, colour = "black",fun.data = "mean_cl_boot", binwidth = 50, position = position_dodge(width = 25)) + 
    stat_summary_bin(data = Bill_Diff, aes(x = Time_Bin, y = Bill_Diff), geom = "point", size = 3, pch = 21, colour = "black",fun = "mean", binwidth = 50, position = position_dodge(width = 25), alpha = 0.8) + theme_bw() +    
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + xlab("Time (s)") + ylab("Difference in Eye Temperature from Baseline (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme +  
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

# Indeed, trend-lines do not follow means at all, probably owing to high rho. Retaining original models.

#------------------------------------------------------------------------------#
## Modeling rate of heat loss by treatment group.
## These results are for interest sake.
#------------------------------------------------------------------------------#

# Assessing sample sizes 

Therm_Bin %>% filter(!is.na(m.qTot)) %>% group_by(Treatment) %>% summarise(n())
Therm_Bin %>% filter(!is.na(m.qTot)) %>% group_by(Treatment, New.ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())
Therm_Bin %>% filter(!is.na(m.qTot)) %>% group_by(Treatment, ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())

# Assessing distribution of heat transfer rate measurements 

ggplot(Therm_Bin, aes(x = m.qTot, fill = Treatment)) + 
    geom_histogram(binwidth = 0.002, colour = "black", alpha = 0.5) + 
    geom_density(aes(y=0.002 * ..count..), colour = "black", alpha = 0.5) + 
    theme_bw()

# By bird

ggplot(Therm_Bin, aes(x = ID, y = m.qTot, fill = Treatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1, 
        colour = "black", width = 0.3, position = position_dodge(width = 0.5)) +
    stat_summary(fun = "mean", geom = "point", size = 3, colour = "black",
        position = position_dodge(width = 0.5), pch = 21) +         
    theme_bw()

# Bimodal in both groups, perhaps according to time? TZ22 also particularly high.
# Proceeding with Gaussian distribution, then assessing residual distribution. 
# First, however, converting W to mW.

Therm_Bin$qTot = Therm_Bin$m.qTot*1000

# Assuming autocorrelation, and preparing to correct for such.

HT_Mod = bam(qTot ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),  
    method = "REML", na.action = na.omit, data = Therm_Bin)

AR_HT = start_event(Therm_Bin, column="Time_Bin", event=c("New.ID"))
Rho = acf(resid_gam(HT_Mod), plot = F)$acf[2]
Rho

# High autocorrelation confirm. Proceeding with corrected model.

HT_Mod_AR = bam(qTot ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_HT$start.event, rho = Rho,
    method = "REML", na.action = na.omit, 
    data = AR_HT
)

ACF_Plot(model = HT_Mod_AR, residual_type = "normalised")

# Autocorrelation largely dissolved. Continuing with model assessments.

{
p1 = ggplot(data = data.frame("Res" = resid_gam(HT_Mod_AR)),
  aes(x = 1:length(resid_gam(HT_Mod_AR)), y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalised Residuals")

p2 = ggplot(data = data.frame("Res" = resid_gam(HT_Mod_AR, incl_na = TRUE), 
  "Fit" = predict(HT_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)))),
  aes(x = Fit, y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalised Residuals")

p3 = ggplot(data.frame("Res" = resid_gam(HT_Mod_AR)), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(resid_gam(HT_Mod_AR)), 
    sd = sd(resid_gam(HT_Mod_AR)))) +
  theme_bw() + xlab("Normalised Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(resid_gam(HT_Mod_AR)) - 
    3*sd(resid_gam(HT_Mod_AR))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(resid_gam(HT_Mod_AR)) + 
    3*sd(resid_gam(HT_Mod_AR))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = resid_gam(HT_Mod_AR, incl_na = TRUE), 
    "Fit" = predict(HT_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

p5 = as.data.frame(Therm_Bin %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)) %>% 
  mutate("Pred" = predict(HT_Mod_AR, newdata = as.data.frame(Therm_Bin %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)))) %>%
    ggplot(aes(x = Pred, y = qTot)) + 
    geom_point(size = 2, colour = "black", fill = "orchid", pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = FALSE) + 
    my.theme + theme_bw()

p6 = ACF_Plot(model = HT_Mod_AR, residual_type = "normalised")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, top = "Model Residuals")
}    

# Significant left-tailing. Trying to identify.

Group_Diag(model = HT_Mod_AR, cutoff = (mean(resid_gam(HT_Mod_AR, incl_na = FALSE)) + 3*sd(resid_gam(HT_Mod_AR, incl_na = FALSE))),
    data = as.data.frame(Therm_Bin %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)), group = "ID", view = "Time_Bin")

Group_Diag(model = HT_Mod_AR, cutoff = (mean(resid_gam(HT_Mod_AR, incl_na = FALSE)) + 3*sd(resid_gam(HT_Mod_AR, incl_na = FALSE))),
    data = as.data.frame(Therm_Bin %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)), label = "oTreatment", group = "ID", view = "Time_Bin")

# Outliers appear evenly dispersed across individuals, times, and groups. Plotting residuals by predictors.

p1 = as.data.frame(Therm_Bin %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = oTreatment, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p2 = as.data.frame(Therm_Bin %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Res, fill = oTreatment)) + 
    geom_density(alpha = 0.5, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p3 = as.data.frame(Therm_Bin %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = Time_Bin, y = Res, fill = oTreatment, linetype = oTreatment)) + 
    stat_summary_bin(fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.5, 
    binwidth = 50) + stat_summary_bin(fun = "mean", geom = "point", size = 2, pch = 21,
        colour = "black", binwidth = 10) + 
    geom_smooth(method = "lm", size = 1, colour = "black") +    
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p4 = as.data.frame(Therm_Bin %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR, incl_na = TRUE)) %>% 
        ggplot(aes(x = ID, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4, position = position_dodge(width = 0.5)) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black",
        position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Residual Distributions")

# TZ19 again highly variable. Removing this individual and proceeding.

Therm_Rev = Therm_Bin %>% filter(ID != "TZ19") %>% 
    mutate("ID" = factor(ID))

# Sample sizes

Therm_Rev %>% filter(!is.na(qTot)) %>% group_by(Treatment) %>% summarise(n())
Therm_Rev %>% filter(!is.na(qTot)) %>% group_by(Treatment, New.ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())
Therm_Rev %>% filter(!is.na(qTot)) %>% group_by(Treatment, ID) %>% 
    summarise(n()) %>% group_by(Treatment) %>% summarise(n())

HT_Mod = bam(qTot ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),  
    method = "REML", na.action = na.omit, 
    data = Therm_Rev)

AR_HT_Rev = start_event(Therm_Rev, column="Time_Bin", event=c("New.ID"))
Rho_Rev = acf(resid_gam(HT_Mod), plot = F)$acf[2]
Rho_Rev

HT_Mod_AR_Rev = bam(qTot ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    AR.start = AR_HT_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_HT_Rev
)

# Residual diagnostics.

{
p1 = ggplot(data = data.frame("Res" = resid_gam(HT_Mod_AR_Rev)),
  aes(x = 1:length(resid_gam(HT_Mod_AR_Rev)), y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "mediumseagreen", alpha = 0.5) + xlab("Row Number") +
  ylab("Normalised Residuals")

p2 = ggplot(data = data.frame("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE), 
  "Fit" = predict(HT_Mod_AR_Rev, newdata = as.data.frame(AR_HT_Rev %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)))),
  aes(x = Fit, y = Res)) + theme_bw() +
  geom_point(size = 3, colour = "cornflowerblue", alpha = 0.5) + xlab("Y hat") +
  ylab("Normalised Residuals")

p3 = ggplot(data.frame("Res" = resid_gam(HT_Mod_AR_Rev)), aes(Res)) +
  geom_histogram(alpha = 0.5, colour = "black", fill = "mediumorchid",
    aes(y=..density.., fill=..count..)) +
  stat_function(fun = dnorm, size = 1,
    args = list(mean = mean(resid_gam(HT_Mod_AR_Rev)), 
    sd = sd(resid_gam(HT_Mod_AR_Rev)))) +
  theme_bw() + xlab("Normalised Residuals") + ylab("Count") +
  geom_vline(xintercept = (mean(resid_gam(HT_Mod_AR_Rev)) - 
    3*sd(resid_gam(HT_Mod_AR_Rev))),
    size = 1, linetype = "dashed") +
  geom_vline(xintercept = (mean(resid_gam(HT_Mod_AR_Rev)) + 
    3*sd(resid_gam(HT_Mod_AR_Rev))),
    size = 1, linetype = "dashed")

p4 = ggplot(data.frame("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE), 
    "Fit" = predict(HT_Mod_AR_Rev, newdata = as.data.frame(AR_HT_Rev %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)))), 
    aes(sample=Res)) + stat_qq(colour = "gold") + 
    stat_qq_line() + my.theme + theme_bw()

p5 = as.data.frame(AR_HT_Rev %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)) %>% 
  mutate("Pred" = predict(HT_Mod_AR_Rev, newdata = as.data.frame(AR_HT_Rev %>% 
  dplyr::select(qTot, oTreatment, Time_Bin, ID) %>% drop_na(.)))) %>%
    ggplot(aes(x = Pred, y = qTot)) + 
    geom_point(size = 2, colour = "black", fill = "orchid", pch = 21) + 
    geom_smooth(method = "lm", size = 1, colour = "black", se = FALSE) + 
    my.theme + theme_bw()

p6 = ACF_Plot(model = HT_Mod_AR_Rev, residual_type = "normalised")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, top = "Model Residuals")
}    

# Tailing persists, but not so severe. Plotting by predictors.

p1 = as.data.frame(AR_HT_Rev %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE)) %>% 
        ggplot(aes(x = oTreatment, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p2 = as.data.frame(AR_HT_Rev %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE)) %>% 
        ggplot(aes(x = Res, fill = oTreatment)) + 
    geom_density(alpha = 0.5, colour = "black") + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p3 = as.data.frame(AR_HT_Rev %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE)) %>% 
        ggplot(aes(x = Time_Bin, y = Res, fill = oTreatment, linetype = oTreatment)) + 
    stat_summary_bin(fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.5, 
    binwidth = 50) + stat_summary_bin(fun = "mean", geom = "point", size = 2, pch = 21,
        colour = "black", binwidth = 10) + 
    geom_smooth(method = "lm", size = 1, colour = "black") +    
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

p4 = as.data.frame(AR_HT_Rev %>% dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE)) %>% 
        ggplot(aes(x = ID, y = Res, fill = oTreatment)) + 
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", size = 1,
        colour = "black", width = 0.4, position = position_dodge(width = 0.5)) + 
    stat_summary(fun = "mean", size = 2, pch = 21, colour = "black",
        position = position_dodge(width = 0.5)) + 
    scale_fill_manual(values = c("black", "grey70")) +             
    theme_bw()

grid.arrange(p1, p2, p3, p4, nrow = 2, top = "Residual Distributions")

# Assessing whether tails are problematic. 

mod_fits = fitted(HT_Mod_AR_Rev, type = "response")

leave_out = function(x) {
  updated_data = AR_HT_Rev[-x, ]
  model = bam(qTot ~ oTreatment +
    s(Time_Bin, bs = "tp", k = 3, fx = TRUE) + 
    s(Time_Bin, by = oTreatment, bs = "tp", k = 3, m = 1, fx = TRUE) + 
    s(ID, bs = "re"),
    AR.start = AR_HT_Rev$start.event, rho = Rho_Rev,
    method = "REML", na.action = na.omit, 
    data = AR_HT_Rev
)
  sum((mod_fits[-x] - fitted(model))^2)
}

influence = sapply(1:nrow(AR_HT_Rev), leave_out)
AR_HT_Rev$Influence = influence

Crit_I = AR_HT_Rev %>% filter(!is.na(qTot) & start.event == "FALSE") %>%
    summarise(Crit = 3*mean(Influence, na.rm = T),
        Weak_Crit = mean(Influence, na.rm = T) + 3*sd(Influence, na.rm = T))

ggplot(as.data.frame(AR_HT_Rev %>% filter(!is.na(qTot) & start.event == "FALSE")), 
    aes(x = Time_Bin, y = Influence, fill = oTreatment, linetype = oTreatment)) + 
    geom_point(pch = 21, size = 2, colour = "black") + 
    geom_smooth(method = "lm", size = 1, colour = "black") + 
    theme_bw() + geom_hline(yintercept = Crit_I$Crit, linetype = "longdash",
        size = 1, colour = "mediumvioletred") + 
        geom_hline(yintercept = Crit_I$Weak_Crit, linetype = "longdash",
          size = 1, colour = "mediumvioletred") +
    scale_fill_manual(values = c("black", "grey70"))

# Little concern; no individual points appear to hold significant influence. 
# Assessing coefficients

Coef_Dist(model = HT_Mod_AR_Rev, incl_int = FALSE)

# Interesting; only the first time bin appears to differ between treatments. 
# Checking whether smooths have sufficient knots.

gam.check(HT_Mod_AR_Rev)

# In appears so, although the time spline itself is largely linear.

# Plotting marginal means.

HT_MMDat = emmip(HT_Mod_AR_Rev, oTreatment ~ Time_Bin, data = AR_HT_Rev,
    cov.reduce = FALSE, nesting = NULL, nesting.order = FALSE,
    CIs = TRUE, plot = F)

ggplot(HT_MMDat, aes(x = Time_Bin, y = yvar, 
    linetype = oTreatment, fill = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cr"), size = 1, se = F) +
    stat_summary_bin(data = AR_HT_Rev, aes(x = Time_Bin, y = qTot),
        geom = "point", fun = "mean", binwidth = 10, pch = 21, size = 2) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    theme_bw() + xlab("Frame Number") + ylab("Heat Transfer (W)") + 
    geom_vline(xintercept = 420, colour = "black", size = 1, linetype = "dashed")

# Plotting with mean centred residuals

ggplot(HT_MMDat, aes(x = Time_Bin, y = yvar, 
    linetype = oTreatment, fill = oTreatment, colour = oTreatment)) + 
    geom_ribbon(aes(x = Time_Bin, ymin = LCL, ymax = UCL), alpha = 0.5) + 
    geom_smooth(method = "gam", formula = y~s(x, k = 4, bs = "cr"), size = 1, se = F) +
    stat_summary_bin(data = as.data.frame(AR_HT_Rev %>% 
        dplyr::select(qTot, Time_Bin, oTreatment, ID) %>% 
        drop_na(.)) %>% mutate("Res" = resid_gam(HT_Mod_AR_Rev, incl_na = TRUE) + 
        mean(AR_HT_Rev$qTot, na.rm = TRUE)), 
        aes(x = Time_Bin, y = Res),
        geom = "point", fun = "mean", binwidth = 10, pch = 21, size = 2) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") + 
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment") + 
    theme_bw() + xlab("Frame Number") + ylab("Heat Transfer (mW)") + 
    geom_vline(xintercept = 420, colour = "black", size = 1, linetype = "dashed")

# No distinct pull. Plotting model coefficients.

var_cov = vcov(HT_Mod_AR_Rev)
diag(var_cov)[7:15] = 0

# Expansion

pred_dat = with(as.data.frame(Therm_Bin %>% filter(ID != "TZ19") %>% mutate("ID" = factor(ID))), 
    expand.grid("ID" = unique(ID), 
    "Time_Bin" = seq(0, 840, by = 10), "oTreatment" = unique(oTreatment)))

pred_vals = predict(HT_Mod_AR_Rev, newdata = pred_dat, se.fit = TRUE)
Cg = predict(HT_Mod_AR_Rev, pred_dat, type = "lpmatrix")
fitted_se = pred_vals$se.fit

set.seed(20)
N_sim = 10000

BUdiff = rmvn_custom(N_sim, mu = rep(0, nrow(var_cov)), sig = var_cov)
simDeviation = Cg %*% t(BUdiff)
absDeviation = abs(sweep(simDeviation, 1, fitted_se, FUN = "/"))
MA_Cnt = apply(absDeviation[c(which(pred_dat$oTreatment == "Control")),], 2L, max)
MA_Hand = apply(absDeviation[c(which(pred_dat$oTreatment == "Handled")),], 2L, max)

par(mfrow=c(1,2))
hist(MA_Cnt) 
hist(MA_Hand) 
par(mfrow=c(1,1))

# Deviance values, again, are right-skewed.

CV_Cnt = quantile(MA_Cnt, prob = 0.95, type = 8)
CV_Hand = quantile(MA_Hand, prob = 0.95, type = 8)
pred_df = transform(cbind(data.frame(pred_vals), pred_dat),
                  uprS = fit + (CV_Cnt * fitted_se),
                  lwrS = fit - (CV_Cnt * fitted_se))
pred_df$uprS[c(which(pred_dat$oTreatment == "Handled"))] = 
    pred_vals$fit[c(which(pred_dat$oTreatment == "Handled"))] + 
    (CV_Hand*fitted_se[c(which(pred_dat$oTreatment == "Handled"))])
pred_df$lwrS[c(which(pred_dat$oTreatment == "Handled"))] = 
    pred_vals$fit[c(which(pred_dat$oTreatment == "Handled"))] - 
    (CV_Hand*fitted_se[c(which(pred_dat$oTreatment == "Handled"))])
pred_df_grouped = pred_df %>% group_by(Time_Bin, oTreatment) %>%
    summarise(fit = mean(fit), se.fit = mean(se.fit), 
    UCL = mean(uprS), LCL = mean(lwrS))

# Posterior simulation

set.seed(20)

HT_Coefs = HT_Mod_AR_Rev$coefficients
var_cov = vcov(HT_Mod_AR_Rev)
sim_gam = rmvn_custom(1000, mu = HT_Coefs, sig = var_cov)
mod_sims_fit = Cg %*% t(sim_gam)
ylims = range(mod_sims_fit)

Pred_Long_Form = vector('list', ncol(mod_sims_fit))

for (i in 1:length(Pred_Long_Form)){
    Pred_Long_Form[[i]] = pred_dat %>% mutate("Pred" = mod_sims_fit[,i], "Sim_Number" = i)
}

Pred_Long_Form = bind_rows(Pred_Long_Form)
Pred_Long_Form$Group = paste(Pred_Long_Form[,"oTreatment"], Pred_Long_Form$Sim_Number, sep = "_")

ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    geom_path(data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin, Group) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = Group),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("mediumvioletred", "mediumvioletred"), name = "Treatment") +
    scale_colour_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw()

# CIs look fitting.

HT_correct_CI = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) + 
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Heat Transfer (mW)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210"))

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - HT Correct CIs.jpeg", 
    HT_correct_CI, dpi = 800, width = 7.76, height = 6.67, units = "in")

Therm_Dat$oTreatment = Therm_Dat$Treatment
Therm_Dat$qTot = Therm_Dat$qTot*1000

HT_with_points = ggplot(pred_df_grouped, aes(x = Time_Bin, 
    fill = oTreatment, linetype = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    stat_summary(geom = "smooth", fun = "mean", method = "gam", formula = y~s(x, k = 4),
        data = as.data.frame(Pred_Long_Form %>% 
        group_by(oTreatment, Time_Bin) %>%
        summarise(yvar = mean(Pred))), 
        aes(x = Time_Bin, y = yvar, colour = oTreatment, group = oTreatment),
        size = 1, alpha = 0.4) + 
    stat_summary_bin(data = Therm_Dat, binwidth = 5, geom = "point",
        fun = "mean", aes(x = Slice, y = qTot, fill = oTreatment),
        pch = 21) +     
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    theme_bw() + xlab("Time (s)") + ylab("Heat Transfer (mW)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    my.theme + scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210"))

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/COLI_Thermal/Figures/COLI - HT Correct CIs with points.jpeg", 
    HT_with_points, dpi = 800, width = 7.76, height = 6.67, units = "in")

# Summarising model

summary(HT_Mod_AR_Rev)

# Calculating marginal means 

HT_Grid = ref_grid(HT_Mod_AR_Rev,  
    at = list(Time_Bin = seq(420, 840, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

HT_Grid_Before = ref_grid(HT_Mod_AR_Rev,  
    at = list(Time_Bin = seq(0, 419, by = 1)), 
    cov.reduce = FALSE, type = "response",
    data = Therm_Bin_Rev)

emmeans(HT_Grid, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

emmeans(HT_Grid_Before, specs = pairwise ~ oTreatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Therm_Bin_Rev)

#-----------------------------------------------------------------------------# 
# To preserve: Dropping of random effects # 
#-----------------------------------------------------------------------------#

vcov_Dropped = vcov(BillTemp_Mod_AR)
vcov_Dropped[,grep("ID", colnames(vcov_Dropped))] = vcov_Dropped[,grep("ID", colnames(vcov_Dropped))]*0
vcov_Dropped[grep("ID", colnames(vcov_Dropped)),] = vcov_Dropped[grep("ID", colnames(vcov_Dropped)),]*0

pred_dat = with(Therm_Bin, expand.grid("ID" = unique(ID), 
    "Time_Bin" = seq(0, 840, by = 10), "oTreatment" = unique(oTreatment)))

pred_vals = predict.gam(BillTemp_Mod_AR, newdata = pred_dat, exclude = "s(ID)", 
    se.fit = TRUE, type = "terms")
pred_fit = rowSums(pred_vals$fit)

pred_fit = data.frame("fit" = pred_fit + as.numeric(attr(pred_vals, "constant")))

# Estimating standard errors from lpmatix with dropped random effects

Lp_Dropped = predict(BillTemp_Mod_AR, newdata = pred_dat, type = "lpmatrix")
Lp_Dropped[,grep("ID", colnames(Lp_Dropped))] = Lp_Dropped[,grep("ID", colnames(Lp_Dropped))]*0

se_yhat = as.numeric(sqrt(diag(Lp_Dropped)%*%Vcov_dropped%*%t(Lp_Dropped)))

# Assessing change in standard errors 

plot(test_se ~ se_yhat)

# Odd trends, but absolute change is limited.

set.seed(20)
N_sim = 10000

# Setting mean of simulated multivariate normal distribution to zero, 
# with variance estimated from our model.

BUdiff = rmvn(N_sim, mu = rep(0, nrow(Vcov_dropped)), sig = Vcov_dropped)

# Calcuting deviation from predicted model across simulations, then calculating
# absolute deviation from true model using predicted standard errors.

simDeviation = Lp_Dropped %*% t(BUdiff)
absDeviation = abs(sweep(simDeviation, 1, se_yhat, FUN = "/"))

# Maximal absolute deviation is pulled from each simulation-model comparison, then 
# critical value estimated from maximal deviations (95% confidence is assumed here).

Max_Abs_Deviation = apply(absDeviation, 2L, max)
hist(Max_Abs_Deviation) 

# Right-skewed; using type 8 quantiles.

Crit_vals = quantile(Max_Abs_Deviation, prob = 0.95, type = 9)

# Creating data frame to plot confidence intervals, then averaging at time intervals
# by treatment.

pred_df = transform(cbind(pred_dat, pred_fit),
    uprS = fit + (Crit_vals * se_yhat),
    lwrS = fit - (Crit_vals * se_yhat))

pred_df_grouped = pred_df %>% group_by(Time_Bin, oTreatment) %>%
    summarise(fit = mean(fit),  
    UCL = mean(uprS), LCL = mean(lwrS))

# Plotting 
dev.new()
ggplot(pred_df_grouped, aes(x = Time_Bin, fill = oTreatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.2) +
    scale_fill_manual(values = c("grey70", "grey30"), name = "Treatment") +
    theme_bw()

# CIs appear even broader! Perhaps simultaneous confidence intervals should
# not be calculated in this way. 



