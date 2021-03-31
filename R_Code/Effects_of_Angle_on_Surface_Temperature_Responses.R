#######################################################################################
# Effects of head orientation on surface temperature measurements in Domestic Pigeons #
#######################################################################################

# Loading packages and additional functions

library('easypackages')
libraries('car', 'colRoz','itsadug','emmeans','mgcv','tidyverse')

getFrom = function (pkg, name){ 
    get(name, envir = asNamespace(pkg), inherits = FALSE)
}

matmultdiag = function (x, y, ty = t(y)) {
    if (ncol(x) != ncol(ty)) 
        stop("non-conformable arguments")
    if (nrow(x) != nrow(ty)) 
        stop("result is not a square matrix")
    return(rowSums(x * ty))
}

pred_gls = function (object, newdata, se.fit = FALSE, na.action = na.omit) {
    if (missing(newdata) && !se.fit) {
        return(fitted(object))
    }
    form <- getFrom("nlme", "getCovariateFormula")(object)
    mfArgs <- list(formula = form, data = newdata, na.action = na.action)
    #mfArgs$drop.unused.levels <- FALSE
    dataMod <- do.call(model.frame, mfArgs)
    contr <- object$contrasts
    for (i in names(dataMod)) {
        if (inherits(dataMod[, i], "factor") && !is.null(contr[[i]])) {
            levs <- levels(dataMod[, i])
            levsC <- dimnames(contr[[i]])[[1]]
            if (any(wch <- is.na(match(levs, levsC)))) {
                stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s", 
                  "levels %s not allowed for %s"), paste(levs[wch], 
                  collapse = ",")), domain = NA)
            }
            attr(dataMod[, i], "contrasts") <- contr[[i]][levs, 
                , drop = FALSE]
        }
    }
    N <- nrow(dataMod)
    if (length(all.vars(form)) > 0) {
        X <- model.matrix(form, dataMod)
    } else {
        X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
    }
    cf <- coef(object)
    val <- X %*% cf
    if (se.fit) {
        se <- sqrt(matmultdiag(X %*% vcov(object), ty = X))
        val <- list(fit = val, se.fit = unname(se))
    }
    lab <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        lab <- paste(lab, aux)
    }
    structure(val, label = lab)
}

# Adding plot theme 

my.theme = theme(
  panel.grid.minor = element_blank(),
  panel.grid.major =
    element_line(colour = "grey75"),
  axis.title = element_text(size = 14, family = "Noto Sans"),
  axis.text = element_text(size = 13, colour = "black", family = "Noto Sans"),
  axis.title.y = element_text(vjust = 1), 
  legend.title = element_text(
    size = 14,
    colour = "black", family = "Noto Sans"
  ), legend.text = element_text(
    size = 14,
    colour = "black", family = "Noto Sans"
  )
)

# Reading data pertaining to orientation estimates

{
    T6 = read.csv("/home/joshk/Desktop/All/T6/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T18 = read.csv("/home/joshk/Desktop/All/T18/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T11 = read.csv("/home/joshk/Desktop/All/T11/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T9 = read.csv("/home/joshk/Desktop/All/T9/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T22 = read.csv("/home/joshk/Desktop/All/T22/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T24 = read.csv("/home/joshk/Desktop/All/T24/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T8 = read.csv("/home/joshk/Desktop/All/T8/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))
    T24_2 = read.csv("/home/joshk/Desktop/All/T24_2/Rotation_Out.csv") %>% 
        mutate(Keep = "TRUE", Slice = seq(1,nrow(.),1))

    T6_S = read.csv("/home/joshk/Desktop/All/T6/TZ6_Handling.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T18_S = read.csv("/home/joshk/Desktop/All/T18/TZ18_Control_Adj.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T11_S = read.csv("/home/joshk/Desktop/All/T11/TZ11_Handling_Adj.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T9_S = read.csv("/home/joshk/Desktop/All/T9/TZ9_Control_Adj.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T22_S = read.csv("/home/joshk/Desktop/All/T22/TZ22_Control.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T24_S = read.csv("/home/joshk/Desktop/All/T24/TZ24_Handling.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T8_S = read.csv("/home/joshk/Desktop/All/T8/TZ8_Handling.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )
    T24_2_S = read.csv("/home/joshk/Desktop/All/T24_2/TZ24_Control.csv") %>% 
        select(-c(X,Initials,Slice)) %>% mutate("Sums" = rowSums(.), 
        "X_Dir" = ifelse(CB_x == -1 | BT_x == -1, NA, CB_x - BT_x),
        "Y_Dir" = ifelse(LEF_y == -1 | LER_y == -1, 
        ifelse(REF_y == -1 | RER_y == -1, NA, RER_y - REF_y),
        LER_y - LEF_y)
        )

    # Labelling rows where orientation was not estimated

    T6$Keep[c(which(T6_S$Sums == 0))] = "FALSE"
    T18$Keep[c(which(T18_S$Sums == 0))] = "FALSE"
    T11$Keep[c(which(T11_S$Sums == 0))] = "FALSE"
    T9$Keep[c(which(T9_S$Sums == 0))] = "FALSE"
    T22$Keep[c(which(T22_S$Sums == 0))] = "FALSE"
    T24$Keep[c(which(T24_S$Sums == 0))] = "FALSE"
    T8$Keep[c(which(T8_S$Sums == 0))] = "FALSE"
    T24_2$Keep[c(which(T24_2_S$Sums == 0))] = "FALSE"

    T6$X_Dir = T6_S$X_Dir
    T18$X_Dir = T18_S$X_Dir
    T11$X_Dir = T11_S$X_Dir
    T9$X_Dir = T9_S$X_Dir
    T22$X_Dir = T22_S$X_Dir
    T24$X_Dir = T24_S$X_Dir
    T8$X_Dir = T8_S$X_Dir
    T24_2$X_Dir = T24_2_S$X_Dir

    T6$Y_Dir = T6_S$Y_Dir
    T18$Y_Dir = T18_S$Y_Dir
    T11$Y_Dir = T11_S$Y_Dir
    T9$Y_Dir = T9_S$Y_Dir
    T22$Y_Dir = T22_S$Y_Dir
    T24$Y_Dir = T24_S$Y_Dir
    T8$Y_Dir = T8_S$Y_Dir
    T24_2$Y_Dir = T24_2_S$Y_Dir
}

# Reading in surface temperature data

{
    Pair_T6 = read.csv("/home/joshk/Desktop/All/T6/TZ6_Handling_ST.csv")    
    Pair_T18 = read.csv("/home/joshk/Desktop/All/T18/TZ18_Control_ST.csv")
    Pair_T11 = read.csv("/home/joshk/Desktop/All/T11/TZ11_Handling_ST.csv")
    Pair_T9 = read.csv("/home/joshk/Desktop/All/T9/TZ9_Control_ST.csv")
    Pair_T22 = read.csv("/home/joshk/Desktop/All/T22/TZ22_Control_ST.csv")
    Pair_T24 = read.csv("/home/joshk/Desktop/All/T24/TZ24_Handling_ST.csv")
    Pair_T8 = read.csv("/home/joshk/Desktop/All/T8/TZ8_Handling_ST.csv")
    Pair_T24_2 = read.csv("/home/joshk/Desktop/All/T24_2/TZ24_Control_ST.csv")
}

# Binding orientation and surface temperature data together

{
    T6_B = left_join(Pair_T6 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T6 %>% select(-c(X)), by = c("Slice"))
    T18_B = left_join(Pair_T18 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T18 %>% select(-c(X)), by = c("Slice"))
    T11_B = left_join(Pair_T11 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T11 %>% select(-c(X)), by = c("Slice"))
    T9_B = left_join(Pair_T9 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T9 %>% select(-c(X)), by = c("Slice"))     
    T22_B = left_join(Pair_T22 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T22 %>% select(-c(X)), by = c("Slice"))
    T24_B = left_join(Pair_T24 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T24 %>% select(-c(X)), by = c("Slice"))   
    T8_B = left_join(Pair_T8 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T8 %>% select(-c(X)), by = c("Slice"))
    T24_2_B = left_join(Pair_T24_2 %>% select(-c(Xmin,Xmax,Xmean,Ymin,Ymax,Ymean)), 
    T24_2 %>% select(-c(X)), by = c("Slice")) %>% mutate(ID = "TZ24_2")     
}

# And ensuring that rotation values for rows wherein orientation was not estimated are "NA".

All = rbind(T6_B,T18_B,T11_B,T9_B,T22_B,T24_B,T8_B,T24_2_B) %>% 
    mutate(R1 = ifelse(Keep == "FALSE", NA, R1),
        R2 = ifelse(Keep == "FALSE", NA, R2),
        R3 = ifelse(Keep == "FALSE", NA, R3),
        P = ifelse(Keep == "FALSE", NA, P),
        Y = ifelse(Keep == "FALSE", NA, Y),
        R = ifelse(Keep == "FALSE", NA, R)
    )

# Assessing spread of euclidean residuals 

ggplot(All, aes(x = MEuc)) + 
    geom_histogram(colour = "black", fill = "slateblue") +
    theme_bw() + facet_wrap(~ID)

# Quite wide. A cut-off between 10 and 25 pixels appears reasonable. Comparing.

ggplot(subset(All, MEuc < 25), aes(x = MEuc)) + 
    geom_histogram(colour = "black", fill = "slateblue") +
    theme_bw() + facet_wrap(~ID)

ggplot(subset(All, MEuc < 10), aes(x = MEuc)) + 
    geom_histogram(colour = "black", fill = "slateblue") +
    theme_bw() + facet_wrap(~ID)

# Fair. Retaining all orientation estimates with mean euclidean residuals below 25 pixels (or ~ 1 cm). Assessing sample sizes according to this cut-off.

All %>% filter(!is.na(Y)) %>%
    mutate("Group" = ifelse(MEuc > 25, "Out", "In")) %>%
    group_by(Group) %>% 
    summarise("Count" = n())

# Only 86 estimates removed. 

All %>% filter(!is.na(Y)) %>%
    mutate("Group" = ifelse(MEuc > 25, "Out", "In")) %>%
    filter(Group == "In") %>% 
    group_by(ID) %>% 
    summarise("Count" = n()) %>% 
    summarise("Mean_Count" = mean(Count),
    "SD Count" = sd(Count))

# Cleaning according to cut-off criterea. Note that when mean euclidean distance exceeds 25 pixels, the relevant rotational estimates are reduced to "NA" values.

All = All %>% 
    mutate(ID = factor(ID), Treatment = factor(Treatment, ordered = TRUE)) %>% 
    mutate(R1 = ifelse(MEuc > 25, NA, R1),
        R2 = ifelse(MEuc > 25, NA, R2),
        R3 = ifelse(MEuc > 25, NA, R3),
        P = ifelse(MEuc > 25, NA, P),
        Y = ifelse(MEuc > 25, NA, Y),
        R = ifelse(MEuc > 25, NA, R),
        X_Dir = ifelse(MEuc > 25, NA, X_Dir),
        Y_Dir = ifelse(MEuc > 25, NA, Y_Dir)
        )

# Next, making sure that onset of handling is rooted at 420 seconds for all handled birds.

subset(All[c(which(All$Handling != lag(All$Handling, 1))),], Treatment == "Handled" & Handling == 1) 

Bird_Split = split(All, f = All$ID)
names(Bird_Split)
{
    Bird_Split[[1]]$Slice = Bird_Split[[1]]$Slice - 71
    Bird_Split[[4]]$Slice = Bird_Split[[4]]$Slice - 8
    Bird_Split[[7]]$Slice = Bird_Split[[7]]$Slice - 10
    Bird_Split[[2]]$Slice = Bird_Split[[2]]$Slice - 88
}

Rooted = bind_rows(Bird_Split) %>% 
    filter(Slice >= 0 & Slice <= 840)

# Checking correction.

subset(Rooted[c(which(Rooted$Handling != lag(Rooted$Handling, 1))),], Treatment == "Handled" & Handling == 1) 

# Good. Double-checking that axes of rotation (pertaining to orientation estimates) are indeed correct.

P1 = Rooted %>% mutate(Y_abs = abs(Y)) %>% 
    ggplot(aes(x = X_Dir, y = Y_abs)) + 
    geom_smooth(colour = "black") + 
    stat_summary_bin(geom = "point", fun = "mean", binwidth = 0.5, 
    size = 3, pch = 21, alpha = 0.7, colour = "black", fill = "lightseagreen") + 
    theme_classic() + xlab("Cyr Base x - Bill Tip x") + ylab("|Yaw|") + 
    theme(axis.text = element_text(family = "Noto Sans", size = 11, colour = "black"),
    axis.title = element_text(family = "Noto Sans", size = 14, colour = "black")) + 
    ylim(c(0,200)) + 
    annotate("text", x = 68, y = 0, label = "Pointing left", size = 4, family = "Noto Sans", colour = "grey20") + 
    annotate("text", x = -55, y = 0, label = "Pointing right", size = 4, family = "Noto Sans", colour = "grey20")

P2 = ggplot(Rooted, aes(x = Y_Dir, y = P)) + 
    geom_smooth(method = "lm", colour = "black") + 
    stat_summary_bin(geom = "point", fun = "mean", binwidth = 0.5, 
    size = 3, pch = 21, alpha = 0.7, colour = "black", fill = "lightseagreen") + 
    theme_classic() + xlab("Eye Rear y - Eye Front y") + ylab("Pitch") + 
    theme(axis.text = element_text(family = "Noto Sans", size = 11, colour = "black"),
    axis.title = element_text(family = "Noto Sans", size = 14, colour = "black")) + 
    annotate("text", x = 50, y = -170, label = "Pointing down", size = 4, family = "Noto Sans", colour = "grey20") + 
    annotate("text", x = -65, y = -170, label = "Pointing up", size = 4, family = "Noto Sans", colour = "grey20")

P3 = ggplot(Rooted, aes(x = X_Dir, y = R)) + 
    geom_smooth(method = "lm", colour = "black") + 
    stat_summary_bin(geom = "point", fun = "mean", binwidth = 0.5, 
    size = 3, pch = 21, alpha = 0.7, colour = "black", fill = "lightseagreen") + 
    theme_classic() + xlab("Cyr Base x - Bill Tip x") + ylab("Roll") + 
    theme(axis.text = element_text(family = "Noto Sans", size = 11, colour = "black"),
    axis.title = element_text(family = "Noto Sans", size = 14, colour = "black")) + 
    annotate("text", x = 65, y = -68, label = "Pointing left", size = 4, family = "Noto Sans", colour = "grey20") + 
    annotate("text", x = -55, y = -68, label = "Pointing right", size = 4, family = "Noto Sans", colour = "grey20")

P4 = ggplot(Rooted, aes(x = Y_Dir, y = R)) + 
    geom_smooth(method = "lm", colour = "black") + 
    stat_summary_bin(geom = "point", fun = "mean", binwidth = 0.5, 
    size = 3, pch = 21, alpha = 0.7, colour = "black", fill = "lightseagreen") + 
    theme_classic() + xlab("Eye Rear y - Eye Front y") + ylab("Roll") + 
    theme(axis.text = element_text(family = "Noto Sans", size = 11, colour = "black"),
    axis.title = element_text(family = "Noto Sans", size = 14, colour = "black")) + 
    annotate("text", x = 50, y = -90, label = "Pointing down", size = 4, family = "Noto Sans", colour = "grey20") + 
    annotate("text", x = -65, y = -90, label = "Pointing up", size = 4, family = "Noto Sans", colour = "grey20")

require(grid)
require(gridExtra)
grid.arrange(P1,P2,P3,P4, nrow = 2)

# Again, good. Absolute yaw values are nicely sigmoidal when plotted against the difference in x value between the bill tip and cyr tip. Next, adjusting yaw values such that -90 represents an individual pointing toward the camera, and 90 represents an individual pointing away from the camera (0 = perpendicular)

Rooted_New = Rooted %>% 
    mutate(Y = ifelse(Y > 90, 180 - Y, 
        ifelse(Y < -90, -180 - Y, Y))
        )

# Re-plotting 

Rooted_New %>%  
ggplot(aes(x = abs(X_Dir), y = abs(Y))) + 
    geom_smooth(method = "lm", colour = "black") + 
    stat_summary_bin(geom = "point", fun = "mean", binwidth = 0.5, 
    size = 3, pch = 21, alpha = 0.7, colour = "black", fill = "lightseagreen") + 
    theme_classic() + xlab("|Cyr Base x - Bill Tip x|") + ylab("Yaw") + 
    theme(axis.text = element_text(family = "Noto Sans", size = 11, colour = "black"),
    axis.title = element_text(family = "Noto Sans", size = 14, colour = "black"))

# Moving forward; removing inaccurate surface temperature estimates (i.e. eye region temperature < 20°C or > 45°C, and bill temperatures < 10°C or > 45°C). Checking for outliers after removal of inaccurate temperature estimates. 

Clean = Rooted_New %>% 
    mutate(Eye.Temp = ifelse(Eye.Temp > 20 & Eye.Temp < 45, Eye.Temp, NA),
    Bill.Temp = ifelse(Bill.Temp > 10 & Bill.Temp < 45, Bill.Temp, NA)) %>%    
    mutate(Eye.Temp = ifelse(Eye.Quality == "Good", Eye.Temp, NA),
    Bill.Temp = ifelse(Bill.Quality == "Good", Bill.Temp, NA))

ggplot(Clean, aes(x = 1:nrow(Clean), y = Eye.Temp)) + 
    geom_point(size = 2, pch = 21, colour = "black", fill = "gold2", alpha = 0.5) + 
    annotate("rect", xmin = 1, xmax = nrow(Clean), 
        ymin = mean(Clean$Eye.Temp, na.rm = T) - 4*sd(Clean$Eye.Temp, na.rm = T),
        ymax = mean(Clean$Eye.Temp, na.rm = T) + 4*sd(Clean$Eye.Temp, na.rm = T),
        fill = "grey20", alpha = 0.3, colour = "black"
        ) + theme_classic() + xlab("Row Number") + 
        ylab("Eye Region Temperature (°C)") + 
        annotate("text", x = 3000, y = 30, label = "Grey box represents mean eye region temperature \n +/- 4 times the standard deviation", size = 4, colour = "black")

# One obvious outlier. Observing, removing, and re-plotting

Clean[c(which(Clean$Eye.Temp < 30)),] # TZ24 during control treatment. Note that no orientation was estimated here.
Clean$Eye.Temp[c(which(Clean$Eye.Temp < 30))] = NA

ggplot(Clean, aes(x = 1:nrow(Clean), y = Eye.Temp)) + 
    geom_point(size = 2, pch = 21, colour = "black", fill = "gold2", alpha = 0.5) + 
    annotate("rect", xmin = 1, xmax = nrow(Clean), 
        ymin = mean(Clean$Eye.Temp, na.rm = T) - 4*sd(Clean$Eye.Temp, na.rm = T),
        ymax = mean(Clean$Eye.Temp, na.rm = T) + 4*sd(Clean$Eye.Temp, na.rm = T),
        fill = "grey20", alpha = 0.3, colour = "black"
        ) + theme_classic() + xlab("Row Number") + 
        ylab("Eye Region Temperature (°C)")

# Checking bill temperature estimates 

ggplot(Clean, aes(x = 1:nrow(Clean), y = Bill.Temp)) + 
    geom_point(size = 2, pch = 21, colour = "black", fill = "gold2", alpha = 0.5) + 
    annotate("rect", xmin = 1, xmax = nrow(Clean), 
        ymin = mean(Clean$Bill.Temp, na.rm = T) - 4*sd(Clean$Bill.Temp, na.rm = T),
        ymax = mean(Clean$Bill.Temp, na.rm = T) + 4*sd(Clean$Bill.Temp, na.rm = T),
        fill = "grey20", alpha = 0.3, colour = "black"
        ) + theme_classic() + xlab("Row Number") + 
        ylab("Bill Temperature (°C)") + 
        annotate("text", x = 3000, y = 20, label = "Grey box represents mean bill temperature \n +/- 4 times the standard deviation", size = 4, colour = "black") 

# Great. Averaging rotational and temperature estimates across ten second intervals for analysis.

bin_Fun = function(x) floor(x/10)*10
Clean$Time_Bin = bin_Fun(Clean$Slice)

All_Bin = Clean %>% 
    group_by(ID, Time_Bin, Treatment) %>%
    summarise(Max.Eye = mean(Eye.Temp, na.rm = T),
    Max.Bill = mean(Bill.Temp, na.rm = T),
    Yaw = mean(Y, na.rm = T),
    Pitch = mean(P, na.rm = T),
    Roll = mean(R, na.rm = T),
    X_Dir = mean(X_Dir, na.rm = T)
    ) %>% ungroup() %>%
    mutate(Treatment = factor(Treatment, ordered = T),
    ID = factor(ID),
    Time_Bin = as.numeric(Time_Bin)
    )

############################################################################################
# Assessing effects of yaw on surface temperature values and surface temperature responses #
# to stress exposure                                                                       #
############################################################################################

# Beginning with temperature response at the eye. Note that only temperature measurements from times at which yaw was estimated are used. This is to ensure that sample sizes are equivalent between models including yaw as a predictor, and models excluding yaw as a predictor. This equivalence of sample sizes will permit statistical comparisons of log-likelihoods between models.

Data = All_Bin %>% 
    mutate(Max.Eye = ifelse(!is.na(Max.Eye) & is.na(Yaw), NA, Max.Eye))

# Lastly, ensuring that the individual observed in both treatments is recognised as the same.

Data = Data %>% 
    mutate(Trial = "A") %>% 
    mutate(Trial = ifelse(ID == "TZ24_2", "B", Trial)) %>% 
    mutate(Trial = factor(Trial)) %>% 
    mutate(ID = as.character(ID)) %>% 
    mutate(ID = ifelse(ID == "TZ24_2", "TZ24", ID)) %>% 
    mutate(ID = factor(ID))

Eye_Mod_Bin = bam(Max.Eye ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = Data)

# Checking for residual autocorrelation

acf(resid_gam(Eye_Mod_Bin, incl_na = TRUE))

# Yes, very high residual autocorrelation. Adjusting for this.

Data$start.event = "FALSE"
Data$Bird_Trial = factor(paste(Data$ID, Data$Trial, sep = "_"))
Data$start.event[c(1, c(which(Data$Bird_Trial != lag(Data$Bird_Trial, 1))))] = "TRUE"
Data$start.event = as.logical(Data$start.event)

Rho = acf(residuals(Eye_Mod_Bin), plot = F)$acf[2]
Rho # Quite high: 0.7226316.

Eye_Mod_Bin_Adj = bam(Max.Eye ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = Data$start.event, rho = Rho,
    data = Data)

# Summarising model. Note that residual checks are bypassed because of comprehensive tests in original analyses. 

summary(Eye_Mod_Bin_Adj)

Base_Coefs = Eye_Mod_Bin_Adj$coefficients
Base_Coefs.SE = as.data.frame(sqrt(diag(vcov(Eye_Mod_Bin_Adj, unconditional = TRUE))))
Base_Coefs.SE$Row = c(1:nrow(Base_Coefs.SE))
Base_Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Base_Coefs.SE))))

data.frame("Coefficient" = Base_Coef_Names, "Estimate" = Base_Coefs, "SE" = Base_Coefs.SE[,1]) %>% group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

# Calculating marginal means

{
    Eye_Grid_Before = ref_grid(Eye_Mod_Bin_Adj,  
        at = list(Time_Bin = seq(1, 419, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Data)
    Eye_Grid_After = ref_grid(Eye_Mod_Bin_Adj,  
        at = list(Time_Bin = seq(420, 840, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Data)
}    
    
EMEye_Before = emmeans(Eye_Grid_Before, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Data)
EMEye_After = emmeans(Eye_Grid_After, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Data)

EMEye_Before
EMEye_After

# Plotting marginal means 

emmip(Eye_Mod_Bin_Adj, Treatment ~ Time_Bin, at = list(Time_Bin = seq(0, 840, 10)),
    CIs = TRUE)

# Results largely with previous model (wherein data from all individuals are included). Now adding yaw to model and comparing likelihood and outcomes with above model. Here, 5 knots are assumed in the yaw spline to account for a possible sigmoidal relationship.

Eye_Mod_Bin_Yaw = bam(Max.Eye ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = Data)

# Again, checking residual autocorrelation.

acf(residuals(Eye_Mod_Bin_Yaw, incl_na = T))

# Yes, residual autocorrelation remains but appears lessened. Correcting for this.

Yaw_Rho = acf(residuals(Eye_Mod_Bin_Yaw), plot = F)$acf[2]
Yaw_Rho # Note that rho has indeed decreased: 0.683

Eye_Mod_Bin_Yaw_Adj = bam(Max.Eye ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = Data$start.event, rho = Yaw_Rho,
    data = Data)

# Extracting residuals to assess model assumptions.

Data$Row.ID = c(1:nrow(Data))
All_Bin_Small = Data %>%
    drop_na(Max.Eye, Treatment, Time_Bin, ID, Yaw) %>%
    mutate(Eye_Residuals = resid_gam(Eye_Mod_Bin_Yaw_Adj, incl_na = TRUE)) %>%
    select(Row.ID, Eye_Residuals)
Data = left_join(Data, All_Bin_Small, by = c("Row.ID"))

# First, checking for heteroskedasticity by yaw and by absolute yaw.

ggplot(Data, aes(x = Yaw, y = Eye_Residuals, fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Yaw (Degrees)") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

ggplot(Data, aes(x = abs(Yaw), y = Eye_Residuals, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("| Yaw | (Degrees)") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic() 

# Possible subtle trends, although visually questionable. Formally assessing this tentative heteroskedasticity using two Breusch-Pagan tests.

mod = lm(Eye_Residuals^2 ~ Yaw, data = Data)
ssize = nrow(subset(Data, !is.na(Eye_Residuals) & !is.na(Yaw)))
lambda_val = (ssize - 1) / ssize * var(Data$Eye_Residuals^2, na.rm = T) / 
    (2 * ((ssize - 1) / ssize * var(Data$Eye_Residuals, na.rm = T))^2)
bp_val = summary(mod)$r.squared*(nrow(subset(Data, !is.na(Eye_Residuals) & !is.na(Yaw))))
print(bp_val*lambda_val)
pchisq(bp_val*lambda_val, df = 1, lower.tail = FALSE)
# X2 = 0.6600
# p = 0.4166

Data$Absolute_Yaw = with(Data, abs(Yaw))
mod = lm(Eye_Residuals^2 ~ Absolute_Yaw, data = Data)
ssize = nrow(subset(Data, !is.na(Eye_Residuals) & !is.na(Absolute_Yaw)))
lambda_val = (ssize - 1) / ssize * var(Data$Eye_Residuals^2, na.rm = T) / 
    (2 * ((ssize - 1) / ssize * var(Data$Eye_Residuals, na.rm = T))^2)
bp_val = summary(mod)$r.squared*(nrow(subset(Data, !is.na(Eye_Residuals) & !is.na(Absolute_Yaw))))
print(bp_val*lambda_val)
pchisq(bp_val*lambda_val, df = 1, lower.tail = FALSE)
# X2 = 0.336
# p = 0.562

# No evidence for heteroskedasticity by yaw or absolute yaw (at alpha = 0.05). Assessing residual patterns across other predictors.

p1 = ggplot(Data, aes(x = Time_Bin, y = Eye_Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Time (s)") + ylab("Eye Region Temperature\nResiduals (°C)") + 
    theme_classic()

p2 = ggplot(Data, aes(x = Treatment, y = Eye_Residuals, fill = Treatment)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Treatment Type") + ylab("Eye Region Temperature\nResiduals (°C)") + 
    theme_classic()

pal = colRoz_pal(name = "c.azureus", n = 8, type = "continuous")
p3 = ggplot(Data, aes(x = ID, y = Eye_Residuals, fill = ID)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = pal) + 
    xlab("Individual ID") + ylab("Eye Region Temperature\nResiduals (°C)") + 
    theme_classic()

grid.arrange(p1, p2, p3) 

# Good. Checking normality of residuals, and plotting by fitted values.

p4 = ggplot(Data %>% mutate("Fit" = predict(Eye_Mod_Bin_Yaw_Adj, type = "response", Data)), 
    aes(x = Fit, y = Eye_Residuals, fill = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Fitted Value") + ylab("Eye Region Temperature\nResiduals (°C)") + 
    theme_classic()

p5 = ggplot(Data, aes(x = Eye_Residuals, fill = Treatment)) + 
    geom_histogram(colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Eye Region Temperature\nResiduals (°C)") + ylab("Frequency") + 
    theme_classic()

grid.arrange(p4, p5)

# Good. Note that there is some partitioning of residuals across fitted estimates (likely owing to different average temperatures among individuals), but no correlation between the two. Summarising model.

summary(Eye_Mod_Bin_Yaw_Adj)

# Yaw is significant (at alpha = 0.05). Comparing log likelihoods of model with yaw included as a predictor, and with yaw excluded as a predictor. Note that a log likelihood ratio test is used here.

test_stat = -2*(logLik(Eye_Mod_Bin_Adj) - logLik(Eye_Mod_Bin_Yaw_Adj))
print(test_stat)
# X2 = 10.22988, df = 1.77815
pchisq(as.numeric(test_stat), df = (12.55085 - 10.7727), lower.tail = FALSE)
# p = 0.0046

# Confirming with lmtest

lmtest::lrtest(Eye_Mod_Bin_Adj, Eye_Mod_Bin_Yaw_Adj)
# X2 = 10.23, df = 2.16
# p = 0.006

# Together, strong evidence that inclusion of yaw as a predictor increasing the log likelihood of the model. Plotting marginal means from model.

emmip(Eye_Mod_Bin_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-90, 90, by = 1)),
    CIs = TRUE, data = Data)

# Interesting and expected decline in surface temperature estimates as birds look away from the camera.

emmip(Eye_Mod_Bin_Yaw_Adj, Treatment ~ Time_Bin, 
    at = list(Time_Bin = seq(0, 840, by = 10)),
    CIs = TRUE, data = Data)

# Note that birds in control groups may be slightly hotter at the onset of experimentation than birds in stress-exposure groups. Assessing patterns using case control, after assessing descriptive statistics
 
Coefs = Eye_Mod_Bin_Yaw_Adj$coefficients
Coef.SE = as.data.frame(sqrt(diag(vcov(Eye_Mod_Bin_Yaw_Adj, unconditional = TRUE))))
Coef.SE$Row = c(1:nrow(Coef.SE))
Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Coef.SE))))

data.frame("Coefficient" = Coef_Names, "Estimate" = Coefs, "SE" = Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

{
    Eye_Grid_Before_Yaw = ref_grid(Eye_Mod_Bin_Yaw_Adj,  
        at = list(Time_Bin = seq(1, 419, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Data)
    Eye_Grid_After_Yaw = ref_grid(Eye_Mod_Bin_Yaw_Adj,  
        at = list(Time_Bin = seq(420, 840, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Data)
}    
    
EMEye_Before = emmeans(Eye_Grid_Before_Yaw, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Data)
EMEye_After = emmeans(Eye_Grid_After_Yaw, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Data)

EMEye_Before
EMEye_After # Note eye region temperature significant differs between treatment groups only after onset of handling.

# Calculating marginal means across yaw

Low_Yaw = emmip(Eye_Mod_Bin_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-90, -45, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Low") %>%
    as.data.frame()
Mid_Yaw = emmip(Eye_Mod_Bin_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-44.5, 44.5, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Mid") %>%
    as.data.frame()
High_Yaw = emmip(Eye_Mod_Bin_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(45, 90, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "High") %>%
    as.data.frame()

rbind(Low_Yaw, Mid_Yaw, High_Yaw)

# Point estimates

min(Data$Yaw, na.rm = T); max(Data$Yaw, na.rm = T)

emmip(Eye_Mod_Bin_Yaw_Adj, ~Yaw, at = list(Yaw = min(Data$Yaw, na.rm = T)), type = "response", CIs = TRUE,
    plot = FALSE)
emmip(Eye_Mod_Bin_Yaw_Adj, ~Yaw, at = list(Yaw = max(Data$Yaw, na.rm = T)), type = "response", CIs = TRUE,
    plot = FALSE)

## Running case-control models.
# Individual that was first observed twice must be split into two treatment groups, then re-labelled. 

CControl = Data %>% mutate(Split = paste(ID, Trial, sep = "_")) %>% 
    group_by(Split) %>% 
    mutate("TDiff" = Max.Eye - Max.Eye[1]) %>% 
    ungroup() %>% 
    select(-Split) %>% 
    mutate(ID = factor(ID))

CControl_No_Yaw = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = CControl)

acf(resid_gam(CControl_No_Yaw, incl_na = T))

# Of course, autocorrelation remains.

CControl_Rho = acf(resid_gam(CControl_No_Yaw, incl_na = T), plot = F)$acf[2]
CControl_Rho # 0.7244013

CControl_No_Yaw_Adj = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = CControl$start.event, rho = CControl_Rho,
    data = CControl)

# Quickly assessing residuals, then summarising model.

Diag = CControl %>% 
    drop_na(TDiff, Treatment, Time_Bin, ID) 
Diag$Residuals = resid_gam(CControl_No_Yaw_Adj, incl_na = T)
Diag$Fit = predict(Eye_Mod_Bin_Yaw_Adj, type = "response", newdata = Diag)

p1 = ggplot(Diag, aes(x = Time_Bin, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Time (s)") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

p2 = ggplot(Diag, aes(x = Treatment, y = Residuals, fill = Treatment)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Treatment Type") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

p3 = ggplot(Diag, aes(x = ID, y = Residuals, fill = ID)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = pal) + 
    xlab("Individual ID") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

p4 = ggplot(Diag, aes(x = Fit, y = Residuals, fill = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Fitted Value") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

p5 = ggplot(Diag, aes(x = Residuals, fill = Treatment)) + 
    geom_histogram(colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Eye Region Temperature Residuals (°C)") + ylab("Frequency") + 
    theme_classic()

grid.arrange(p1,p2,p3,p4,p5, nrow = 2)

# Good. Summarising model.

summary(CControl_No_Yaw_Adj)
emmip(CControl_No_Yaw_Adj, Treatment ~ Time_Bin, data = CControl, CIs = TRUE, 
    at = list(Time_Bin = seq(0, 840, by = 10)))
 
CC_Coefs = CControl_No_Yaw_Adj$coefficients
CC_Coef.SE = as.data.frame(sqrt(diag(vcov(CControl_No_Yaw_Adj, unconditional = TRUE))))
CC_Coef.SE$Row = c(1:nrow(CC_Coef.SE))
CC_Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(CC_Coef.SE))))

data.frame("Coefficient" = CC_Coef_Names, "Estimate" = CC_Coefs, "SE" = CC_Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))
    
# Results are consistent with non case-control model. Now including yaw as a predictor.

CControl_Yaw = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = CControl)

Yaw_Rho = acf(residuals(CControl_Yaw), plot = F)$acf[2]
Yaw_Rho # Autocorrelation remains: 0.679

CControl_Yaw_Adj = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = CControl$start.event, rho = Yaw_Rho,
    data = CControl)

# Checking model residuals.

Diag = CControl %>% 
    drop_na(TDiff, Treatment, Time_Bin, ID, Yaw) 
Diag$Residuals = resid_gam(CControl_No_Yaw_Adj, incl_na = T)
Diag$Fit = predict(CControl_No_Yaw_Adj, type = "response", newdata = Diag)

{ 
    p1 = ggplot(Diag, aes(x = Time_Bin, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Time (s)") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

    p2 = ggplot(Diag, aes(x = Yaw, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Yaw (Degrees)") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

    p3 = ggplot(Diag, aes(x = Treatment, y = Residuals, fill = Treatment)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Treatment Type") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

    p4 = ggplot(Diag, aes(x = ID, y = Residuals, fill = ID)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = pal) + 
    xlab("Individual ID") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

    p5 = ggplot(Diag, aes(x = Fit, y = Residuals, fill = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Fitted Value") + ylab("Eye Region Temperature Residuals (°C)") + 
    theme_classic()

    p6 = ggplot(Diag, aes(x = Residuals, fill = Treatment)) + 
    geom_histogram(colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Eye Region Temperature Residuals (°C)") + ylab("Frequency") + 
    theme_classic()
}

grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)

# Residuals well behaved. Summarising and plotting model results.

summary(CControl_Yaw_Adj) # Very similar results to non-case-control model.
 
CC_Coefs_Yaw = CControl_Yaw_Adj$coefficients
CC_Coef_Yaw.SE = as.data.frame(sqrt(diag(vcov(CControl_Yaw_Adj, unconditional = TRUE))))
CC_Coef_Yaw.SE$Row = c(1:nrow(CC_Coef_Yaw.SE))
CC_Coef_Yaw_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(CC_Coef_Yaw.SE))))

data.frame("Coefficient" = CC_Coef_Yaw_Names, "Estimate" = CC_Coefs_Yaw, "SE" = CC_Coef_Yaw.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

{
    Eye_Grid_Before_Yaw_CC = ref_grid(CControl_Yaw_Adj,  
        at = list(Time_Bin = seq(1, 419, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = CControl)
    Eye_Grid_After_Yaw_CC = ref_grid(CControl_Yaw_Adj,  
        at = list(Time_Bin = seq(420, 840, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = CControl)
}    
    
EMEye_Before_CC = emmeans(Eye_Grid_Before_Yaw_CC, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = CControl)
EMEye_After_CC = emmeans(Eye_Grid_After_Yaw_CC, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = CControl)

EMEye_Before_CC
EMEye_After_CC

Low_Yaw_CC = emmip(CControl_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-90, -45, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Low") %>%
    as.data.frame()
Mid_Yaw_CC = emmip(CControl_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-44.5, 44.5, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Mid") %>%
    as.data.frame()
High_Yaw_CC = emmip(CControl_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(45, 90, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "High") %>%
    as.data.frame()

rbind(Low_Yaw_CC, Mid_Yaw_CC, High_Yaw_CC)

min(CControl$Yaw, na.rm = T); max(CControl$Yaw, na.rm = T)

emmip(CControl_Yaw_Adj, ~Yaw, at = list(Yaw = min(CControl$Yaw, na.rm = T)), type = "response", CIs = TRUE,
    plot = FALSE)
emmip(CControl_Yaw_Adj, ~Yaw, at = list(Yaw = max(CControl$Yaw, na.rm = T)), type = "response", CIs = TRUE,
    plot = FALSE)

# Consistent with non case-control model. Plotting for a visual comparison.

emmip(CControl_Yaw_Adj, Treatment ~ Time_Bin, 
    at = list(Time_Bin = seq(0, 840, by = 10)),
    CIs = TRUE, data = CControl, nesting = NULL, nesting.order = FALSE)

emmip(CControl_Yaw_Adj, ~ Yaw, 
    at = list(Yaw = seq(-90, 90, by = 1)),
    CIs = TRUE, data = CControl, nesting = NULL, nesting.order = FALSE)

# Visually consistent with previous findings.
# Producing final plots.

Eye_EMs = emmip(Eye_Mod_Bin_Yaw_Adj, Treatment ~ Time_Bin, CIs = TRUE,
    at = list(Time_Bin = seq(1, 840, by = 1)), plot = FALSE, 
    data = Data) # Note that confidence intervals below are Wald because of uncertainty surrounding simulataneous confidence interval estimates from a weighted model.

Eye_YawAdj_Plot = ggplot(Eye_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    geom_point(data = Data, aes(x = Time_Bin, y = Max.Eye), pch = 21, size = 2, alpha = 0.6) + 
    geom_line() + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + my.theme + 
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_YawAdj_Plot

# Smaller y-limits

Eye_YawAdj_Plot = ggplot(Eye_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    geom_point(data = Data, aes(x = Time_Bin, y = Max.Eye), pch = 21, size = 2, alpha = 0.6) + 
    geom_line() + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + my.theme +
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) +
    scale_y_continuous(limits = c(32.5,37.5)) + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_YawAdj_Plot

# With dots averaged across birds 

Eye_YawAdj_Plot_Av = ggplot(Eye_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "point", fun = "mean", binwidth = 10,
        data = Data, aes(x = Time_Bin, y = Max.Eye), pch = 21, size = 2, alpha = 0.6) + 
    geom_line() + annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 38.75) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Eye Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + my.theme + 
    #scale_y_continuous(limits = c(33,36.5)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_YawAdj_Plot_Av

# And lastly, with CIs around points

Eye_YawAdj_Plot_Av_CIs = ggplot(Eye_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Data, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Data, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 37.85) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Eye Region Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    scale_y_continuous(limits = c(32,38), breaks = c(32, 34, 36, 38)) + my.theme + 
    #scale_y_continuous(limits = c(33,36.5)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_YawAdj_Plot_Av_CIs

ggsave("/home/joshk/Desktop/All/Figures/Eye_Yaw_Corrected_Model.jpeg",  Eye_YawAdj_Plot_Av_CIs,
       height = 7, width = 8, dpi = 800
)

# And with grey boxes signifying duration of handling

Eye_YawAdj_Plot_Av_CIs = ggplot(Eye_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Data, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Data, aes(x = Time_Bin, y = Max.Eye), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 800, y = 37.85) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Eye Region Temperature (°C)") + 
    annotate("rect", xmin = 420, xmax = 840, ymin = -Inf, ymax = Inf, colour = "black", fill = "grey40", alpha = 0.4) + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    scale_y_continuous(limits = c(32,38), breaks = c(32, 34, 36, 38)) + my.theme + 
    #scale_y_continuous(limits = c(33,36.5)) + my.theme + 
    theme(legend.position = c(0.21,0.1), legend.background = NULL)    

Eye_YawAdj_Plot_Av_CIs

ggsave("/home/joshk/Desktop/All/Figures/Eye_Yaw_Corrected_Model_Grey.jpeg",  Eye_YawAdj_Plot_Av_CIs,
       height = 7, width = 8, dpi = 800
)

# Case control plot.

Eye_CC_EMs = emmip(CControl_Yaw_Adj, Treatment ~ Time_Bin, CIs = TRUE,
    at = list(Time_Bin = seq(1, 840, by = 1)), plot = FALSE, 
    nesting = NULL, nesting.order = FALSE,
    data = CControl) 

Eye_CC_CIs = ggplot(Eye_CC_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = CControl, aes(x = Time_Bin, y = TDiff), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = CControl, aes(x = Time_Bin, y = TDiff), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 830, y = 1.45) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Change in Eye Region Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    scale_y_continuous(limits = c(-1.5,1.5), breaks = c(-1.50, -0.75, 0, 0.75, 1.50)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_CC_CIs

ggsave("/home/joshk/Desktop/All/Figures/Case_Control_Eye.jpeg",  Eye_CC_CIs,
       height = 7, width = 8, dpi = 800
)

# And finally, producing yaw plot.

Eye_Yaw_EM = emmip(Eye_Mod_Bin_Yaw_Adj, ~ Yaw, CIs = TRUE,
    at = list(Yaw = seq(-90, 90, by = 1)), plot = FALSE, 
    data = Data)

Eye_Yaw_CIs = ggplot(Eye_Yaw_EM, 
    aes(x = Yaw, y = yvar)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55, fill = "#6f196c", colour = "black") +
    #geom_point(data = Data, aes(x = Yaw, y = Max.Eye), pch = 21, size = 3, alpha = 0.5, colour = "black", fill = "#a42a64") + 
    #stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = CControl, aes(x = Yaw, y = Max.Eye), binwidth = 5, size = 0.5, alpha = 0.8, colour = "black", width = 1) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Data, aes(x = Yaw, y = Max.Eye), binwidth = 1, pch = 21, size = 3, alpha = 0.6, fill = "#6f196c") + 
    annotate(geom = "text", label = "A", size = 10, family = "Noto Sans", x = 85, y = 37.85) + geom_line(colour = "black", linetype = "solid") + 
    theme_bw() + xlab("Yaw (°)") + ylab("Eye Region Temperature (°C)") + 
    geom_vline(xintercept = 0, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(-90, -45, 0, 45, 90),
        labels = c("-90", "-45", "0", "45", "90")) + 
    scale_y_continuous(limits = c(32,38), breaks = c(32, 34, 36, 38)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Eye_Yaw_CIs

ggsave("/home/joshk/Desktop/All/Figures/Eye_Yaw.jpeg",  Eye_Yaw_CIs,
       height = 7, width = 8, dpi = 800
)

####################################################################
# Repeating analyses and model comparison at the level of the bill #
####################################################################

Bill_Data = All_Bin %>% 
    mutate(Max.Bill = ifelse(!is.na(Max.Bill) & is.na(Yaw), NA, Max.Bill))

Bill_Data = Bill_Data %>% 
    mutate(Trial = "A") %>% 
    mutate(Trial = ifelse(ID == "TZ24_2", "B", Trial)) %>% 
    mutate(Trial = factor(Trial)) %>% 
    mutate(ID = as.character(ID)) %>% 
    mutate(ID = ifelse(ID == "TZ24_2", "TZ24", ID)) %>% 
    mutate(ID = factor(ID))

Bill_Mod_Bin = bam(Max.Bill ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"), 
    method = "REML", na.action = na.omit, 
    data = Bill_Data)

acf(resid_gam(Bill_Mod_Bin)) # Severe autocorrelation.

Bill_Data$start.event = "FALSE"
Bill_Data$Bird_Trial = factor(paste(Data$ID, Data$Trial, sep = "_"))
Bill_Data$start.event[c(1, c(which(Data$Bird_Trial != lag(Data$Bird_Trial, 1))))] = "TRUE"
Bill_Data$start.event = as.logical(Bill_Data$start.event)

Bill_Rho = acf(resid_gam(Bill_Mod_Bin), plot = F)$acf[2]
Bill_Rho #0.8648689. Degree of residual autocorrelation quite notable, relative to eye region temperature estimates. Accounting for this in the model below.

Bill_Mod_Bin_Adj = bam(Max.Bill ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = Bill_Data$start.event, rho = Bill_Rho,
    data = Bill_Data)

# Summarising model and model coefficients. Again, mote that residual checks are bypassed because of comprehensive tests in original analyses. 

summary(Bill_Mod_Bin_Adj) 

Coefs = Bill_Mod_Bin_Adj$coefficients
Coef.SE = as.data.frame(sqrt(diag(vcov(Bill_Mod_Bin_Adj, unconditional = TRUE))))
Coef.SE$Row = c(1:nrow(Coef.SE))
Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Coef.SE))))

data.frame("Coefficient" = Coef_Names, "Estimate" = Coefs, "SE" = Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

{
    Bill_Grid_Before = ref_grid(Bill_Mod_Bin_Adj,  
        at = list(Time_Bin = seq(1, 419, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Bill_Data)
    Bill_Grid_After = ref_grid(Bill_Mod_Bin_Adj,  
        at = list(Time_Bin = seq(420, 840, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Bill_Data)
}    
    
EMBill_Before = emmeans(Bill_Grid_Before, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Bill_Data)
EMBill_After = emmeans(Bill_Grid_After, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Bill_Data)

EMBill_Before
EMBill_After

# Quickly plotting

emmip(Bill_Mod_Bin_Adj, Treatment ~ Time_Bin, CIs = TRUE,
    at = list(Time_Bin = seq(1, 800, by = 1)),
    data = Bill_Data)

# All together, results largely mirror those derived from complete data set. Adding yaw in as a predictor.

Bill_Mod_Yaw = bam(Max.Bill ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"), 
    method = "REML", na.action = na.omit, 
    data = Bill_Data)

Bill_Rho_Yaw = acf(resid_gam(Bill_Mod_Yaw), plot = F)$acf[2]
Bill_Rho_Yaw #0.863107. Rho still quite high, even after addition of yaw.

Bill_Mod_Yaw_Adj = bam(Max.Bill ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"), 
    method = "REML", na.action = na.omit, 
    AR.start = Bill_Data$start.event, rho = Bill_Rho_Yaw,
    data = Bill_Data)

# Checking for heteroskedasticity by yaw

Bill_Data$Row.ID = c(1:nrow(Bill_Data))
Resid_extract = Bill_Data %>%
    drop_na(Max.Bill, Treatment, Time_Bin, ID, Yaw) %>%
    mutate(Bill_Residuals = resid_gam(Bill_Mod_Bin_Adj, incl_na = TRUE)) %>%
    select(Row.ID, Bill_Residuals)
Bill_Data = left_join(Bill_Data, Resid_extract, by = c("Row.ID"))

# First, checking for heteroskedasticity by yaw and by absolute yaw.

ggplot(Bill_Data, 
    aes(x = Yaw, y = Bill_Residuals, fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Yaw (Degrees)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

ggplot(Bill_Data, aes(x = abs(Yaw), y = Bill_Residuals, 
    fill = Treatment, linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("| Yaw | (Degrees)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic() 

# Possible subtle trends, although visually questionable. Formally assessing this tentative heteroskedasticity using two Breusch-Pagan tests.

mod = lm(Bill_Residuals^2 ~ Yaw, data = Bill_Data)
ssize = nrow(subset(Bill_Data, !is.na(Bill_Residuals) & !is.na(Yaw)))
lambda_val = (ssize - 1) / ssize * var(Bill_Data$Bill_Residuals^2, na.rm = T) / 
    (2 * ((ssize - 1) / ssize * var(Bill_Data$Bill_Residuals, na.rm = T))^2)
bp_val = summary(mod)$r.squared*(nrow(subset(Bill_Data, !is.na(Bill_Residuals) & !is.na(Yaw))))
print(bp_val*lambda_val)
pchisq(bp_val*lambda_val, df = 1, lower.tail = FALSE)
# X2 = 1.826001
# p = 0.1766006

Bill_Data$Absolute_Yaw = with(Bill_Data, abs(Yaw))
mod = lm(Bill_Residuals^2 ~ Absolute_Yaw, data = Bill_Data)
ssize = nrow(subset(Bill_Data, !is.na(Bill_Residuals) & !is.na(Absolute_Yaw)))
lambda_val = (ssize - 1) / ssize * var(Bill_Data$Bill_Residuals^2, na.rm = T) / 
    (2 * ((ssize - 1) / ssize * var(Bill_Data$Bill_Residuals, na.rm = T))^2)
bp_val = summary(mod)$r.squared*(nrow(subset(Bill_Data, !is.na(Bill_Residuals) & !is.na(Absolute_Yaw))))
print(bp_val*lambda_val)
pchisq(bp_val*lambda_val, df = 1, lower.tail = FALSE)
# X2 = 0.06389092
# p = 0.8004486

# Plotting residual trends across other predictors and fitted variables.

p1 = ggplot(Bill_Data, aes(x = Time_Bin, y = Bill_Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Time (s)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

p2 = ggplot(Bill_Data, aes(x = Treatment, y = Bill_Residuals, fill = Treatment)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Treatment Type") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

p3 = ggplot(Bill_Data, aes(x = ID, y = Bill_Residuals, fill = ID)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = pal) + 
    xlab("Individual ID") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

p4 = ggplot(Bill_Data %>% 
    mutate("Fit" = predict(Eye_Mod_Bin_Yaw_Adj, type = "response", Bill_Data)), 
    aes(x = Fit, y = Bill_Residuals, fill = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Fitted Value") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

p5 = ggplot(Bill_Data, aes(x = Bill_Residuals, fill = Treatment)) + 
    geom_histogram(colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Bill Temperature Residuals (°C)") + ylab("Frequency") + 
    theme_classic()

grid.arrange(p1,p2,p3,p4,p5, nrow = 2)

# Summarising model

summary(Bill_Mod_Yaw_Adj)

Bill_Coefs = Bill_Mod_Yaw_Adj$coefficients
Bill_Coef.SE = as.data.frame(sqrt(diag(vcov(Bill_Mod_Yaw_Adj, unconditional = TRUE))))
Bill_Coef.SE$Row = c(1:nrow(Bill_Coef.SE))
Bill_Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Bill_Coef.SE))))

data.frame("Coefficient" = Bill_Coef_Names, "Estimate" = Bill_Coefs, "SE" = Bill_Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

{
    Bill_Grid_Before_Yaw = ref_grid(Bill_Mod_Yaw_Adj,  
        at = list(Time_Bin = seq(1, 419, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Bill_Data)
    Bill_Grid_After_Yaw = ref_grid(Bill_Mod_Yaw_Adj,  
        at = list(Time_Bin = seq(420, 840, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = Bill_Data)
}    
    
EMBill_Before = emmeans(Bill_Grid_Before_Yaw, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Bill_Data)
EMBill_After = emmeans(Bill_Grid_After_Yaw, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = Bill_Data)

EMBill_Before
EMBill_After

Bill_Low_Yaw = emmip(Bill_Mod_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-90, -45, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Low") %>%
    as.data.frame()
Bill_Mid_Yaw = emmip(Bill_Mod_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-44.5, 44.5, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Mid") %>%
    as.data.frame()
Bill_High_Yaw = emmip(Bill_Mod_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(45, 90, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "High") %>%
    as.data.frame()

rbind(Bill_Low_Yaw, Bill_Mid_Yaw, Bill_High_Yaw)

# Point estimates

min(Bill_Data$Yaw, na.rm = T); max(Bill_Data$Yaw, na.rm = T)

emmip(Bill_Mod_Yaw_Adj, ~Yaw, at = list(Yaw = min(Bill_Data$Yaw, na.rm = T)), type = "response", CIs = TRUE,
    plot = FALSE)
emmip(Bill_Mod_Yaw_Adj, ~Yaw, at = list(Yaw = max(Bill_Data$Yaw, na.rm = T)), type = "response", CIs = TRUE,
    plot = FALSE)

# Interesting. No effect of yaw, however, effect of stress exposure remains. 
# Comparing log likelihood of yaw-adjusted, and non-yaw adjusted models

test_stat = -2*(logLik(Bill_Mod_Bin_Adj) - logLik(Bill_Mod_Yaw_Adj))
print(test_stat)
# X2 = -0.1991817, df = 0.99018
pchisq(as.numeric(test_stat), df = (12.76519 - 11.77496), lower.tail = FALSE)
# p = 1

# Confirming with lmtest

lmtest::lrtest(Bill_Mod_Bin_Adj, Bill_Mod_Yaw_Adj)
# X2 = 0.1992, df = 0.99023
# p = 0.6554

# Plotting results.

emmip(Bill_Mod_Yaw_Adj, Treatment ~ Time_Bin, at = list(Time_Bin = seq(0, 840, by = 10)),
    CIs = TRUE, data = Bill_Data)
emmip(Bill_Mod_Yaw_Adj, ~ Yaw, at = list(Yaw = seq(-90, 90, by = 1)),
    CIs = TRUE, data = Bill_Data)

# Cleaning plots from yaw-adjusted model.

Bill_EMs = emmip(Bill_Mod_Yaw_Adj, Treatment ~ Time_Bin, CIs = TRUE,
    at = list(Time_Bin = seq(1, 840, by = 1)), plot = FALSE, 
    data = Bill_Data) 

Bill_YawAdj_Plot_Av_CIs = ggplot(Bill_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Bill_Data, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Bill_Data, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 38.75) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    #scale_y_continuous(limits = c(27,36), breaks = c(28,30,32,34,36)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL) + 
    geom_line(colour = "black", linetype = "solid") 
    
Bill_YawAdj_Plot_Av_CIs

ggsave("/home/joshk/Desktop/All/Figures/Bill_Yaw_Corrected_Model.jpeg",  Bill_YawAdj_Plot_Av_CIs,
       height = 7, width = 8, dpi = 800
)

# With grey box indicating duration of handling

Bill_YawAdj_Plot_Av_CIs = ggplot(Bill_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Bill_Data, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Bill_Data, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 790, y = 38.75) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    annotate("rect", xmin = 420, xmax = 840, ymin = -Inf, ymax = Inf, colour = "black", fill = "grey40", alpha = 0.4) + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    #scale_y_continuous(limits = c(27,36), breaks = c(28,30,32,34,36)) + my.theme + 
    theme(legend.position = c(0.21,0.1), legend.background = NULL) + 
    geom_line(colour = "black", linetype = "solid") 
    
Bill_YawAdj_Plot_Av_CIs

ggsave("/home/joshk/Desktop/All/Figures/Bill_Yaw_Corrected_Model_Grey.jpeg",  Bill_YawAdj_Plot_Av_CIs,
       height = 7, width = 8, dpi = 800
)

# And yaw from the same model.

Bill_Yaw_EM = emmip(Bill_Mod_Yaw_Adj, ~ Yaw, CIs = TRUE,
    at = list(Yaw = seq(-90, 90, by = 1)), plot = FALSE, 
    data = Bill_Data)

Bill_Yaw_CIs = ggplot(Bill_Yaw_EM, 
    aes(x = Yaw, y = yvar)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55, fill = "#6f196c", colour = "black") +
    #geom_point(data = Bill_Data, aes(x = Yaw, y = Max.Bill), pch = 21, size = 3, alpha = 0.5, colour = "black", fill = "#a42a64") + 
    #stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = CControl, aes(x = Yaw, y = Max.Eye), binwidth = 5, size = 0.5, alpha = 0.8, colour = "black", width = 1) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Bill_Data, aes(x = Yaw, y = Max.Bill), binwidth = 1, pch = 21, size = 3, alpha = 0.6, fill = "#6f196c") + 
    annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 85, y = 38.75) + geom_line(colour = "black", linetype = "solid") + 
    theme_bw() + xlab("Yaw (°)") + ylab("Bill Temperature (°C)") + 
    geom_vline(xintercept = 0, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(-90, -45, 0, 45, 90),
        labels = c("-90", "-45", "0", "45", "90")) + my.theme + 
    #scale_y_continuous(limits = c(31,39), breaks = c(31, 33, 35, 37, 39)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

Bill_Yaw_CIs

ggsave("/home/joshk/Desktop/All/Figures/Bill_Yaw.jpeg",  Bill_Yaw_CIs,
       height = 7, width = 8, dpi = 800
)

# Running case control.

CControl_Bill = Bill_Data %>% mutate(Split = paste(ID, Trial, sep = "_")) %>% 
    group_by(Split) %>% 
    mutate("TDiff" = Max.Bill - Max.Bill[1]) %>% 
    ungroup() %>% 
    select(-Split) %>% 
    mutate(ID = factor(ID))

CControl_Bill_No_Yaw = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = CControl_Bill)

acf(resid_gam(CControl_No_Yaw, incl_na = T))

# Autocorrelation persists. Accounting for this.

CControl_Bill_Rho = acf(resid_gam(CControl_Bill_No_Yaw, incl_na = T), plot = F)$acf[2]
CControl_Bill_Rho # 0.8361072

CControl_Bill_No_Yaw_Adj = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = CControl_Bill$start.event, rho = CControl_Bill_Rho,
    data = CControl_Bill)

# Extracting residuals and assessing distributions.

Diag_Bill = CControl_Bill %>% 
    drop_na(TDiff, Treatment, Time_Bin, ID, Yaw) 
Diag_Bill$Residuals = resid_gam(CControl_Bill_No_Yaw_Adj, incl_na = TRUE)
Diag_Bill$Fit = predict(CControl_Bill_No_Yaw_Adj, type = "response", newdata = Diag_Bill)

{ 
    p1 = ggplot(Diag_Bill, aes(x = Time_Bin, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Time (s)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p2 = ggplot(Diag_Bill, aes(x = Yaw, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Yaw (Degrees)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p3 = ggplot(Diag_Bill, aes(x = Treatment, y = Residuals, fill = Treatment)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Treatment Type") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p4 = ggplot(Diag_Bill, aes(x = ID, y = Residuals, fill = ID)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = pal) + 
    xlab("Individual ID") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p5 = ggplot(Diag_Bill, aes(x = Fit, y = Residuals, fill = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Fitted Value") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p6 = ggplot(Diag_Bill, aes(x = Residuals, fill = Treatment)) + 
    geom_histogram(colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Bill Temperature Residuals (°C)") + ylab("Frequency") + 
    theme_classic()
}

grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)

# Quite homogenous. Summarising, plotting, then adding yaw as a predictor.

summary(CControl_Bill_No_Yaw_Adj)

Bill_CC_Coefs = CControl_Bill_No_Yaw_Adj$coefficients
Bill_CC_Coef.SE = as.data.frame(sqrt(diag(vcov(CControl_Bill_No_Yaw_Adj, unconditional = TRUE))))
Bill_CC_Coef.SE$Row = c(1:nrow(Bill_CC_Coef.SE))
Bill_CC_Coef_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(Bill_CC_Coef.SE))))

data.frame("Coefficient" = Bill_CC_Coef_Names, "Estimate" = Bill_CC_Coefs, "SE" = Bill_CC_Coef.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))
   
emmip(CControl_Bill_No_Yaw_Adj, Treatment ~ Time_Bin,
    at = list(Time_Bin = seq(0, 840, by = 10)), CIs = TRUE,
    data = Diag_Bill)

# Similar to initial models. Adding yaw as a predictor.

CControl_Bill_Yaw = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 3, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    data = CControl_Bill)

acf(resid_gam(CControl_Bill_Yaw, incl_na = T))

# Autocorrelation similar to previous model. Accounting for this.

CControl_Bill_Yaw_Rho = acf(resid_gam(CControl_Bill_Yaw, incl_na = T), plot = F)$acf[2]
CControl_Bill_Yaw_Rho # 0.8358594

CControl_Bill_No_Yaw_Adj = bam(TDiff ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 3, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = CControl_Bill$start.event, rho = CControl_Bill_Yaw_Rho,
    data = CControl_Bill)

# Again, assessing residuals.

Diag_Bill = CControl_Bill %>% 
    drop_na(TDiff, Treatment, Time_Bin, ID, Yaw) 
Diag_Bill$Residuals = resid_gam(CControl_Bill_No_Yaw_Adj, incl_na = TRUE)
Diag_Bill$Fit = predict(CControl_Bill_No_Yaw_Adj, type = "response", newdata = Diag_Bill)

{ 
    p1 = ggplot(Diag_Bill, aes(x = Time_Bin, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Time (s)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p2 = ggplot(Diag_Bill, aes(x = Yaw, y = Residuals, fill = Treatment, 
    linetype = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    geom_smooth(method = "lm", colour = "black") + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Yaw (Degrees)") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p3 = ggplot(Diag_Bill, aes(x = Treatment, y = Residuals, fill = Treatment)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Treatment Type") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p4 = ggplot(Diag_Bill, aes(x = ID, y = Residuals, fill = ID)) + 
    geom_boxplot(colour = "black") + 
    geom_jitter(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = pal) + 
    xlab("Individual ID") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p5 = ggplot(Diag_Bill, aes(x = Fit, y = Residuals, fill = Treatment)) + 
    geom_point(size = 2, pch = 21, colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Fitted Value") + ylab("Bill Temperature Residuals (°C)") + 
    theme_classic()

    p6 = ggplot(Diag_Bill, aes(x = Residuals, fill = Treatment)) + 
    geom_histogram(colour = "black", alpha = 0.7) + 
    scale_fill_manual(values = c("lightseagreen", "grey70")) + 
    xlab("Bill Temperature Residuals (°C)") + ylab("Frequency") + 
    theme_classic()
}

grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)

# Great. Summarising and producing case-control plot.

summary(CControl_Bill_No_Yaw_Adj)

Bill_CC_EMs = emmip(CControl_Bill_No_Yaw_Adj, Treatment ~ Time_Bin, CIs = TRUE,
    at = list(Time_Bin = seq(1, 840, by = 1)), plot = FALSE, 
    nesting = NULL, nesting.order = FALSE,
    data = CControl_Bill) 

Bill_CC_CIs = ggplot(Bill_CC_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = CControl_Bill, aes(x = Time_Bin, y = TDiff), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = CControl_Bill, aes(x = Time_Bin, y = TDiff), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 3.25) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Change in Bill Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + my.theme + 
    #scale_y_continuous(limits = c(-1.5,1.5), breaks = c(-1.50, -0.75, 0, 0.75, 1.50)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL)    

ggsave("/home/joshk/Desktop/All/Figures/Case_Control_Bill.jpeg",  Bill_CC_CIs,
       height = 7, width = 8, dpi = 800
)

# Calculating final descriptive statistics

CC_Coefs_Bill_Yaw = CControl_Bill_No_Yaw_Adj$coefficients
CC_Coef_Bill_Yaw.SE = as.data.frame(sqrt(diag(vcov(CControl_Bill_No_Yaw_Adj, unconditional = TRUE))))
CC_Coef_Bill_Yaw.SE$Row = c(1:nrow(CC_Coef_Bill_Yaw.SE))
CC_Coef_Bill_Yaw_Names = gsub("\\)", "", gsub("\\.[[:digit:]]*.", "", gsub("s\\(", "", rownames(CC_Coef_Bill_Yaw.SE))))

data.frame("Coefficient" = CC_Coef_Bill_Yaw_Names, "Estimate" = CC_Coefs_Bill_Yaw, "SE" = CC_Coef_Bill_Yaw.SE[,1]) %>%
    group_by(Coefficient) %>% summarise("m_Beta" = mean(Estimate), "m_SE" = mean(SE))

{
    Bill_Grid_Before_Yaw_CC = ref_grid(CControl_Bill_No_Yaw_Adj,  
        at = list(Time_Bin = seq(1, 419, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = CControl_Bill)
    Bill_Grid_After_Yaw_CC = ref_grid(CControl_Bill_No_Yaw_Adj,  
        at = list(Time_Bin = seq(420, 840, by = 1)), 
        cov.reduce = FALSE, type = "response",
        data = CControl_Bill)
}    
    
EMBill_Before_CC = emmeans(Bill_Grid_Before_Yaw_CC, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = CControl_Bill)
EMBill_After_CC = emmeans(Bill_Grid_After_Yaw_CC, specs = pairwise ~ Treatment, cov.reduce = mean,
    type = "response", p.adjust.method = "bonferroni",
    data = CControl_Bill)

EMBill_Before_CC
EMBill_After_CC

Low_Yaw_Bill_CC = emmip(CControl_Bill_No_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-90, -45, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Low") %>%
    as.data.frame()
Mid_Yaw_Bill_CC = emmip(CControl_Bill_No_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(-44.5, 44.5, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "Mid") %>%
    as.data.frame()
High_Yaw_Bill_CC = emmip(CControl_Bill_No_Yaw_Adj, ~Yaw, 
    at = list(Yaw = seq(45, 90, by = 1)), type = "response", CIs = TRUE,
    plot = FALSE) %>% summarise("Mean" = mean(yvar), "LCL" = mean(LCL), UCL = mean(UCL)) %>% 
    mutate("Group" = "High") %>%
    as.data.frame()

rbind(Low_Yaw_Bill_CC, Mid_Yaw_Bill_CC, High_Yaw_Bill_CC)

# Lastly, re-running initial bill model (without yaw as a predictor, given the lack of evidence for an influence on log likelihood) with all available data, then plotting.

Bill_Data_Full = All_Bin %>% 
    mutate(Trial = "A") %>% 
    mutate(Trial = ifelse(ID == "TZ24_2", "B", Trial)) %>% 
    mutate(Trial = factor(Trial)) %>% 
    mutate(ID = as.character(ID)) %>% 
    mutate(ID = ifelse(ID == "TZ24_2", "TZ24", ID)) %>% 
    mutate(ID = factor(ID))

Bill_Data_Full$New_ID = paste(Bill_Data_Full$ID, Bill_Data_Full$Trial, sep = "_")
Bill_Data_Full$New_ID = factor(Bill_Data_Full$New_ID)

Bill_Mod_Bin_Full = bam(Max.Bill ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"), 
    method = "REML", na.action = na.omit, 
    data = Bill_Data_Full)

Bill_Data_Full$start.event = "FALSE"
Bill_Data_Full$start.event[c(which(Bill_Data_Full$New_ID != lag(Bill_Data_Full$New_ID, 1)))] = "TRUE"
Bill_Data_Full$start.event[1] = "TRUE"
Bill_Data_Full$start.event = as.logical(Bill_Data_Full$start.event)

Bill_Rho_Full = acf(resid_gam(Bill_Mod_Bin_Full), plot = F)$acf[2]
Bill_Rho_Full # 0.8748577. 

Bill_Mod_Bin_Full_Adj = bam(Max.Bill ~ Treatment + 
    s(Time_Bin, bs = "tp", k = 3) +
    s(Time_Bin, by = Treatment, bs = "tp", k = 3, m = 1) + 
    s(Yaw, k = 5, bs = "tp") + 
    s(ID, bs = "re"),
    method = "REML", na.action = na.omit, 
    AR.start = Bill_Data_Full$start.event, rho = Bill_Rho_Full,
    data = Bill_Data_Full)

summary(Bill_Mod_Bin_Full_Adj)

Bill_Full_EMs = emmip(Bill_Mod_Bin_Full_Adj, Treatment ~ Time_Bin, CIs = TRUE,
    at = list(Time_Bin = seq(1, 840, by = 1)), plot = FALSE, 
    data = Bill_Data_Full) 

Bill_Full = ggplot(Bill_Full_EMs, 
    aes(x = Time_Bin, y = yvar, fill = Treatment, colour = Treatment, linetype = Treatment)) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.55) +
    stat_summary_bin(geom = "errorbar", fun.data = "mean_cl_boot", data = Bill_Data_Full, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, size = 0.5, alpha = 0.8, colour = "black", width = 1, position = position_dodge(width = 25)) + 
    stat_summary_bin(geom = "point", fun = "mean", data = Bill_Data_Full, aes(x = Time_Bin, y = Max.Bill), binwidth = 50, pch = 21, size = 2, alpha = 0.6, position = position_dodge(width = 25)) + 
    geom_line() + annotate(geom = "text", label = "B", size = 10, family = "Noto Sans", x = 830, y = 35.9) + 
    scale_fill_manual(values = c("#f5e364", "#220a4b"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_colour_manual(values = c("grey10", "grey10"), name = "Treatment", labels = c("Control", "Stress-Exposed")) +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Treatment", labels = c("Control", "Stress-Exposed")) + 
    theme_bw() + xlab("Time (s)") + ylab("Bill Temperature (°C)") + 
    geom_vline(xintercept = 420, colour = "grey40", size = 1, linetype = "dashed") + 
    guides(linetype = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_x_continuous(breaks = c(0, 210, 420, 630, 840),
        labels = c("-210", "-105", "0", "105", "210")) + 
    scale_y_continuous(limits = c(27,36), breaks = c(28,30,32,34,36)) + my.theme + 
    theme(legend.position = c(0.18,0.1), legend.background = NULL) + my.theme + 
    geom_line(colour = "black", linetype = "solid")
    
Bill_Full

ggsave("/home/joshk/Desktop/All/Figures/Bill_Full_Model.jpeg",  Bill_YawAdj_Plot_Av_CIs,
       height = 7, width = 8, dpi = 800
)
