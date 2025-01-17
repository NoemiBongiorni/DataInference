# -------------------- Statistical Inference Project 23/24 --------------------
# Authors: Eleonora Banterle, Noemi Bongiorni, Matteo Chiesa, Luca Tagliabue

#BODY FAT

# -------------------- Libraries and Dataset Import --------------------

# libraries
library(RColorBrewer); library(corrplot); library(ellipse); library(faraway); library(GGally); library(olsrr); library(leaps); library(boot); library(MASS); library(BAS); library(rgl); library(car); library(Matrix)

# to install missing packages:
install.packages("boot")
install.packages("olsrr")

# import dataset
data <- read.csv("bodyfat_dataset.csv")
data <- data[-182, ] # remove data with 0 fat
data <- data[,-1] # remove Density column
data$Gender <- as.character(data$Gender)

data <- data[-39, ] # remove extremely fat man
data <- data[-42, ] # remove dwarf

# dataset dimensions
dim(data)

# main statistics for each variable
summary(data)

# check if there are any NA values
print(sapply(data, function(x) any(is.na(x))))
# all FALSE, no NAs, otherwise we would have used na.omit(dataset)

# data types
print(sapply(data, typeof))

# preliminary plot to visually inspect the data
# ggpairs(data = data_num, title ="Relationships between predictors & response", lower = list(continuous=wrap("points", alpha = 0.5, size=0.1)))

# colored plot for male/female differentiation
#x11();
ggpairs(data = data, title ="Relationships between predictors & response", lower = list(continuous=wrap("points", alpha = 0.5, size=0.1), combo = "box_no_facet"), aes(col=as.factor(data$Gender)))

# build the linear model not on all data but on a part (train set), test the model on the remaining part (test set)
# select 25 random indices (25 random observations) from the dataset
train = sample(248,190) # about 75% of the data

# Split the dataframe into two for convenience
data_train_cat = data[train,]
data_test_cat = data[-train,]

# get column names with numeric type
data_num_columns <- unlist(lapply(data, is.numeric), use.names = FALSE)
# data_num is the dataset with only numeric covariates
data_num <- data[ , data_num_columns]
data_train <- data_train_cat[ , data_num_columns]
data_test <- data_test_cat[ , data_num_columns]

# column names
column_names <- colnames(data)
# numeric column names
column_num_names <- colnames(data_num)

# -------------------- First Linear Model --------------------
g_0 = lm(class ~ ., data = data_train) # the dot includes all columns from the dataframe
summary(g_0)

# check assumptions of homoscedasticity and normality
plot(g_0$fit, g_0$res, xlab = "Fitted", ylab = "Residuals", 
     main = "Residuals vs Fitted Values", pch = 16)
abline(h = 0, lwd = 2, lty = 2, col = 'red')

qqnorm(g_0$res, ylab = "Raw Residuals", pch = 16)
qqline(g_0$res)
shapiro.test(g_0$res)

# this plot provides a unified view of all methods to analyze and remove outliers and leverage points

# influential plot
par(mfrow = c(1, 1))
influencePlot(g_0, id.method = "identify", main = "Influential Plot", sub = "Circle size is proportional to Cook's Distance")

# residuals vs leverages
plot(g_0, which = 5)

# points with asterisks are likely leverages
influence.measures(g_0)

# another nice plot, alternative to the bubble plot
ols_plot_cooksd_bar(g_0)


# -------------------- Leverages --------------------
X = model.matrix(g_0)
lev = hatvalues(g_0) # returns elements and also the associated row names

p = g_0$rank
n = dim(data_train)[1]

# leverage plot
plot(g_0$fitted.values, lev, ylab = "Leverages", main = "Plot of Leverages", pch = 16, col = 'black')
abline(h = 2 * p/n, lty = 2, col = 'red')
watchout_points_lev = lev[which(lev > 2 * p/n)]
watchout_ids_lev = seq_along(lev)[which(lev > 2 * p/n)]
points(g_0$fitted.values[watchout_ids_lev], watchout_points_lev, col = 'red', pch = 16)

# position of leverages for each pair of covariates
colors = rep('black', nrow(data_train))
colors[watchout_ids_lev] = c('red', 'blue', 'green', 'orange')
pairs(data_train[, c(column_num_names)], pch = 16, col = colors, cex = 1 + 0.5 * as.numeric(colors != 'black'))


# -------------------- Linear Model without Leverages --------------------

# repeat the process in all 3 ways
# make a plot with x = Cook's distance, y = studentized residuals?
# analyze the very small and very large individuals

# Re-generate the linear model after cleaning:
# studentized residuals
# leverages
# Cook's distance

# LEVERAGES
g_no_leverages <- lm(class ~ ., data_train, subset = (lev < 2 * p/n))
summary(g_no_leverages)
AIC(g_no_leverages)

# STUDENTIZED RESIDUALS
stud = rstandard(g_0)

watchout_ids_stud = which(abs(stud) > 2)
watchout_stud = stud[watchout_ids_stud]
watchout_stud

g_no_stud <- lm(class ~ ., data_train, subset = (abs(stud) < 2))
summary(g_no_stud)
AIC(g_no_stud)

# COOK'S DISTANCE
Cdist = cooks.distance(g_0)
watchout_ids_Cdist = which(Cdist > 4/(n - p))
watchout_Cdist = Cdist[watchout_ids_Cdist]
id_to_keep = !(1:n %in% watchout_ids_Cdist)

g_cook = lm(class ~ ., data_train[id_to_keep, ])
summary(g_cook)
AIC(g_cook)


# note that the best model (g_no_stud) is the one without the influential points found by studentized residuals,
# as it has the highest R-squared and the lowest AIC


# -------------------- Normality Assumption Check --------------------

# Residuals vs fitted values
plot(g_no_stud$fit, g_no_stud$res, xlab = "Fitted", ylab = "Residuals", main = "Residuals vs Fitted Values", pch = 16)
abline(h = 0, lwd = 2, lty = 2, col = 'red')

# another plot
plot(g_no_stud, which = 1)

# QQ plot
qqnorm(g_no_stud$res, ylab = "Raw Residuals", pch = 16)
qqline(g_no_stud$res)

# Shapiro-Wilk normality test
shapiro.test(g_no_stud$res)

# reject normality because the p-value is too low, so perform a Box-Cox transformation

# other useful plots
hist(g_no_stud$res, 10, probability = TRUE, col = 'lavender', main = 'residuals')

boxplot(g_no_stud$res, main = "Boxplot of Data Residuals", pch = 16, col = 'lavender')


# -------------------- BOX COX --------------------
b = boxcox(class ~ ., data = data_train)

# x lambda evaluated
best_lambda_ind = which.max(b$y)
best_lambda = b$x[best_lambda_ind]
best_lambda

# lambda = 1 -> no need for Box-Cox transformation, so change approach:
# redo everything excluding the points with high Cook's distance because it's the model with
# the second best R-squared

# -------------------- Recheck Normality Assumption --------------------

# Residuals vs fitted values
plot(g_cook$fit, g_cook$res, xlab = "Fitted", ylab = "Residuals", main = "Residuals vs Fitted Values", pch = 16)
abline(h = 0, lwd = 2, lty = 2, col = 'red')

# another plot
plot(g_cook, which = 1)

# QQ plot
qqnorm(g_cook$res, ylab = "Raw Residuals", pch = 16)
qqline(g_cook$res)

# Shapiro-Wilk normality test
shapiro.test(g_cook$res)

# normality confirmed :)

# other useful plots
hist(g_cook$res, 10, probability = TRUE, col = 'lavender', main = 'residuals')

boxplot(g_cook$res, main = "Boxplot of Data Residuals", pch = 16, col = 'lavender')


# -------------------- Removing Insignificant Covariates --------------------

# Automatic backward selection
AIC(g_cook)
BIC(g_cook)

g_backward_removed = step(g_cook, direction = "backward", trace = T)
summary(g_backward_removed)

