# AtlanticBats_DePaul
#Final Project Atlantic Bats 
#Isabelle Sparreo


#Load Dataset ####
#Import dataset into environment 
attach(ATLANTIC_BATS_SPARREO_MODIFIED)
ProtectedData<-ATLANTIC_BATS_SPARREO_MODIFIED
df <- data.frame(ATLANTIC_BATS_SPARREO_MODIFIED)

#Instal packages ####
install.packages("terra")
install.packages("geodata")
install.packages("predicts")
install.packages("dismo")

library(terra)
library(geodata)
library(predicts)
library(dismo)

require(graphics)
library(graphics)


# Remove NAs from just the Reserve_Area column
reservearea <- na.omit(df$Reserve_Area)

print(reservearea)

df_area <- df[complete.cases(df$Reserve_Area), ]



#Download Bioclim data, pick the correct resolution 10, 5, 2.5, 0.5 ####
bioclim_data <- worldclim_global(var = "bio",
                                 res = 2.5,
                                 path = "data/")



#renaming variables
Lat <-df_area$Latitude
Long<- df_area$Longitude
BioClim1<-df_area$BioClim1
BioClim2<-df_area$BioClim2
BioClim3<-df_area$BioClim3
BioClim4<-df_area$BioClim4
BioClim5<-df_area$BioClim5
BioClim6<-df_area$BioClim6
BioClim7<-df_area$BioClim7
BioClim8<-df_area$BioClim8
BioClim9<-df_area$BioClim9
BioClim10<-df_area$BioClim10
BioClim11<-df_area$BioClim11
BioClim12<-df_area$BioClim12
BioClim13<-df_area$BioClim13
BioClim14<-df_area$BioClim14
BioClim15<-df_area$BioClim15
BioClim16<-df_area$BioClim16
BioClim17<-df_area$BioClim17
BioClim18<-df_area$BioClim18
BioClim19<-df_area$BioClim19



#Map ####
#BioClim3 - Isothermality
ramp <- colorRamp(c("blue","yellow","red"))
plot(bioclim_data[[3]], main="Isothermality", sub="sub-title",
     xlab="Latitude", ylab="Longitude", xlim= c(-62, -30), ylim = c(-35,0), axes =TRUE, 
     col = (rgb(ramp(seq(0, 1, length = 200)), max = 255)))
points(Long, Lat, pch=21, bg= "black", col="white", cex=.7)


install.packages("raster")
install.packages("ggplot2")
install.packages("tmap")  # optional, for better map visualization
install.packages("sf") 
install.packages("sp") 
library(raster)
library(ggplot2)
library(tmap) 
library(sf)
library(sp)




#Shannon Diversity ####
# Step 1: Calculate total number of individuals per site (sum all species columns from 30 to 133)


#Shannon Diversity ####
# Step 1: Calculate total number of individuals per site (sum all species columns from 30 to 133)
df_area$total_individuals <- rowSums(df_area[, 30:133])  # Sum species columns, starting from column 30 to 133

# Step 2: Calculate the proportion of each species at each site
proportions_df <- df_area[, 30:133] / df_area$total_individuals  # Species proportions for each site

# Step 3: Define the Shannon Diversity function
shannon_diversity <- function(proportions) {
  -sum(proportions * log(proportions), na.rm = TRUE)  # Shannon formula, ignoring NA values
}

# Step 4: Apply the Shannon diversity function to each site (row-wise)
shannon_values <- apply(proportions_df, 1, shannon_diversity)

# Step 5: Add the Shannon diversity values to the original data
df_area$Shannon_Diversity <- shannon_values

# View the results
print(df_area)




# Simpson's Diversity Index calculation ####
library(vegan)
total_individuals <- df_area$total_individuals
Count <- df_area$Count

# Function to calculate Simpson's Diversity Index based on summary data
calculate_simpson_index <- function(total_individuals, Count) {
  # Check if the number of individuals is greater than 1
  if (total_individuals <= 1) {
    return(0)  # No diversity if there are fewer than 2 individuals
  }
  D = Count / total_individuals
  return(D)
}

# Apply the function to each site to calculate Simpson's Index
simpson_values <- mapply(calculate_simpson_index, df_area$total_individuals, df_area$Count)

# Store the results in a new data frame
simpson_results <- data.frame(Site = df_area$ID, SimpsonIndex = simpson_values)

# Add the Simpson's Index values as a new column to the data frame
df_area$SimpsonIndex <- simpson_values

# Check the updated data frame
head(df_area)  # View the first few rows to verify



#Eveness ####
# Load necessary libraries
library(vegan)

# Subset the dataframe to include only species columns (30 to 133)
species_df <- df_area[, 30:133]

# Remove rows with no data (if applicable)
species_df_clean <- species_df[rowSums(species_df) > 0, ]

# Calculate Shannon-Wiener diversity index (H') for each sample
shannon_index <- diversity(species_df_clean, index = "shannon")

# Calculate the number of species (S) for each sample
num_species <- specnumber(species_df_clean)

# Calculate Pielou's Evenness index (J') for each sample
pielou_evenness <- shannon_index / log(num_species)

# Combine results into a data frame for each sample
evenness_results <- data.frame(
  Species_Richness = num_species, 
  Pielou_Evenness = pielou_evenness
)

# View the results
print(evenness_results)

# If df still contains the same sample order, merge these results back to df
# Make sure the SampleID matches the row names of df
df_area <- cbind(df_area, evenness_results)

# View the updated df with the new columns
head(df_area)


# Change column 18 (protected- Y or N) to "Unprotected" and "Protected"
df_area <- df_area %>%
  mutate(Protection_Status = ifelse(uc == 0, "unprotected", "protected"))

# Check the result
head(df_area)



#Boxplots for variables by Protection Status####
# Boxplot for Species_Richness
install.packages("ggplot2")
library(ggplot2)

ggplot(df_area, aes(x = Protection_Status, y = Species_Richness, color = Protection_Status)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Species Richness by Protection Status", y = "Species Richness", x = "Protection Status")


# Boxplot for Shannon_Diversity
ggplot(df_area, aes(x = Protection_Status, y = Shannon_Diversity, color = Protection_Status)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Shannon Diversity by Protection Status", y = "Shannon Diversity", x = "Protection Status")


# Boxplot for SimpsonIndex
ggplot(df_area, aes(x = Protection_Status, y = SimpsonIndex, color = Protection_Status)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Simpson Index by Protection Status", y = "Simpson Index", x = "Protection Status")


# Boxplot for Pielou Evenness
ggplot(df_area, aes(x = Protection_Status, y = Pielou_Evenness, color = Protection_Status)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Pielou Evenness by Protection Status", y = "Pielou Evenness", x = "Protection Status")



# Check the distribution of data first (normality test)
shapiro.test(df_area$Species_Richness) # p-value = 2.61e-07 <0.05, not normal
shapiro.test(df_area$Shannon_Diversity) #p-value = 0.3311 >0.05, NORMAL
shapiro.test(df_area$SimpsonIndex) #p-value = 3.639e-11 <0.05, not normal
shapiro.test(df_area$Pielou_Evenness) #p-value = 0.04431 <0.05, not normal
shapiro.test(df_area$Reserve_Area) #p-value < 2.2e-16 <0.05, not normal
 


#Mann-Whitney U Test (if non-normal distribution):####
# Species Richness
wilcox_test_species_richness <- wilcox.test(Species_Richness ~ Protection_Status, data = df_area)
print(wilcox_test_species_richness) #p-value = 0.01749

# Simpson Index
wilcox_test_simpson_index <- wilcox.test(SimpsonIndex ~ Protection_Status, data = df_area)
print(wilcox_test_simpson_index) #p-value = 0.003903

# Pielou Evenness
wilcox_test_pielou_evenness <- wilcox.test(Pielou_Evenness ~ Protection_Status, data = df_area)
print(wilcox_test_pielou_evenness) #p-value = 0.01329



#T-test #####
# Shannon Diversity
t_test_shannon_diversity <- t.test(Shannon_Diversity ~ Protection_Status, data = df_area)
print(t_test_shannon_diversity) #p-value = 0.9717




#Correlations ####
# Load required libraries
install.packages("Hmisc")
library(Hmisc)


# Pearson's correlation for Species_Richness, SimpsonIndex, and Pielou_Evenness
# Pearson correlation for Reserve_Area and Species_Richness
cor_species_richness <- rcorr(df_area$Reserve_Area, df_area$Species_Richness, type="pearson")
cat("Pearson's correlation for Reserve_Area and Species_Richness: ", cor_species_richness$r, "\n")
#Corr = -0.04222033

# Pearson correlation for Reserve_Area and SimpsonIndex
cor_simpson_index <- rcorr(df_area$Reserve_Area, df_area$SimpsonIndex, type="pearson")
cat("Pearson's correlation for Reserve_Area and SimpsonIndex: ", cor_simpson_index$r, "\n")
#Corr = 0.1159787

# Pearson correlation for Reserve_Area and Pielou_Evenness
cor_pielou_evenness <- rcorr(df_area$Reserve_Area, df_area$Pielou_Evenness, type="pearson")
cat("Pearson's correlation for Reserve_Area and Pielou_Evenness: ", cor_pielou_evenness$r, "\n")
#Corr = -0.00402356

#Spearman correlation 
# Spearman's correlation for Reserve_Area and Shannon_Diversity
cor_shannon_diversity <- rcorr(df_area$Reserve_Area, df_area$Shannon_Diversity, type="spearman")
cat("Spearman's correlation for Reserve_Area and Shannon_Diversity: ", cor_shannon_diversity$r, "\n")
#Corr = -0.1617037



#Plot Correlations 
par(mfrow = c(2, 2))
# Plot for Species Richness
plot(df_area$Reserve_Area, df_area$Species_Richness,
     xlab = "Reserve Area (hectares)", ylab = "Species Richness")

# Plot for Shannon Diversity
plot(df_area$Reserve_Area, df_area$Shannon_Diversity,
     xlab = "Reserve Area (hectares)", ylab = "Shannon Diversity")

# Plot for Simpson Index
plot(df_area$Reserve_Area, df_area$SimpsonIndex,
     xlab = "Reserve Area (hectares)", ylab = "Simpson Index")

# Plot for Pielou Evenness
plot(df_area$Reserve_Area, df_area$Pielou_Evenness,
     xlab = "Reserve Area (hectares)", ylab = "Pielou Evenness")



# Scatter plot matrix####
# Visualization using ggpairs, placing 'Reserve_Area' on the x-axis
install.packages("GGally")
library(GGally)
ggpairs(df_area[, c("Reserve_Area", "Species_Richness", "Shannon_Diversity", "SimpsonIndex", "Pielou_Evenness", "Protection_Status")], 
        aes(color = Protection_Status, alpha = 0.5))



#Regression ####
reservearea <- df_area$Reserve_Area
evenness <- df_area$Pielou_Evenness
shannon <- df_area$Shannon_Diversity
simpson <- df_area$SimpsonIndex
richness <- df_area$Species_Richness

par(mfrow = c( 2, 2)) #highlight all four plots and run
plot(reservearea,shannon,
     xlab = "Reserve Area (hectares)", ylab = "Shannon Diversity")
modelShan<-lm(reservearea~shannon)
summary(modelShan) #p-value: 0.9768
abline(lm(shannon~reservearea), col= "blue")

plot(reservearea,evenness,
     xlab = "Reserve Area (hectares)", ylab = "Pielou Evenness")
modelEven<-lm(reservearea~evenness)
summary(modelEven) #p-value: 0.6333
abline(lm(evenness~reservearea), col= "blue")

plot(reservearea,simpson,
     xlab = "Reserve Area (hectares)", ylab = "Simpson Index")
modelSimp<-lm(reservearea~simpson)
summary(modelSimp) #p-value: 0.3598
abline(lm(simpson~reservearea), col= "blue")

plot(reservearea,richness,
     xlab = "Reserve Area (hectares)", ylab = "Species Richness")
modelRich<-lm(reservearea~richness)
summary(modelRich) #p-value: 0.683
abline(lm(richness~reservearea), col= "blue")




#Logarithmic Regressions #####
#If the data is heavily skewed or has a wide range of values, a logarithmic transformation can normalize the distribution. By applying a logarithm to the data, you can compress the scale of large values and spread out the smaller values, making it easier to model.

par(mfrow = c(2, 2))

# Logarithmic Regression for Species Richness
plot(log(df_area$Reserve_Area), log(df_area$Species_Richness),
     xlab = "Log(Reserve Area) (hectares)", ylab = "Log(Species Richness)",
     main = "Log-Log Regression: Species Richness")
modelSpecies <- lm(log(Species_Richness) ~ log(Reserve_Area), data = df_area)
abline(modelSpecies, col = "blue")
summary(modelSpecies) #p-value = 0.3951

# Logarithmic Regression for Shannon Diversity
plot(log(df_area$Reserve_Area), log(df_area$Shannon_Diversity),
     xlab = "Log(Reserve Area) (hectares)", ylab = "Log(Shannon Diversity)",
     main = "Log-Log Regression: Shannon Diversity")
modelShannon <- lm(log(Shannon_Diversity) ~ log(Reserve_Area), data = df_area)
abline(modelShannon, col = "blue")
summary(modelShannon) #p-value = 0.3997

# Logarithmic Regression for Simpson Index
plot(log(df_area$Reserve_Area), log(df_area$SimpsonIndex),
     xlab = "Log(Reserve Area) (hectares)", ylab = "Log(Simpson Index)",
     main = "Log-Log Regression: Simpson Index")
modelSimpson <- lm(log(SimpsonIndex) ~ log(Reserve_Area), data = df_area)
abline(modelSimpson, col = "blue")
summary(modelSimpson) #p-value = 0.3886

# Logarithmic Regression for Pielou Evenness
plot(log(df_area$Reserve_Area), log(df_area$Pielou_Evenness),
     xlab = "Log(Reserve Area) (hectares)", ylab = "Log(Pielou Evenness)",
     main = "Log-Log Regression: Pielou Evenness")
modelEvenness <- lm(log(Pielou_Evenness) ~ log(Reserve_Area), data = df_area)
abline(modelEvenness, col = "blue")
summary(modelEvenness) #p-value = 0.6437





#Discriminate Function Analysis ####
install.packages("dplyr")
install.packages("MASS")
library(dplyr)
library(MASS)
library(ggplot2)


# Ensure that protected_type is a factor
df_area$Protection_Status <- factor(df_area$Protection_Status)

# Perform Discriminant Function Analysis (DFA) with two groups (binary outcome)
dfa_model <- lda(Protection_Status ~ Shannon_Diversity + SimpsonIndex + Count, data = df_area)

# Summary of the DFA model
summary(dfa_model)

# Plot the discriminant analysis results
par(mar = c(5, 5, 2, 2))  # Adjust margins (bottom, left, top, right)
plot(dfa_model)

# Get the discriminant function scores (only Discriminant1 for two groups)
df_area$Discriminant1 <- predict(dfa_model)$x[, 1]


# Use ggplot2 to visualize the results
ggplot(df_area, aes(x = Discriminant1, color = Protection_Status)) +
  geom_density() +  # Using a density plot to show distribution
  theme_minimal() +
  labs(title = "Discriminant Function Analysis (DFA)", x = "Discriminant 1", y = "Density") +
  scale_color_manual(values = c("blue", "red"))  # Customize colors for each group


# Load dplyr package
library(dplyr)

# Perform PCA on the dependent variables ####
colnames(df_area)
pca <- prcomp(df_area[, c("Reserve_Area", "Species_Richness", "Shannon_Diversity", "SimpsonIndex", "Pielou_Evenness")], scale. = TRUE)


# Add the PCA results to the data frame
df_area$pca1 <- pca$x[, 1]
df_area$pca2 <- pca$x[, 2]


# Plot PCA results with ggplot2
library(ggplot2)
ggplot(df_area, aes(x = pca1, y = pca2, color = Protection_Status)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Dependent Variables", x = "Principal Component 1", y = "Principal Component 2") +
  scale_color_manual(values = c("blue", "red"))
#This graph shows how the groups ("Protection_Status") are separated based on the combined information from "Species_Richness", Pielou_Evenness", "Shannon_Diversity", and "SimpsonIndex". The points are colored by the "Protection_Status" variable, which indicates whether the species are protected or unprotected.


#Alternate PCA

#PCA ####
# Load required libraries
install.packages("FactoMineR")
install.packages("factoextra")
library(ggplot2)
library(FactoMineR)
library(factoextra)

# Select the relevant variables for PCA (biodiversity metrics and Reserve_Area)
df_pca_data <- df_area[, c("Species_Richness", "Shannon_Diversity", "SimpsonIndex", "Pielou_Evenness", "Reserve_Area")]

# Check for missing data and remove rows with NAs
df_pca_data_clean <- na.omit(df_pca_data)

# Standardizing the data (important for PCA)
df_pca_data_scaled <- scale(df_pca_data_clean)

# Perform PCA
pca_result <- prcomp(df_pca_data_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA results (explained variance, loadings, etc.)
summary(pca_result)

# Scree plot to show variance explained by each principal component
fviz_eig(pca_result)

# Biplot to visualize the principal components and variable loadings
fviz_pca_biplot(pca_result, geom = c("point", "text"), repel = TRUE)

# Color points by Protection Status
# Make sure Protection_Status is included in your data frame
df_pca_data_clean$Protection_Status <- df_area$Protection_Status[!is.na(df_area$Species_Richness)]

# Biplot with color by Protection Status
fviz_pca_biplot(pca_result, geom = c("point", "text"), habillage = df_pca_data_clean$Protection_Status, repel = TRUE)






