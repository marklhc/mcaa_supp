## Drasgow et al (2012) predictive validity coefficients for eight criteria

# Order: Extraversion, Agreeableness,
#        Emotional Stability (Opposite of Neuroticism),
#        Conscientiousness, Openness

## task performance
ta_p <- c(0.02, 0.04, 0.06, 0.07, 0.07)
## contextual performance
c_p <- c(0.15, 0.14, 0.18, 0.20, 0.08)
## counterproductivity
c <- c(-0.04, -0.19, -0.11, -0.23, -0.02)
## training performance
tr_p <- c(0.07, 0.03, 0.10, 0.10, 0.10)
## leadership effectiveness
l_e <- c(0.20, 0.11, 0.15, 0.16, 0.17)
## Turnover
t <- c(-0.08, -0.10, -0.21, -0.17, -0.08)
## adaptability
a <- c(0.17, 0.06, 0.19, 0.14, 0.11)
## Fitness Performance
f_p <- c(0.27, 0.11, 0.13, 0.20, -0.04)

# Combine them together
validity <- rbind(ta_p, c_p, c, tr_p, l_e, t, a, f_p)
# add facet names
colnames(validity) <- c(
    "Extraversion", "Agreeableness",
    "Emotional Stability", "Conscientiousness", "Openness"
)
# Calculate average of these criteria
mean_val <- colMeans(validity)
# Scale latent weight (sum to 1)
latent_weight <- (mean_val)^2 / (sum(mean_val^2))
# Change emotional stability to neuroticism
latent_weight["Emotional Stability"] <- -(latent_weight["Emotional Stability"])
names(latent_weight)[3] <- "Neuroticism"
latent_weight
#  Extraversion     Agreeableness       Neuroticism Conscientiousness          Openness 
#    0.46932640        0.03250183       -0.19509222        0.17949135        0.12358820 
## check sum of absolute value is 1
sum(abs(latent_weight))
