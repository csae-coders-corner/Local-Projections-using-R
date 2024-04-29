
![CC Graphics 2024_LocationProjections](https://github.com/csae-coders-corner/Local-Projections-using-R/assets/148211163/ff24a7fb-2c5d-4a66-9267-4fd75089b6e0)

# Local-Projections-using-R

{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Local Projections using R"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "##### By Karthik Narayan"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Local Projections (Jordà 2005, AER) are a convenient and flexible way to estimate dynamic, causal effects of a policy of interest several time periods after it comes into effect. Suppose, we have time series data on a variable $y_{t}$ and want to estimate the dynamic causal effects of a policy $\\varepsilon_{t}$ up to $H$ time periods ahead. A convenient, parsimonious way to estimate these effects is using a Local Projection model of the following form:\n",
    "\n",
    "$$\n",
    "ln(y_{t + h}) = \\alpha^{(h)} + \\beta^{(h)} \\varepsilon_{t} + \\gamma^{(h)} X_{t} + u^{(h)}_{t} \n",
    "$$\n",
    "where $h = 0,1,2,3,...,H$ and $X_{t}$ are a set of control variables (including potentially lagged values of the dependent variable). The causal effects of interest are summarized by the series of coefficients $\\{ \\beta^{(h)} \\}$. These coefficients are typically visualised through an $\\textit{Impulse Response Function}$ plot which plots the horizon on the x-axis and the coefficient corresponding to the horizon on the y-axis along with confidence intervals. If we are instead interested in the cumulative impulse responses, we can replace $ln(y_{t+h})$ with $ln(y_{t+h}) - ln(y_{t-1})$.\n",
    "\n",
    "While Impulse Response Functions are most common in macroeconomic applications using time series data, they are increasingly being used to study the propagation of shocks through time in panel datasets (for a recent example, https://www.nber.org/papers/w31184). This blogpost will provide a brief introduction to estimating and plotting local projection impulse response functions in R."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {
    "vscode": {
     "languageId": "r"
    }
   },
   "outputs": [],
   "source": [
    "# Import necessary packages\n",
    "library(tidyverse)\n",
    "library(data.table)\n",
    "library(sandwich)\n",
    "library(lmtest)\n",
    "library(ggplot2)\n",
    "library(grid)\n",
    "library(gridExtra)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "While there is an R package $\\textit{lpirfs}$ which implements local projections impulse response functions, it lacks the flexibility that one typically desires in specifying standard errors and in plotting styles. Therefore, this blog will guide you through an implementation of local projections-based impulse response functions using just a few standard R packages. Suppose that we have a dataset consisting of time series stored in an object called $data$."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "vscode": {
     "languageId": "r"
    }
   },
   "outputs": [],
   "source": [
    "n_horizon <- 12 # Setting Horizon\n",
    "lpirfresults <- data.frame(rep(0, n_horizon + 1))\n",
    "colnames(lpirfresults) <- c(\"Horizon\")\n",
    "lpirfresults[\"Horizon\"] <- 0:n_horizon\n",
    "lpirfresults[\"IRF\"] <- rep(0, n_horizon + 1)\n",
    "lpirfresults[\"IRFLower90\"] <- rep(0, n_horizon + 1)\n",
    "lpirfresults[\"IRFUpper90\"] <- rep(0, n_horizon + 1)\n",
    "lpirfresults[\"IRFLower68\"] <- rep(0, n_horizon + 1)\n",
    "lpirfresults[\"IRFUpper68\"] <- rep(0, n_horizon + 1)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "In the code snippet above, $\\textit{n\\_horizon}$ sets the horizon $h$ of interest, and $\\textit{lprifs}$ is initialized as the results dataframe which will hold the impulse response coefficients corresponding to each horizon. In addition, columns are created to hold confidence interval bounds at the 90% and 68% level. Next, we write a for-loop which iterates through each horizon of interest, estimates the local projection and stores the coefficients and confidence intervals."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "vscode": {
     "languageId": "r"
    }
   },
   "outputs": [],
   "source": [
    "# Estimating Local Projections\n",
    "for (i in 0:n_horizon){\n",
    "    regressiondata <- cbind(data, shift(data[\"y\"],\n",
    "    n = i, type = \"lead\", give.names = TRUE))\n",
    "    regressiondata <- regressiondata[complete.cases(regressiondata), ]\n",
    "    y <- monthlydata$y\n",
    "    x <- as.matrix(monthlydata[, c(\"Policy\", \"Control_1\", \"Control_2\", \"Control_3\")]) #nolint\n",
    "    localprojection <- lm(y ~ x)\n",
    "    lpirfresults[i + 1, \"IRF\"] <- summary(localprojection)$coefficients[2,1] #nolint\n",
    "    localprojection_robustse <- coeftest(localprojection, vcov = vcovHC(localprojection, type = \"HC1\")) #nolint\n",
    "    lpirfresults[i + 1, \"IRFLower90\"] <- lpirfresults[i + 1, \"IRF\"] - qt(0.95, df = df.residual(localprojection)) * localprojection_robustse[2,2] #nolint\n",
    "    lpirfresults[i + 1, \"IRFUpper90\"] <- lpirfresults[i + 1, \"IRF\"] + qt(0.95, df = df.residual(localprojection)) * localprojection_robustse[2,2] #nolint\n",
    "    lpirfresults[i + 1, \"IRFLower68\"] <- lpirfresults[i + 1, \"IRF\"] - qt(0.84, df = df.residual(localprojection)) * localprojection_robustse[2,2] #nolint\n",
    "    lpirfresults[i + 1, \"IRFUpper68\"] <- lpirfresults[i + 1, \"IRF\"] + qt(0.84, df = df.residual(localprojection)) * localprojection_robustse[2,2] #nolint\n",
    "    regressiondata <- data\n",
    "}\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Here, for illustrative purposes, the explanatory variables include a policy variable and three control variables, which can be modified appropriately based on the available data. Finally, once we have the updated $\\textit{lpirfs}$ dataframe with the results from the regressions, we can plot them in the following way:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "vscode": {
     "languageId": "r"
    }
   },
   "outputs": [],
   "source": [
    "# Plot the IRF\n",
    "irfplot <- ggplot(lpirfresults) +\n",
    "geom_line(aes(x = Horizon,\n",
    "y = IRF)) +\n",
    "geom_ribbon(aes(x = Horizon, ymin = IRFLower90,\n",
    "ymax = IRFUpper90), fill = \"red\", alpha = 0.2) +\n",
    "geom_ribbon(aes(x = Horizon, ymin = IRFLower68,\n",
    "ymax = IRFUpper68), fill = \"red\", alpha = 0.7) +\n",
    "geom_hline(yintercept = 0, color = \"black\", linetype = \"dashed\") +\n",
    "theme_light() +\n",
    "ggtitle(\"Local Projection Impulse Response Function\") +\n",
    "ylab(TeX(\"$log(y{t + h})$\")) +\n",
    "xlab(\"Horizons after Impulse\") +\n",
    "scale_x_continuous(breaks = 0:(nrow(lpirfresults) - 1))\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This code when implemented with actual time series data, will produce the impulse response function from the local projection model of interest. \n",
    "\n",
    "It is important to note here that the above code can be easily modified to accommodate panel data, in which case $\\textit{lm()}$ can be replaced with $\\textit{plm()}$ in the for-loop. \n",
    "\n",
    "Similarly, the for-loop can be modified to consider standard errors constructed in other ways, including through bootstrapping. This can be implemented in a manner identical to how they are constructed for linear regression models."
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "R",
   "language": "R",
   "name": "ir"
  },
  "language_info": {
   "codemirror_mode": "r",
   "file_extension": ".r",
   "mimetype": "text/x-r-source",
   "name": "R",
   "pygments_lexer": "r",
   "version": "4.2.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
