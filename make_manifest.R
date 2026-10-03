# Generate manifest.json

#install.packages("rsconnect")
library(rsconnect)

rsconnect::writeManifest(
  appDir = ".",
  appPrimaryDoc = "app.R",
  appMode = "shiny"
)

# Connect Cloud requires manifest.json for R deployments, and it must be in the same directory as app.R. Connect Cloud currently uses the R version and package information recorded in that manifest. (docs.posit.co)
# Because dependency detection scans .R files recursively, package calls in phylothin/phylothin.R should be detected as long as that file is included in the app bundle. (rstudio.github.io)

# # Check what will be included:
# rsconnect::listBundleFiles(".")
# 
# # Check detected dependencies:
# rsconnect::appDependencies(
#   appDir = ".",
#   appMode = "shiny"
# )
