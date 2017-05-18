# Config file for variant prioritisation pipeline

# Config: Gene specifications ---------------------------------------------

# Specify the Gene Symbol, Transcript (refseq), and chromosome number 
# for each run
# I would suggest commenting out but leaving here details of old runs for posterity's sake

# gene.i <- data.frame(name = "PTEN",
#                      transcript = "NM_000314",
#                      chromosome = 10
#                      )

# gene.i <- data.frame(name = "SYNGAP1",
#                      transcript = "NM_006772",
#                      chromosome = 6
#                      )
#


gene.i <- data.frame(name = "DYRK1A",
                     transcript = "NM_001396",
                     chromosome = 21
                     )

# Config: Options ----------------------------------------------------------

opt.annovar.cache = FALSE # If TRUE, check for cached copy of annovar results in temp before rerunning 
opt.generanks.cache = TRUE # If TRUE, check for existing final gene ranking, if false recalculate gene rankings
opt.session.local = FALSE # If TRUE, changes input / output directories 

# If you're running locally (on a desktop or laptop) you need to run R or Rstudio with sudo priviliges
# and make sure you're running in a folder in which you've got read / write privileges
# For me: I run sudo sshfs to mount my Projects folder and run Rstudio on mounted project directory
