# WOWA Term Paper
This is the model I designed for WOWA dataset. To read more about the research context and dataset relevance, please refer to [WOWA website](https://multicast.aspra.uni-bamberg.de/resources/wowa/). 

You can run the whole code on your own computer, in case you have Stan installed and set up, along with [rethinking](https://github.com/rmcelreath/rethinking) package by Richard McElreath used for predictions. Note though that some chunks have certain requirements. For instance, [gmapsdistance](https://cran.r-project.org/web/packages/gmapsdistance/readme/README.html) package makes use of Google Maps and hence requires the respective API being set.

## Reproducing the analysis

To reproduce the analysis, one should open and run the relevant R Markdown notebooks in the following order:

1) *glottoTrees.Rmd* and *walking_distance.Rmd* in any order
2) *pre-processing.Rmd* 
3) *running_stan.Rmd*, which requires *model.stan*, so you might want to inspect the latter first

All the files have used `set.seed()` function to allow the full reproducibility, although one can try to erase (comment out) those lines to check if the model yields similar results and that seed number was not cherry-picked.

## Model fit
To avoid any inconvenience related to that, I have stored all the outputs and relevant files in RDS format. The fitted model can also be downloaded from [this](https://drive.google.com/drive/folders/1It6uNf-HYrlO8n6laIutXjmJBqjZE6WB?usp=sharing) Drive folder. As GitHub does not allow files to be over certain size, Google Drive was used. The rest of the files can be found on this GitHub page. 

## NOTE
The fitting of the model may take up to 8 hours on an average laptop, I suggest either running it on a good hardware or using the RDS file. In other cases, you might want to leave the model being processed during the night or some other time, when you do not need your computer. 
