# Masters thesis | WOWA 
The model was designed for the analysis of WOWA dataset. To read more about the research context and dataset relevance, please refer to [WOWA website](https://multicast.aspra.uni-bamberg.de/resources/wowa/). 

You can run the whole code on your own computer, in case you have Stan installed and set up. Note though that some chunks have certain requirements. For instance, [gmapsdistance](https://cran.r-project.org/web/packages/gmapsdistance/readme/README.html) package makes use of Google Maps and hence requires the respective API being set.

## Reproducing the analysis

To reproduce the analysis, one should open and run the relevant R Markdown notebooks in the following order:

1) *glottoTrees.Rmd* and *walking_distance.Rmd* in any order
2) *pre-processing.Rmd* 
3) *run_upd_prior.Rmd*, which requires *model_diff_prior1.stan*, so you might want to inspect the latter first
4) *final_predictions.Rmd* and *making_map.R* involve the stage of predicting the outcome on testing set and inspecting/plotting parameter values

All the files have used `set.seed()` function to allow the full reproducibility, although one can try to erase (comment out) those lines to check if the model yields similar results and that seed number was not cherry-picked.

## Model fit
To avoid any inconvenience related to model fitting, I have stored the fit of the model in an RDS file. The fitted model can thus be downloaded from [this](https://drive.google.com/drive/folders/1It6uNf-HYrlO8n6laIutXjmJBqjZE6WB?usp=sharing) Drive folder. As GitHub does not allow files to be over certain size, Google Drive was used. To reproduce the model comparison, some bits of code (e.g. folders' names) would need to be changed. Alternatively, you could get in contact with me and I will provide the models' fit as an RDS file so you could just compare and reproduce the comparison on your computer. Otherwise, I made sure that the rest of the computations are reproducible by using materials available in this repository. 

## NOTE
The fitting of the model may take up to 8 hours on an average laptop, I suggest either running it on a good hardware or using the RDS file. In other cases, you might want to leave the model being processed during the night or some other time, when you do not need your computer. 
