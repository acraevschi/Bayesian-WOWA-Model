# WOWA Term Paper
This is the preliminary model I design for WOWA dataset. To read more about the research context and dataset relevance, please refer to [WOWA website](https://multicast.aspra.uni-bamberg.de/resources/wowa/). 

You can run the whole code on your own computer, in case you have Stan installed and set up, along with [rethinking](https://github.com/rmcelreath/rethinking) package by Richard McElreath. Note though that some chunks have certain requirements. For instance, [gmapsdistance](https://cran.r-project.org/web/packages/gmapsdistance/readme/README.html) package makes use of Google Maps and hence requires the respective API being set.

## Model fit
To avoid any inconvenience related to that, I have stored all the outputs and relevant files in RDS format. The fitted model can also be downloaded from [this](https://drive.google.com/drive/folders/1It6uNf-HYrlO8n6laIutXjmJBqjZE6WB?usp=sharing) Drive folder. As GitHub does not allow files to be over certain size, Google Drive was used. The rest of the files can be found on this GitHub page. 

## NOTE
The fitting of the model may take up to 8 hours on an average laptop, I suggest either running it on a good hardware or using the RDS file. In other cases, you might want to leave the model being processed during the night or some other time, when you do not need your computer. 
