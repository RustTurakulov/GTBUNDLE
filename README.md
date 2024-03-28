# Find My SNP -- The Genotyping and Methylation Array Content Search

**Find my SNP** is a Shiny web application designed to provide users with a friendly interface for searching, comparing, and analyzing the content of two Illumina's mainstream genotyping and methylation chips: GDA and EPICv2.

### Deployment and Running Instructions:

Clone the Repository: Begin by cloning the GT CHIPS GitHub repository to your local machine using the following command:

Copy code
```
git clone https://github.com/RustTurakulov/GTBUNDLE.git
```

Navigate to the App Directory: After cloning the repository, navigate to the GTBUNDLE directory, which contains the Shiny web application files.

Install Required Packages: Before running the application, ensure that you have the necessary R packages installed. You can install the required packages by running the following command in your R environment:

```
R
install.packages(c("dplyr", "tidyr", "shiny", "shinyjs", "ggplot2", "gridExtra", "data.table", "DT", "kableExtra", "formattable", "UpSetR", "ComplexHeatmap"))
```
Download large uncompressed csv files with annotation from here

```
http://rcombination.com:8888/findmysnp_supplements/

cyto.csv                                           28-Mar-2024 03:51            78888142
epic.csv                                           28-Mar-2024 03:53           537777815
gdav1.csv                                          28-Mar-2024 03:53           131391524
gdav1cyto.csv                                      28-Mar-2024 03:54           144552374
gsav3.csv                                          28-Mar-2024 03:54            44492080
gsav3as.csv                                        28-Mar-2024 03:54            41204403
gsav3md.csv                                        28-Mar-2024 03:54            43557237
gsav3ps.csv                                        28-Mar-2024 03:54            4735432

```
Those files are too big to put straight on gitHub even compressed. Compressed files makes loading of aplication way too slow even with mutithreaded decompression.  

### Run the Application: 

Once the required packages are installed, you can run the Shiny web application using the following command in your R environment:


```
R
shiny::runApp("GTBUNDLE")
```


Access the Application: After running the above command, the application will start running locally. You can access the application by opening a web browser and navigating to the following URL:

```
Copy code
http://127.0.0.1:XXXX
```

Replace XXXX with the port number displayed in your R console where the application is running.

Interact with the Application: You can now interact with the GT CHIPS application in your web browser. The application provides three main tabs: 'Genes', 'SNPs', and 'Info', each offering different functionalities for searching, comparing, and analyzing the content of Illumina's genotyping and methylation microarrays.

Explore and Analyze: Explore the various features of the application, including searching for gene and SNP IDs, generating bar plots and heatmaps, and accessing summary information about the microarray content.

**Note:** _This Shiny application is provided for generalized informational and demonstration purposes. Users should independently verify the accuracy and reliability of the data presented and exercise their own judgment when interpreting the results._



![Image](https://github.com/users/RustTurakulov/projects/1/assets/72537644/534ca8f4-2a1d-4ee4-8205-dcaededdb96b)


---
## Docker container is available
I have docker container ready to go for this aplication. You have to get `docker` running on local system then just call it:
```
docker run -p 80:3838 --rm --name shiny trust1/shiny:v0.1 sh -c '/init'
```
Then gop to the internet browser and on the page address simply type
```
localhost
```
Docker repository info is here: [trust1/shiny](https://hub.docker.com/repository/docker/trust1/shiny/general)
