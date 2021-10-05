ProteoTracker
================
ProteoTracker is also available directly through this web page: http://ProteoTracker.genexplain.com.

# System Requirements
## Hardware Requirements
The ProteoTracker package requires a standard office computer.


Processing of original ProteoTracker data is instant, however the processing of custome uploaded data is dependnet on the datasets size.


## Software Requirements
### OS Requirements
This package was develop under *Windows 10* running *R version 4.0.2* and was succesfully tested also on *Mac OSX* and *Windows 8*.

### R Dependencies
The ProTargetMiner package depends on the follwing R packages (with version):
```
svglite_1.2.3.2
metap_1.4
cowplot_1.1.0
gridExtra_2.3
broom_0.7.2
Cairo_1.5-12.2
DT_0.16
shinycssloaders_1.0.0
ggrepel_0.8.2
ggalluvial_0.12.2
tidyverse_1.3.0
shinydashboard_0.7.1
shiny_1.5.0

```

# Installation Guide
Please install all necessary packages from CRAN or Bioconductor which should take less than 1 min:
    
    install.packages(c('shiny', 'tidyverse', 'DT', 'shinycssloaders', 'shinydashboard', 'ggalluvial', 'ggrepel', 'Cairo', 'broom', 'gridExtra', 'cowplot', 'metap', 'svglite'))  

Once all packages are installed, please start the ProteoTracker R shiny with the following command:

    shiny::runGitHub('ProteoTracker', 'RZlab')
    
# Description
Detailed characterization of cell type transitions is essential for cell biology in general and particularly for the development of stem cell-based therapies in regenerative medicine. To systematically study such transitions, we introduce a method that simultaneously measures protein expression and thermal stability changes in cells and provide a web-based visualization tool ProteoTracker. We apply our method to study the differences between human pluripotent stem cells (PSC) and several cell types including their parental cell line and differentiated progeny. We detect alterations of protein properties in numerous cellular pathways and components including ribosome biogenesis and demonstrate that modulation of ribosome maturation through SBDS protein can be helpful for manipulating cell stemness in vitro. Using our integrative proteomics approach and the web-based tool, we uncover a molecular basis for the uncoupling of robust transcription from parsimonious translation in stem cells and propose a method for maintaining pluripotency in vitro.

For more information, please refer to the paper: <xx>

<img src='www/workflow.jpg' width='100%' />

# License

This project is covered under the **Apache 2.0 License**.



