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
Pluripotency is a unique cellular state that holds enormous promise in regenerative medicine but is still not fully understood. To systematically study human cell transitions between the pluripotent stem cell (PSC) and differentiated cell types, we developed a plurifaceted experimental procedure  PISA -Express simultaneously measuring protein expression and thermal stability changes in several consecutive cell types, and a web-based visualization tool ProteoTracker for the analysis. Alteration of protein properties occurring in numerous cellular pathways and components provided clues for relative timing and molecular mechanisms of events in the transitions. Strikingly, we found thermal destabilization of most ribosomal proteins in PSCs compared to all the tested somatic cells and a strongly  reduced expression in one ribosome maturation factor, SBDS. Knock-down of SBDS maintained pluripotency and inhibited spontaneous differentiation of PSCs. Here we demonstrate that obstruction of ribosome maturation can be useful for sustaining cell stemness in vitro, and offer a resource for comprehensive interrogation of proteome alterations during cell transitions.

For more information, please refer to the paper: <xx>

<img src='www/workflow.jpg' width='100%' />

# License

This project is covered under the **Apache 2.0 License**.



