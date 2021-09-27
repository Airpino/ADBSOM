# ADBSOM
Code for the paper draft __Batch Self Organizing maps for distributional data using adaptive distances__  
(In submission)
Authors:  
_Francisco De Carvalho, Antonio Irpino, Rosanna Verde, Antonio Balzanella_

## The draft of the paper  
The draft of the paper is available here https://arxiv.org/abs/1811.06980

## Abstract
The paper deals with a Batch Self Organizing Map algorithm (DBSOM) for data described by distributional-valued variables. This kind of variables is characterized to take as values one-dimensional probability or frequency distributions on a numeric support. The objective function optimized in the algorithm depends on the choice of the distance measure. According to the nature of the date, the $L2$ Wasserstein distance is proposed as one of the most suitable metrics to compare distributions. It is widely used in several contexts of analysis of distributional data. Conventional batch SOM algorithms consider that all variables are equally important for the training of the SOM. However, it is well known that some variables are less relevant than others for this task. In order to take into account the different contribution of the variables we propose an adaptive version of the DBSOM algorithm that tackles this problem with an additional step: a relevance weight is automatically learned for each distributional-valued variable. Moreover, since the L2 Wasserstein distance allows a decomposition into two components: one related to the means and one related to the size and shape of the distributions, also relevance weights are automatically learned for each of the measurement components to emphasize the importance of the different estimated parameters of the distributions. Examples of real and synthetic datasets of distributional data illustrate the usefulness of the proposed DBSOM algorithms.

### Keywords  
Distributional data, clustering, feature weighting, Batch Self Organizing Maps

