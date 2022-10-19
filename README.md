# Random Forest Prediction Intervals for Spatially Dependent Data

<b>Abstract: </b><br>
The whole world experienced a pandemic due to the outbreak of the COVID-19 dis-ease, caused by a virus called the Severe Acute Respiratory Syndrome Coronavirus2 (SARS-CoV-2). The United States is also severely experienced extensive mortalityfrom this disease. Several studies suggested that many of the pre-existing conditionsincreased the risk of death in those with COVID-19. The long-term exposure to theair pollutants is also proven to be one of the major causes behind this mortality. Inthis study, we have explored the relationship between long-term exposure to one ofthe well-known air pollutants, PM2.5and the mortality from the COVID-19 diseasewhile adjusting for several social and environmental factors. Since the COVID-19pandemic is highly spatial in nature of its association and spread pattern, we believetaking into account the spatial dependence of the disease spread process will be veryuseful. Also, as the disease is spread across time, there is a strong temporal natureof the data. Hence our belief is a spatio-temporal model which can take into accountthe complex spatial and temporal dependencies is essential to understand and inves-tigate the COVID-19. In that respect we are using a Bayesian Zero Inflated NegativeBinomial regression model where the spatial and temporal associations are modeledthrough random effects. Our model can capture the overall uncertainty pattern acrosstime and space. This uncertainty quantification taking into account the spatial corre-lation is an attractive novelty of our approach. We have applied our model on a fewstates where the mortality rates are high compared to the other states of the US. Asa part of this study, we have developed a user friendly R tool which can be used toget inference from any state of the US. Along with that, application can also be car-ried out on the entire US data with aggregated state level information (i.e. work withstate level COVID-19 counts).

<b>Summary Results: </b><br>
![](./Paper_work/Table_1.png) <br>
Table 1 shows summary statistics of our data for the four regions in the USA.

![](./Paper_work/Table_3.png) <br>
Table 3 describes the estimated coefficients and their 95% credible intervals of the covariates on COVID-19 weekly death counts.

![](./Paper_work/Figure_2.png) <br>
Figure 2 shows the time-averaged COVID-19 weekly death probability at risk.


<b>Code: </b><br>
[`Result_general.R`](https://github.com/junpeea/COVID-PM-STZINB/blob/main/Paper_work/Code/Result_general.R) includes the code to provide main tables and figues in our Result session.

[`Result_210401_Div2(NJ,NY,PA).R`](https://github.com/junpeea/COVID-PM-STZINB/blob/main/Papaer_work/Code/Result_210401_Div2(NJ,NY,PA).R) includes the code to provide model outputs in the Mid-Atlantic (New Jersey, New York, and Pennsylvania) study.

[`Result_210401_Div4(IA,KS,MO,NE,ND,SD).R`](https://github.com/junpeea/COVID-PM-STZINB/blob/main/Papaer_work/Code/Result_210401_Div4(IA,KS,MO,NE,ND,SD).R) includes the code to provide main results in the Midwest (Iowa, Kansas, Missouri, Nebraska, North Dakota, and South Dakota) study.

[`Result_210401_Div5(GA,FL,NC,SC).R`](https://github.com/junpeea/COVID-PM-STZINB/blob/main/Papaer_work/Code/Result_210401_Div5(GA,FL,NC,SC).R) includes the code to provide main results in the South Atlantic (Florida, Georgia, North Carolina, and South Carolina) study.

[`Result_210401_Div9(CA,OR,WA).R`](https://github.com/junpeea/COVID-PM-STZINB/blob/main/Papaer_work/Code/Result_210401_Div9(CA,OR,WA).R) includes the code to provide main results in the Pacific (California, Oregon, and Washington) study.

<b>Data: </b><br>

Adjacency.csv: Adjacency information across states and counties in US.

County_details.csv

<b>Reference pages: </b><br>
Reference Dashboard: [Johns Hopkins University COVID-19 dashboard](https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6)<br>
COVID-19 data source 1: [The COVID Tracking Project](https://covidtracking.com/)<br>
COVID-19 data source 2: [2019 Novel Coronavirus COVID-19 (2019-nCoV) Data Repository by Johns Hopkins CSSE ](https://github.com/CSSEGISandData/COVID-19)<br>
PM2.5 data source: [Public available code and data to Reproduce Analyses in <Exposure to air pollution and COVID-19 mortality in the United States>](https://github.com/wxwx1993/PM_COVID) <br>
US Hospitalization: [ArcGIS Hub](https://hub.arcgis.com/search) <br>
