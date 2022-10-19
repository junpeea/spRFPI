# Random Forest Prediction Intervals for Spatially Dependent Data

<b>Abstract: </b><br>
Random forest is a popular machine learning technique that is often used on spatial data for prediction. However, in practice, the spatial dependence of a response variable has typically been ignored in constructing prediction intervals with random forest algorithms, which may result in either low accuracy or low efficiency in such applications. We propose a generalized version of out-of-bag guided random forest prediction intervals, which is well-adapted for spatially dependent data. We use dependency-adjusted regression tree (DART) node-splitting and a novel non-parametric kernel out-of-bag estimator to estimate the underlying conditional prediction error distribution. Theoretical results on the asymptotic consistency of our approach are obtained. Empirical simulation studies and the analysis of global earthquake data indicate that our proposed prediction interval provides good coverage, and is generally more efficient than  existing approaches when observations are spatially dependent. 

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
