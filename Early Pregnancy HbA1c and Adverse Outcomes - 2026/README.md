# Project Description
This project assessed the association between HbA1c measured before 20 weeks of gestation in pregnancy and adverse maternal and newborn outcomes. 

# Study setting 
This analysis used data from the Pregnancy Risk, Infant Surveillance, and Measurement Alliance (PRISMA) Maternal and Newborn Health (MNH) study, a cohort study of pregnant women and their babies from six sites in five countries (s (Kintampo, Ghana; Kisumu, Kenya; Lusaka, Zambia; Karachi, Pakistan; Vellore, southern India; and Hodal, northern India)

# Statistical methods
The independent variable of interest is early pregnancy HbA1c measured before 20 weeks of gestation. The primary outcome is GDM; the secondary outcomes include hypertensive disorder of pregnancy, emergent CS, LGA, and stillbirth after GA 28 weeks. 
We conducted a few analyses
1.	We summarized participant characteristics by GDM status. 
2.	We evaluated the predictive performance of early HbA1c for identifying GDM using receiver operating characteristic curves and area under the curve analyses in both pooled and site-specific datasets. 
3.	We estimated the associations between early HbA1c, as a continuous predictor, and each outcome using three modified Poisson regressions with robust standard errors for each study site separately. A meta-analysis approach was then adopted to pool site-specific estimates to generate overall effect measures.  
  a.	Model 1 only includes HbA1c.  	
  b.	Model 2 includes HbA1c, maternal age, BMI, maternal education, multiparity, and GA at HbA1c measurement.  
  c.	Model 3 includes variables in model 2 except for multiparity and is done among participants with prior pregnancies.
4.	Several sensitivity analyses were performed:
  a.	Associations were assessed by stratifying participants’ anemia status 
  b.	Associations were assessed among those with normal early HbA1c (<5.7%)   
  c.	Associations were assessed excluding those with HbA1c measured after GA 18 weeks 
  d.	Associations were assessed using rule-out and rule-in cutoffs identified from ROC and AUC analysis   
5.	GAM analysis was done to explore the nonlinear relationships between early HbA1c and adverse outcomes.

# R code 
The R code has content table at the beginning listing the analysis for each section, which mainly covers the following:
1.	Descriptive analysis
2.	Primary analysis (meta-analysis to pool site specific estimates from modified Poisson regressions with robust standard errors)
3.	Sensitivity analysis 
4.	GAM analysis 
5.	ROC analysis

# Others 
1.	The study is currently under review in BMC Pregnancy and Childbirth (2026 Sep. 19) 
