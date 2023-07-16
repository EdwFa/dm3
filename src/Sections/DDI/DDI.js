import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate, Link } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'
//import "./network.css";

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';
import '../ag-theme-acmecorp.css';


import { variables } from '../Variables.js';

import Slider from 'react-input-slider';


const test_articles = [
  {
    "uid": 34630317,
    "aid": "10.3389/fendo.2021.671225",
    "titl": "The Optimal Time of Ovarian Reserve Recovery After Laparoscopic Unilateral Ovarian Non-Endometriotic Cystectomy.",
    "mesh": "Adult;;; Case-Control Studies;;; Cystectomy/*methods;;; Female;;; Follow-Up Studies;;; Humans;;; Laparoscopy/*methods;;; Ovarian Cysts/pathology/*surgery;;; Ovarian Follicle/*physiology;;; Ovarian Reserve/*physiology;;; Prognosis;;; Prospective Studies;;; *Recovery of Function;;; Young Adult",
    "majr": "",
    "subh": "",
    "auth": "Li, Huaping; Yan, Bin; Wang, Yanli; Shu, Zhiming; Li, Ping; Liu, Yahong; Wang, Ying; Ni, Xiaohong; Liu, Zhou",
    "jour": "Frontiers in endocrinology",
    "affl": "Department of Obstetrics and Gynecology, Shanghai University of Medicine & Health Sciences Affiliated Zhoupu Hospital, Shanghai, China.;;; Department of Obstetrics and Gynecology, Shanghai Punan Hospital, Shanghai, China.;;; Department of Obstetrics and Gynecology, Ren Ji Hospital School of Medicine, Shanghai Jiao Tong University, Shanghai, China.;;; Department of Obstetrics and Gynecology, The First People's Hospital of Zhengzhou, Zhengzhou, China.;;; Department of Obstetrics and Gynecology, Shanghai University of Medicine & Health Sciences Affiliated Zhoupu Hospital, Shanghai, China.;;; Shanghai University of Medicine & Health Sciences, Shanghai, China.;;; Department of Obstetrics and Gynecology, Shanghai University of Medicine & Health Sciences Affiliated Zhoupu Hospital, Shanghai, China.;;; Department of Obstetrics and Gynecology, Shanghai University of Medicine & Health Sciences Affiliated Zhoupu Hospital, Shanghai, China.;;; Department of Obstetrics and Gynecology, Shanghai University of Medicine & Health Sciences Affiliated Zhoupu Hospital, Shanghai, China.;;; Department of Obstetrics and Gynecology, Shanghai University of Medicine & Health Sciences Affiliated Zhoupu Hospital, Shanghai, China.",
    "pdat": "2021",
    "tiab": "BACKGROUND: Laparoscopic ovarian cystectomy is established as the standard surgical approach for the treatment of benign ovarian cysts. However, previous studies have shown that potential fertility can be directly impaired by laparoscopic ovarian cystectomy, diminished ovarian reserve (DOR), and even premature ovarian failure. Therefore, fertility-preserving interventions are required for benign gynecologic diseases. However, there are still little data on the time period required for recovery of ovarian reserve after the laparoscopic unilateral ovarian cystectomy, which is very important for the individualization of treatment protocols. This study aimed at investigating the time needed for the ovarian reserve to recover after laparoscopic unilateral ovarian non-endometriotic cystectomy. MATERIALS AND METHODS: Sixty-seven patients with unilateral ovarian non-endometriotic cyst from Zhoupu and Punan Hospitals who underwent laparoscopic unilateral ovarian cystectomy were recruited as a postoperative observation group (POG). Also, 69 healthy age-matched women without ovarian cyst who did not undergo surgery were recruited as a referent group (RFG). Ovarian reserve with the serum anti-Mullerian hormone (AMH), follicle-stimulating hormone (FSH), estradiol (E2) levels, ovarian arterial resistance index (OARI), and antral follicle counts (AFCs) were measured on the third to fifth days of the same menstrual cycle. A postoperative 6-month follow-up of cases was performed. RESULTS: Compared with RFG, AFC of cyst side in the POG group showed no difference in the first, third, and sixth postoperative month (F = 0.03, F = 0.02, F = 0.55, respectively; p = 0.873, p = 0.878, p = 0.460, respectively). The OARI of cyst side in the POG group revealed no differences in the first, third, and sixth postoperative month (F = 0.73, F = 3.57, F = 1.75, respectively; p = 0.395, p = 0.061, p = 0.701, respectively). In the first month, the postoperative AMH levels significantly declined, reaching 1.88 ng/ml [interquartile range (IQR): 1.61-2.16 ng/ml] in POG and 2.57 ng/ml (IQR: 2.32-2.83 ng/ml) in RFG (F = 13.43, p = 0.000). For the data of AMH levels stratified by age, the same trend was observed between less than 25 and more than 26 years old. At this same time interval, the postoperative rate of decline was significantly lower compared to the preoperative one in POG (32.75%). The same trend was observed between the POG and RFG groups (26.67%). CONCLUSIONS: The optimal time for recovery of ovarian reserve after laparoscopic unilateral ovarian cystectomy is estimated to be 6 months.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/34630317/",
    "urlaid": "https://sci-hub.do/10.3389/fendo.2021.671225",
    "pt": "Journal Article; Research Support, Non-U.S. Gov't",
    "pl": "Switzerland",
    "annotations": null,
    "query_number": 1,
    "score": 0.92,
    "section": "INTRO",
    "text": "Anti-Mullerian hormone (AMH) level testing is a useful screening test for assessing ovarian reserve in women at high risk of diminished ovarian reserve, especially for young women with cancer."
  },
  {
    "uid": 29259490,
    "aid": "RMB212055 10.1002/rmb2.12055",
    "titl": "Age-specific serum anti-Mullerian hormone concentration in Japanese women and its usefulness as a predictor of the ovarian response.",
    "mesh": "без MH",
    "majr": "",
    "subh": "",
    "auth": "Asada, Yoshimasa; Morimoto, Yoshiharu; Nakaoka, Yoshiharu; Yamasaki, Takahiro; Suehiro, Yutaka; Sugimoto, Hikaru; Yoshida, Masayuki; Irahara, Minoru",
    "jour": "Reproductive medicine and biology",
    "affl": "Asada Ladies Clinic Medical Corporation Aichi Japan.;;; HORAC Grand Front Osaka Clinic Osaka Japan.;;; IVF Namba Clini cOsaka Japan.;;; Department of Oncology and Laboratory Medicine Yamaguchi University Graduate School of Medicine Yamaguchi Japan.;;; Department of Oncology and Laboratory Medicine Yamaguchi University Graduate School of Medicine Yamaguchi Japan.;;; Medical Science Department Roche Diagnostics Tokyo Japan.;;; Medical Science Department Roche Diagnostics Tokyo Japan.;;; Department of Obstetrics and Gynecology Graduate School of Biomedical Sciences Tokushima University Tokushima Japan.",
    "pdat": "2017 Oct",
    "tiab": "PURPOSE: To compare the ovarian response predictive ability of anti-Mullerian hormone (AMH), follicle-stimulating hormone (FSH), and estradiol (E2) and to determine the age-specific distribution of serum AMH concentrations of Japanese women. METHODS: This was a multicenter (four-site), observational, analytic, cross-sectional Japanese study consisting of two parts: Study 1 (the prediction of the ovarian response of 236 participants who were undergoing controlled ovarian stimulation [COS]) and Study 2 (the distribution of the AMH concentration with an assay of 417 healthy women who were aged 20-49 years and who had normal menstrual cycles). RESULTS: The AMH had a significantly higher predictive value for the normal and high responders than FSH and E2 as a stronger correlation between the ovarian response and AMH was observed than for FSH and E2. The serum AMH concentration decreased proportionally with age. CONCLUSION: The AMH concentration correlated well with the oocyte count in the patients who were undergoing COS for in vitro fertilization and was shown to predict the risk of ovarian hyperstimulation syndrome among these patients.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/29259490/",
    "urlaid": "https://sci-hub.do/RMB212055 https://sci-hub.do/10.1002/rmb2.12055",
    "pt": "Journal Article",
    "pl": "Japan",
    "annotations": null,
    "query_number": 1,
    "score": 0.92,
    "section": "title",
    "text": "Age-specific serum anti-Mullerian hormone concentration in Japanese women and its usefulness as a predictor of the ovarian response"
  },
  {
    "uid": 33623676,
    "aid": "sfz164 10.1093/ckj/sfz164",
    "titl": "Anti-Mullerian hormone concentrations in women with chronic kidney disease.",
    "mesh": "без MH",
    "majr": "",
    "subh": "",
    "auth": "Wiles, Kate; Anckaert, Ellen; Holden, Francesca; Grace, Jan; Nelson-Piercy, Catherine; Lightstone, Liz; Chappell, Lucy C; Bramham, Kate",
    "jour": "Clinical kidney journal",
    "affl": "Department of Women and Children's Health, King's College London, London, UK.;;; Guy's and St Thomas NHS Foundation Trust, London, UK.;;; Laboratory of Hormonology and Tumour Markers, Universitair Ziekenhuis Brussel, Free University of Brussels, Brussels, Belgium.;;; Department of Renal Medicine, King's College Hospital NHS Foundation Trust, London, UK.;;; Guy's and St Thomas NHS Foundation Trust, London, UK.;;; Guy's and St Thomas NHS Foundation Trust, London, UK.;;; Imperial Healthcare NHS Trust, London, UK.;;; Imperial Healthcare NHS Trust, London, UK.;;; Faculty of Medicine, Centre for Inflammatory Disease, Imperial College London, London, UK.;;; Department of Women and Children's Health, King's College London, London, UK.;;; Guy's and St Thomas NHS Foundation Trust, London, UK.;;; Department of Women and Children's Health, King's College London, London, UK.;;; Department of Renal Medicine, King's College Hospital NHS Foundation Trust, London, UK.",
    "pdat": "2021 Feb",
    "tiab": "BACKGROUND: Serum anti-Mullerian hormone (AMH) is a biomarker of ovarian reserve. There are limited data to guide the clinical interpretation of AMH in women with chronic kidney disease (CKD). The purpose of this study was to examine AMH concentrations in women with CKD compared with women without CKD. METHODS: We conducted a prospective cohort study of serum AMH concentrations in 163 non-pregnant women with CKD. Serum AMH concentrations were compared with age-specific AMH centiles from 887 healthy female controls. RESULTS: Participants included 30 women with Stage 1 CKD, 37 women with Stage 2 CKD, 26 women with Stage 3a CKD, 31 women with Stage 3b CKD and 39 women with Stages 4 and 5 CKD. The median estimated glomerular filtration rate (eGFR) was 51 (interquartile range 31-80) mL/min/1.73 m(2). Serum AMH concentrations were lower in all CKD stages compared with women without CKD. Women ages 20-24 years with CKD had comparable serum AMH concentrations (median 1.959 ng/mL) to women ages 35-39 years without CKD (median 1.995 ng/mL). There was no evidence that eGFR was an independent modifier of serum AMH concentrations. More than half of women with CKD (58%) were predicted to have a low response to gonadotrophin stimulation. CONCLUSIONS: Women with CKD have a lower ovarian reserve and are predicted to have a lower ovarian response to gonadotrophin stimulation compared with women without CKD of a similar age. Women with CKD who fail to conceive within 6 months of regular unprotected intercourse should be considered for fertility assessment and intervention.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/33623676/",
    "urlaid": "https://sci-hub.do/sfz164 https://sci-hub.do/10.1093/ckj/sfz164",
    "pt": "Journal Article",
    "pl": "England",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "title",
    "text": "Anti-Mullerian hormone concentrations in women with chronic kidney disease"
  },
  {
    "uid": 28367303,
    "aid": "10.22074/ijfs.2016.4645",
    "titl": "Follicle Stimulating Hormone and Anti-Mullerian Hormone among Fertile and Infertile Women in Ile-Ife, Nigeria: Is there A Difference?",
    "mesh": "без MH",
    "majr": "",
    "subh": "",
    "auth": "Okunola, Temitope; Olusegun Ajenifuja, Kayode; Morebise Loto, Olabisi; Salawu, Afolabi; Omitinde, Stephen Oluseyi",
    "jour": "International journal of fertility & sterility",
    "affl": "Department of Obstetrics and Gynecology, Obafemi Awolowo University Teaching Hospitals Complex, Ile-Ife, Nigeria.;;; Obafemi Awolowo University, Ile-Ife, Nigeria.;;; Obafemi Awolowo University, Ile-Ife, Nigeria.;;; Ladoke Akintola University of Technology Teaching Hospital, Ogbomoso, Nigeria.;;; Department of Obstetrics and Gynecology, Obafemi Awolowo University Teaching Hospitals Complex, Ile-Ife, Nigeria.",
    "pdat": "2017 Apr-Jun",
    "tiab": "BACKGROUND: Reduced ovarian reserve predicts poor ovarian response and poor suc- cess rates in infertile women who undergo assisted reproductive technology (ART). Ovarian reserve also decreases with age but the rate of decline varies from one woman to another. This study aims to detect differences in ovarian reserve as measured by basal serum follicle stimulating hormone (FSH) and anti-Mullerian hormone (AMH) between a matched cohort of fertile and infertile regularly menstruating women, 18-45 years of age. MATERIALS AND METHODS: This case-control study involved 64 fertile and 64 subfertile women matched by age at recruitment. Peripheral blood samples were taken from the women recruited from the Gynecological and Outpatient Clinics of Obafemi Awolowo University Teaching Hospital, Ile-Ife, Nigeria. Serum FSH and AMH were quantified using ELISA at the Metabolic Research Laboratory of LAUTECH Teaching Hospital, Ogbomoso, Nigeria. RESULTS: A significant difference existed in the mean FSH of fertile (6.97 +/- 3.34) and infertile (13.34 +/- 5.24, P=0.013) women. We observed a significant difference in AMH between fertile (2.71 +/- 1.91) and infertile (1.60 +/- 2.51, P=0.029) women. There was a negative correlation between FSH and AMH in both fertile (r=-0.311, P=0.01) and infertile (r=-0.374, P=0.002) women. CONCLUSION: The difference in ovarian reserve observed in this study suggests that reduced ovarian reserve in regularly menstruating women may be associated with early ovarian ageing or subfertility.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/28367303/",
    "urlaid": "https://sci-hub.do/10.22074/ijfs.2016.4645",
    "pt": "Journal Article",
    "pl": "Iran",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "INTRO",
    "text": "Longitudinal studies in fertile women have shown declines in anti-Mullerian hormone (AMH) levels with age; it is the earliest marker of decline in ovarian reserve in young women."
  },
  {
    "uid": 34044635,
    "aid": "10.1177_03000605211016161 10.1177/03000605211016161",
    "titl": "Factors affecting the accuracy and reliability of the measurement of anti-Mullerian hormone concentration in the clinic.",
    "mesh": "*Anti-Mullerian Hormone;;; Enzyme-Linked Immunosorbent Assay;;; Gonadotropin-Releasing Hormone;;; Humans;;; *Ovarian Reserve;;; Reproducibility of Results",
    "majr": "",
    "subh": "",
    "auth": "Fu, Yun-Xing; Wang, Hui; Hu, Ting; Wang, Fei-Miao; Hu, Rong",
    "jour": "The Journal of international medical research",
    "affl": "Ningxia Medical University, General Hospital of Ningxia Medical University, Yinchuan, Ningxia, China.;;; Reproductive Medicine Center, Yinchuan Maternal and Child Health Hospital, Yinchuan, Ningxia, China.;;; Gansu Province Maternity and Child-care hospital, Lan zhou, Gansu, China.;;; Reproductive Medicine Center, Key Laboratory of Fertility Preservation and Maintenance of Ministry of Education, Ningxia Medical University, General Hospital of Ningxia Medical University, Yinchuan, Ningxia, China.;;; Reproductive Medicine Center, Key Laboratory of Fertility Preservation and Maintenance of Ministry of Education, Ningxia Medical University, General Hospital of Ningxia Medical University, Yinchuan, Ningxia, China.",
    "pdat": "2021 May",
    "tiab": "OBJECTIVE: We aimed to identify the factors that influence serum anti-Mullerian hormone (AMH) concentration measurements. METHODS: We collected serum samples between May and September 2018 and compared the effect on AMH concentration measured by ELISA of conditions including venepuncture, storage time, storage temperature, locations of the reaction microplate, and the use of the oral contraceptive pill and gonadotrophin-releasing hormone (GnRH). RESULTS: AMH concentration was not affected by food intake but was affected by haemolysis. It was also much higher in samples on the edge of the ELISA microtitre plate. AMH concentration increased after incubation at room temperature for 1 day, 4 degrees C for 3 days, -20 degrees C for 1 month and -40 degrees C for 4 months, but no change occurred during storage at -80 degrees C for 9 months. AMH concentration was high in patients following GnRH agonist treatment but was not affected by oral contraceptives. CONCLUSIONS: No fasting is required prior to AMH measurement. Placement of serum samples on the edge of microtitre plates affects the results of the AMH ELISA. If serum samples cannot be assayed immediately, it is best to store them at -80 degrees C. Basal AMH concentration cannot be used as a measure of ovarian reserve after GnRH agonist treatment.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/34044635/",
    "urlaid": "https://sci-hub.do/10.1177_03000605211016161 https://sci-hub.do/10.1177/03000605211016161",
    "pt": "Journal Article",
    "pl": "England",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "title",
    "text": "Factors affecting the accuracy and reliability of the measurement of anti-Mullerian hormone concentration in the clinic"
  },
  {
    "uid": 30396559,
    "aid": "S0015-0282(18)30547-8 10.1016/j.fertnstert.2018.06.037",
    "titl": "Dietary factors and serum antimullerian hormone concentrations in late premenopausal women.",
    "mesh": "Adult;;; Anti-Mullerian Hormone/*blood;;; Biomarkers/blood;;; Case-Control Studies;;; Cohort Studies;;; Cross-Sectional Studies;;; Dietary Fats/*administration & dosage/*blood;;; Energy Intake/physiology;;; Feeding Behavior/*physiology;;; Female;;; Humans;;; Middle Aged;;; *Nuns;;; Premenopause/*blood",
    "majr": "",
    "subh": "",
    "auth": "Anderson, Chelsea; Mark Park, Yong-Moon; Stanczyk, Frank Z; Sandler, Dale P; Nichols, Hazel B",
    "jour": "Fertility and sterility",
    "affl": "Department of Epidemiology, University of North Carolina, Chapel Hill, North Carolina.;;; Epidemiology Branch, National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina.;;; Department of Obstetrics and Gynecology and Department of Preventive Medicine, University of Southern California Keck School of Medicine, Los Angeles, California.;;; Epidemiology Branch, National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina.;;; Department of Epidemiology, University of North Carolina, Chapel Hill, North Carolina. Electronic address: hazel.nichols@unc.edu.",
    "pdat": "2018 Nov",
    "tiab": "OBJECTIVE: To study the associations between dietary factors and circulating antimullerian hormone (AMH) concentrations among late premenopausal women. DESIGN: AMH concentrations were measured in serum samples collected at enrollment from 296 women (aged 35-45 years) in the Sister Study cohort. Usual dietary intakes in the past 12 months were assessed using a validated food frequency questionnaire. Dietary exposures of interest included macronutrients, dietary fat subtypes, fiber, and glycemic index. Multivariable linear regression was used to evaluate associations between dietary variables and serum AMH concentrations. We also used nutrient density models to examine isocaloric replacement of macronutrients. SETTING: Not applicable. PATIENTS: Women aged 35-45 years. INTERVENTIONS: Not applicable. MAIN OUTCOME MEASURES: Serum AMH concentrations in nanograms per milliliter (ng/mL). RESULTS: AMH concentrations were positively associated with percentage of energy from carbohydrates (beta per 5% calories = 0.141 [95% CI 0.023, 0.259]; P trend = .019), and inversely associated with percentage of energy from fat (beta per 5% calories = -0.152 [95% CI -0.299, -0.004]; P trend = .044). In analyses of dietary fat subtypes, AMH decreased with increasing monounsaturated fatty acids (P trend = .082) and polyunsaturated fatty acids (P trend = .043), particularly omega-6 fatty acids (P trend = .044), whereas no strong trend was observed for saturated fatty acids. Protein and alcohol intake were not strongly associated with AMH. CONCLUSIONS: Our cross-sectional analyses in a sample of late premenopausal women suggest that dietary fat intake may be inversely associated with circulating AMH concentrations. Further research in prospective studies is warranted to evaluate dietary factors as potential modifiers of ovarian reserve.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/30396559/",
    "urlaid": "https://sci-hub.do/S0015-0282(18)30547-8 https://sci-hub.do/10.1016/j.fertnstert.2018.06.037",
    "pt": "Journal Article; Research Support, N.I.H., Extramural; Research Support, N.I.H., Intramural; Research Support, Non-U.S. Gov't",
    "pl": "United States",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "title",
    "text": "Dietary factors and serum Anti-Mullerian hormone concentrations in late premenopausal women"
  },
  {
    "uid": 24140311,
    "aid": "S1472-6483(13)00452-5 10.1016/j.rbmo.2013.07.014",
    "titl": "Fertility preservation in patients with haematological disorders: a retrospective cohort study.",
    "mesh": "Cohort Studies;;; Drug-Related Side Effects and Adverse Reactions/metabolism;;; Estradiol/metabolism;;; Female;;; Fertility Preservation/*methods/*statistics & numerical data;;; Gonadotropins/administration & dosage/pharmacology;;; Hematologic Diseases/*physiopathology/radiotherapy;;; Humans;;; Ovarian Follicle/drug effects;;; Regression Analysis;;; *Reproductive Techniques, Assisted;;; Retrospective Studies",
    "majr": "",
    "subh": "",
    "auth": "Senapati, Suneeta; Morse, Christopher B; Sammel, Mary D; Kim, Jayeon; Mersereau, Jennifer E; Efymow, Brenda; Gracia, Clarisa R",
    "jour": "Reproductive biomedicine online",
    "affl": "Department of Obstetrics and Gynecology, University of Pennsylvania, United States. Electronic address: Suneeta.Senapati@uphs.upenn.edu.;;; Department of Obstetrics and Gynecology, University of Pennsylvania, United States.;;; Department of Obstetrics and Gynecology, University of Pennsylvania, United States.;;; Department of Obstetrics and Gynecology, University of North Carolina, Chapel Hill, United States.;;; Department of Obstetrics and Gynecology, University of North Carolina, Chapel Hill, United States.;;; Department of Obstetrics and Gynecology, University of Pennsylvania, United States.;;; Department of Obstetrics and Gynecology, University of Pennsylvania, United States.",
    "pdat": "2014 Jan",
    "tiab": "This study investigated the factors associated with utilization of fertility preservation and the differences in treatments and outcomes by prior chemotherapy exposure in patients with haematological diseases. This study included all 67 women with haematological diseases seen for fertility preservation consultation at two university hospitals between 2006 and 2011. Of the total, 49% had lymphoma, 33% had leukaemia, 7% had myelodysplastic syndrome and 4% had aplastic anaemia; 46% had prior chemotherapy; and 33% were planning for bone marrow transplantation, 33% pursued ovarian stimulation and 7% used ovarian tissue banking; and 48% of patients did not pursue fertility preservation treatment. All five cycle cancellations were in the post-chemotherapy group: three patients with leukaemia and two with lymphoma. Patients with prior chemotherapy had lower baseline antral follicle count (10 versus 22) and received more gonadotrophins to achieve similar peak oestradiol concentrations, with no difference in oocyte yield (10.5 versus 10) after adjustment for age. Embryo yield was similar between those who had prior chemotherapy and those who had not. Half of the patients with haematological diseases who present for fertility preservation have been exposed to chemotherapy. While ovarian reserve is likely impaired in this group, oocyte yield may be acceptable.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/24140311/",
    "urlaid": "https://sci-hub.do/S1472-6483(13)00452-5 https://sci-hub.do/10.1016/j.rbmo.2013.07.014",
    "pt": "Journal Article",
    "pl": "Netherlands",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "DISCUSS",
    "text": "Serum anti-Mullerian hormone concentrations may be helpful in assessing ovarian reserve in this patient population in order to determine the optimal stimulation protocol, but these data are not available for the current study."
  },
  {
    "uid": 24396344,
    "aid": "10.1155/2013/125080",
    "titl": "Anti-mullerian hormone as a sensitive marker of ovarian function in young cancer survivors.",
    "mesh": "без MH",
    "majr": "",
    "subh": "",
    "auth": "Krawczuk-Rybak, Maryna; Leszczynska, Elzbieta; Poznanska, Marta; Zelazowska-Rutkowska, Beata; Wysocka, Jolanta",
    "jour": "International journal of endocrinology",
    "affl": "Department of Pediatric Oncology and Hematology, Medical University of Bialystok, Waszyngtona 17, 15-274 Bialystok, Poland.;;; Department of Pediatric Oncology and Hematology, Medical University of Bialystok, Waszyngtona 17, 15-274 Bialystok, Poland.;;; Department of Pediatric Oncology and Hematology, Medical University of Bialystok, Waszyngtona 17, 15-274 Bialystok, Poland.;;; Department of Pediatric Laboratory Diagnostics, Medical University of Bialystok, Poland.;;; Department of Pediatric Laboratory Diagnostics, Medical University of Bialystok, Poland.",
    "pdat": "2013",
    "tiab": "We evaluated ovarian function by measuring the levels of anti-Mullerian hormone (AMH), estradiol, and gonadotropins in 83 young women treated for cancer during childhood and adolescence, and classified according to post-treatment gonadal toxicity versus 38 healthy females. Results. The mean AMH values were lower in the entire cohort independently of the risk group as compared to the control, whereas FSH was elevated only in the high risk group. The lowest AMH values were noted in patients after bone marrow transplantation (BMT) and those treated for Hodgkin lymphoma (HL). Nineteen patients (22.9%) had elevated FSH. They all had low AMH values. Lowered AMH values (but with normal FSH and LH) were observed in 43 patients (51.8%). There was no effect of age at the time of treatment (before puberty, during or after puberty) on AMH levels. Conclusion. Our results show the utility of AMH measurement as a sensitive marker of a reduced ovarian reserve in young cancer survivors. Patients after BMT and patients treated for HL, independently of age at treatment (prepuberty or puberty), are at the highest risk of gonadal damage and early menopause.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/24396344/",
    "urlaid": "https://sci-hub.do/10.1155/2013/125080",
    "pt": "Journal Article",
    "pl": "Egypt",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "abstract",
    "text": "We evaluated ovarian function by measuring the levels of anti-Mullerian hormone (AMH), estradiol, and gonadotropins in 83 young women treated for cancer during childhood and adolescence, and classified according to post-treatment gonadal toxicity versus 38 healthy females."
  },
  {
    "uid": 34375942,
    "aid": "S0160-4120(21)00434-7 10.1016/j.envint.2021.106809",
    "titl": "Urinary phthalate metabolite concentrations are negatively associated with follicular fluid anti-mullerian hormone concentrations in women undergoing fertility treatment.",
    "mesh": "Adult;;; *Anti-Mullerian Hormone;;; Cross-Sectional Studies;;; Female;;; Follicular Fluid;;; Humans;;; *Phthalic Acids",
    "majr": "",
    "subh": "",
    "auth": "Sacha, Caitlin R; Souter, Irene; Williams, Paige L; Chavarro, Jorge E; Ford, Jennifer; Mahalingaiah, Shruthi; Donahoe, Patricia K; Hauser, Russ; Pepin, David; Minguez-Alarcon, Lidia",
    "jour": "Environment international",
    "affl": "Massachusetts General Hospital Fertility Center, Boston, MA, United States; Massachusetts General Hospital Pediatric Surgical Research Laboratories, Boston, MA, United States; Harvard Medical School, Boston, MA, United States. Electronic address: CSACHA@partners.org.;;; Massachusetts General Hospital Fertility Center, Boston, MA, United States; Harvard Medical School, Boston, MA, United States.;;; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA, United States.;;; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA, United States; Department of Nutrition and Department of Epidemiology, Harvard T.H. Chan School of Public Health, Boston, MA, United States; Channing Division of Network Medicine, Harvard Medical School & Brigham and Women's Hospital, United States.;;; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA, United States.;;; Massachusetts General Hospital Fertility Center, Boston, MA, United States; Harvard Medical School, Boston, MA, United States; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA, United States.;;; Massachusetts General Hospital Pediatric Surgical Research Laboratories, Boston, MA, United States; Harvard Medical School, Boston, MA, United States.;;; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA, United States.;;; Massachusetts General Hospital Pediatric Surgical Research Laboratories, Boston, MA, United States; Harvard Medical School, Boston, MA, United States.;;; Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA, United States; Channing Division of Network Medicine, Harvard Medical School & Brigham and Women's Hospital, United States.",
    "pdat": "2021 Dec",
    "tiab": "Exposure to phthalates, endocrine-disrupting chemicals commonly used as plasticizers and in consumer products, has been associated with infertility and premature ovarian failure. Our objective was to investigate whether urinary phthalate metabolite concentrations were associated with pre-ovulatory follicular fluid (FF) anti-mullerian hormone (AMH) concentrations in women undergoing fertility treatment. This cross-sectional analysis included 138 women with urinary phthalate data available in the Environment and Reproductive Health (EARTH) Study (2010-2016) in whom FF AMH concentrations were quantified using a sandwich enzyme-linked immunosorbent assay (ELISA). We also quantified 8 phthalate metabolite concentrations using tandem mass spectrometry in 1-2 urine samples per cycle (total 331 urines) and calculated the cycle-specific geometric mean for each metabolite. We applied cluster-weighted generalized estimating equation models (CWGEE) to evaluate the associations of tertiles of urinary phthalate metabolite concentrations with log-transformed FF AMH concentrations adjusting for potential confounders. Study participants had median age of 34.0 years (IQR 32.0, 37.0), 83% were white, and median BMI of 23.1 kg/m(2) (IQR 21.2, 26.1). The following stimulation protocols were used: luteal phase agonist (70%), antagonist (14%), or flare (16%). Urinary concentrations of select phthalate metabolites were negatively associated with FF AMH. For example, women whose urinary mEOHP was in the lowest tertile (range 0.30-4.04 ng/ml) had an adjusted mean FF AMH of 0.72 ng/mL (95% CI = 0.36, 1.44), compared to women in the highest tertile (range 9.90-235), who had an adjusted mean of 0.24 ng/mL (95% CI = 0.12-0.48, p < 0.05). The negative association between urinary concentrations of certain phthalate metabolites with FF AMH concentrations may have implications for antral follicle recruitment and fertility treatment outcomes.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/34375942/",
    "urlaid": "https://sci-hub.do/S0160-4120(21)00434-7 https://sci-hub.do/10.1016/j.envint.2021.106809",
    "pt": "Journal Article; Research Support, N.I.H., Extramural",
    "pl": "Netherlands",
    "annotations": null,
    "query_number": 1,
    "score": 0.91,
    "section": "title",
    "text": "Urinary phthalate metabolite concentrations are negatively associated with follicular fluid anti-mullerian hormone concentrations in women undergoing fertility treatment"
  },
  {
    "uid": 25873369,
    "aid": "1940-6207.CAPR-14-0377 10.1158/1940-6207.CAPR-14-0377",
    "titl": "Anti-Mullerian hormone concentrations in premenopausal women and breast cancer risk.",
    "mesh": "Adult;;; Anti-Mullerian Hormone/*blood;;; Biomarkers, Tumor/*blood;;; Breast;;; Breast Neoplasms/*blood/*diagnosis/epidemiology;;; Case-Control Studies;;; Enzyme-Linked Immunosorbent Assay;;; Female;;; Follow-Up Studies;;; Humans;;; Logistic Models;;; Middle Aged;;; Neoplasm Staging;;; Premenopause;;; Prognosis;;; Risk Factors;;; United States/epidemiology",
    "majr": "",
    "subh": "",
    "auth": "Nichols, Hazel B; Baird, Donna D; Stanczyk, Frank Z; Steiner, Anne Z; Troester, Melissa A; Whitworth, Kristina W; Sandler, Dale P",
    "jour": "Cancer prevention research (Philadelphia, Pa.)",
    "affl": "Department of Epidemiology, University of North Carolina Gillings School of Global Public Health, Chapel Hill, North Carolina. hazel.nichols@unc.edu.;;; Epidemiology Branch, National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina.;;; Departments of Obstetrics and Gynecology, and Preventive Medicine, University of Southern California Keck School of Medicine, Los Angeles, California.;;; Department of Obstetrics and Gynecology, University of North Carolina School of Medicine, Chapel Hill, North Carolina.;;; Department of Epidemiology, University of North Carolina Gillings School of Global Public Health, Chapel Hill, North Carolina.;;; The University of Texas School of Public Health, San Antonio Regional Campus, San Antonio, Texas.;;; Epidemiology Branch, National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina.",
    "pdat": "2015 Jun",
    "tiab": "Laboratory models support an inverse association between anti-Mullerian hormone (AMH) and breast tumor development. Human studies are lacking; one study (N = 105 cases, 204 controls) with prospectively collected serum reported the opposite-an approximate 10-fold increase in breast cancer risk comparing fourth with first quartile AMH levels. We investigated the relation between serum AMH levels and breast cancer risk in a case-control (N = 452 cases, 902 controls) study nested within the prospective Sister Study cohort of 50,884 women. At enrollment, participants were ages 35 to 54, premenopausal, and completed questionnaires on medical and family history, lifestyle factors, and demographics. AMH (ng/mL) was measured by ultrasensitive ELISA in serum collected at enrollment and log-transformed for analysis. Multivariate conditional logistic regression was used to calculate ORs and 95% confidence intervals (CI) to account for matching on age and enrollment year. Mean age at enrollment was 46.8 years with an average 2.9 years from blood draw to breast cancer diagnosis (SD = 1.9). AMH concentrations were below the limit of detection (0.003 ng/mL) for approximately 25% of samples. Compared with samples below the LOD, women with AMH >2.84 ng/mL (90th percentile among controls) had a 2-fold increase in breast cancer odds (OR, 2.25; 95% CI, 1.26-4.02). For each 1-unit increase in lnAMH, overall breast cancer odds increased by 8% (OR, 1.08; 95% CI, 1.02-1.15) and odds of estrogen receptor-positive, invasive disease increased by 15% (OR, 1.15; 95% CI, 1.05-1.25). Our findings demonstrate an overall positive relation between AMH and breast cancer.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/25873369/",
    "urlaid": "https://sci-hub.do/1940-6207.CAPR-14-0377 https://sci-hub.do/10.1158/1940-6207.CAPR-14-0377",
    "pt": "Comparative Study; Journal Article; Research Support, N.I.H., Extramural; Research Support, Non-U.S. Gov't",
    "pl": "United States",
    "annotations": null,
    "query_number": 1,
    "score": 0.9,
    "section": "title",
    "text": "Anti-Mullerian Hormone Concentrations in Premenopausal Women and Breast Cancer Risk"
  },
  {
    "uid": 30930424,
    "aid": "2019-003 10.1262/jrd.2019-003",
    "titl": "Efficacy of a single measurement of plasma anti-Mullerian hormone concentration for ovum pick-up donor selection of Japanese Black heifers in herd breeding programs.",
    "mesh": "Animals;;; Anti-Mullerian Hormone/analysis/*blood;;; Blood Chemical Analysis/methods/veterinary;;; *Breeding/methods;;; *Cattle;;; Cell Count;;; Donor Selection/*methods;;; Embryo Transfer/methods/veterinary;;; Female;;; Japan;;; Male;;; *Oocyte Retrieval/methods/veterinary;;; Oocytes/cytology;;; Predictive Value of Tests;;; Superovulation/blood",
    "majr": "",
    "subh": "",
    "auth": "Fushimi, Yasuo; Monniaux, Danielle; Takagi, Mitsuhiro",
    "jour": "The Journal of reproduction and development",
    "affl": "Guardian Co. Ltd., Kagoshima 890-0033, Japan.;;; Physiologie de la Reproduction, Centre INRA Val de Loire, 37380 Nouzilly, France.;;; Joint Faculty of Veterinary Medicine, Yamaguchi University, Yamaguchi 753-8515, Japan.",
    "pdat": "2019 Aug 9",
    "tiab": "In this study, we evaluated the efficiency of a single measurement of plasma anti-Mullerian hormone (AMH) concentration in heifers in determining the number of oocytes recoverable by ovum pick-up (OPU), and compared AMH concentrations among sister heifers from the same parents. For this, blood samples from 50 embryo-transfer-derived female Japanese Black (JB) heifers (mean: 8.7 age in months) were collected and plasma AMH concentration was measured. At 13-15 months of age, both the number of follicles (2-9 mm) and the number of collected oocytes after OPU were counted and compared. Results indicated that the heifers with the highest AMH concentration had the highest number of follicles in their ovaries and gave the highest number of collected oocytes with OPU, thereby indicating that a single measurement of plasma AMH concentration is informative for the selection of OPU-donor heifers in herd breeding programs. The practice of performing a single AMH measurement may accelerate the intensive breeding of JB herds.",
    "ptyp": "",
    "url": "https://pubmed.ncbi.nlm.nih.gov/30930424/",
    "urlaid": "https://sci-hub.do/2019-003 https://sci-hub.do/10.1262/jrd.2019-003",
    "pt": "Evaluation Study; Journal Article",
    "pl": "Japan",
    "annotations": null,
    "query_number": 1,
    "score": 0.9,
    "section": "METHODS",
    "text": "Measurement of plasma anti-Mullerian hormone concentration and classification of the heifers"
  }
]

const obj_color = {
  'disease': "#fdbbbb",
  'drug': "#ECC58B",
  'gene': "#E2DB8C",
  'chemical': "#21c354",
  'species': "#A6EFDC",
  'mutation': "#B2DDEA",
  'cell_type': "#C6DEF5",
  'cell_line': "#A3B3D2",
  'DNA': "#C9B9E8",
  'RNA': "#D7DBE8",
}

function markup_text(text, annotations) {
  if (!annotations) {
    return text
  }
  let markup_text = ''
  let last_position = 0
  for (let annotation of annotations) {
    let start = annotation.span.begin
    let end = annotation.span.end
    if (!annotation.prop) { console.log('This') }
    let markup_str = `<span style=\"color: ${obj_color[annotation.obj]}\">${text.slice(start, end)}<sub>${annotation.prob ? annotation.prob.toFixed(2) : ""}</sub></span>`
    markup_text = `${markup_text}${text.slice(last_position, start)}${markup_str}`

    last_position = end
  }
  return markup_text
}


export class DDIReview extends Component {

  constructor(props) {
    super(props);

    this.gridRef = createRef();
    this.gridAnaliseRef = createRef();
    this.state = {
      token: variables.token,
      loading: false,

      //queries
      query_list: [],
      message: null,
      articles: test_articles,
      articlesInfo: [
        { field: 'text', filter: 'agTextColumnFilter', editable: true, },
        { field: 'score', filter: 'agNumberColumnFilter', sortable: true, editable: true, },
        { field: 'query_number', editable: true, },
        { field: 'section', filter: 'agTextColumnFilter', editable: true, },
      ],
      summarise: null,
      task_id: null,

      // Filters
      queryText: 'What methods are available to measure anti-mullerian hormone concentrations in young women?',
      queryDate: 1,
      queryScore: 0.9,
      queryTypes: new Set(),
    }
  }

  getArticles = (task_id, query_number = 0, interval = 1000) => {
    fetch(variables.API_URL + `/api/ddi_review`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then(response => {
        console.log(query_number)
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          throw Error(response.statusText)
        }
      })
      .then(data => {
        if (data.data === null) {
          this.setState({ loading: true, message: data.message });
          setTimeout(() => {
            return this.getArticles(task_id, query_number, interval)
          }, interval);
        } else {
          this.setState({
            articles: [...this.state.articles, ...data.data], DetailArticle: data.data[0], loading: false, message: 'Выполненно'
          });
          if (query_number !== 0) {
            this.state.query_list[query_number - 1].status = 1;
          }
        }
      })
      .catch(error => {
        console.log(error);
        this.setState({ articles: [], DetailArticle: null, loading: false, message: 'Что-то пошло не так' });
        if (query_number !== 0) {
          this.state.query_list[query_number - 1].status = 2;
        }
      })
  }

  createTask() {
    // Отправляем запрос на сервер для получения статей
    this.state.query_list.push({ query: this.state.queryText, status: 0 });
    const query_number = this.state.query_list.length;
    fetch(variables.API_URL + '/api/ddi_review', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        query: this.state.queryText,
        score: this.state.queryScore,
        number_of_query: query_number,
        date: this.state.queryDate,
        type: [...this.state.queryTypes],
      })
    })
      .then(response => response.json())
      .then(data => {
        this.setState({
          task_id: data.data,
          message: 'Запрос начал обрабатываться'
        });
        alert("Ваш запрос в очереди. Пожайлуста дождитесь результата");
        this.getArticles(data.data, query_number)
      })
      .catch(error => {
        console.log(error);
        this.setState({ task: null, message: 'Что-то пошло не так' });
      }
      )
  }

  clearTask() {
    this.setState({ query_list: [], articles: [], DetailArticle: null })
    alert("Таблица очищена!");
  }

  componentDidMount() {
    this.getArticles();
    console.log('start');
  }

  onSelectionChanged = () => {
    const selectedRows = this.gridRef.current.api.getSelectedRows();
    this.setState({ DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  changeQueryText = (e) => {
    this.setState({ queryText: e.target.value });
  }

  changeQueryDate = (e) => {
    this.setState({ queryDate: e });
  }

  changeQueryTypes(type) {
    if (this.state.queryTypes.has(type)) {
      this.state.queryTypes.delete(type)
    } else {
      this.state.queryTypes.add(type)
    }
    this.setState({ updateOr: !this.state.updateOr })
  }

  // Summarise

  getSummarise = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/summarise_emb?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          this.setState({ loading: true })
          setTimeout(() => {
            return this.getSummarise(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        this.setState({
          summarise: data.data,
          message: 'Суммаризация прошла успешно'
        });
      })
      .catch((err) => {
        console.log(err);
        this.setState({ summarise: null, message: 'Произошла ошибка при суммаризации' });
      });
  }

  createSummariseQuery() {
    fetch(variables.API_URL + '/api/summarise_emb', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: 'some'
      })
    })
      .then((res) => {
        if (res.status == 200) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа' })
        this.getSummarise(task_id);
      })
      .catch((error) => {
        alert('Ошибка')
      })
  }

  // MarkUp article

  getMarkUp = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/markup?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          setTimeout(() => {
            return this.getMarkUp(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        try {
          this.setState({
            DetailArticle: data.data,
            message: 'Разметка прошла успешно',
            loading: false,
          });
        } catch {
          console.log('access')
        }
      })
      .catch((err) => {
        console.log(err);
        this.setState({ message: 'Произошла ошибка при разметке', loading: false, });
      });
  }

  markUpArticle(DetailArticle) {
    this.setState({ loading: true })
    fetch(variables.API_URL + '/api/markup', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        article: DetailArticle
      })
    })
      .then((res) => {
        console.log(res.status)
        if (res.ok) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа' })
        this.getMarkUp(task_id);
      })
      .catch((err) => {
        console.log(err);
        this.setState({
          message: 'ошибка при разметке',
          loading: false,
        });
      });
  }

  suppressCutToClipboard = false;

  onRemoveSelected = () => {

    const selectedData = this.gridRef.current.api.getSelectedRows();
    console.log(selectedData)
    const res = this.gridRef.current.api.applyTransaction({ remove: selectedData });
  }

  onCellValueChanged = (params) => {
    console.log('Callback onCellValueChanged:', params);
    console.log(params.node)
    const res = this.gridRef.current.api.applyTransaction({ remove: [params.node.data] });
  }

  onCutStart = (params) => {
    console.log('Callback onCutStart:', params);
  }

  onCutEnd = (params) => {
    console.log('Callback onCutEnd:', params);
  }

  getRowId = () => {
    return (params) => {
      console.log(params)
      return params.data.code;
    };
  }


  render() {
    const {
      token,
      loading,
      query_list,
      articlesInfo,
      articles,
      DetailArticle,
      message,
      summarise,

      queryText,
      queryDate,
      queryScore,
    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />
    } else {
      return (
        <>
          <header>
            <nav class="bg-white border-gray-200 px-4 lg:px-6 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img src="https://flowbite.s3.amazonaws.com/logo.svg" class="mr-3 h-8" alt="FlowBite Logo" />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">EBM DаtaMed</span>
                  </a>
                  <ul class="flex font-medium flex-row space-x-8">
                    <Link to="/tematic_review">
                      <li>
                        <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 bg-blue-700 rounded md:bg-transparent md:text-blue-700 md:p-0" aria-current="page">Тематический анализ</a>
                      </li>
                    </Link>
                    <Link to="/ddi_review">
                      <li>
                        <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Факты для EBM</a>
                      </li>
                    </Link>
                  </ul>
                </div>
                <div class="flex items-center lg:order-2">
                  <button type="button" class="hidden sm:inline-flex items-center justify-center text-white bg-primary-700 hover:bg-primary-800 focus:ring-4 focus:ring-primary-300 font-medium rounded-lg text-xs px-3 py-1.5 mr-2 dark:bg-primary-600 dark:hover:bg-primary-700 focus:outline-none dark:focus:ring-primary-800"><svg aria-hidden="true" class="mr-1 -ml-1 w-5 h-5" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M10 5a1 1 0 011 1v3h3a1 1 0 110 2h-3v3a1 1 0 11-2 0v-3H6a1 1 0 110-2h3V6a1 1 0 011-1z" clip-rule="evenodd"></path></svg> Действие</button>
                </div>
              </div>
            </nav>
            <nav class="bg-white border-gray-200 px-6">
              <div class="max-w-screen-xl py-3">
                <div class="flex items-start">
                  <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar" class="hidden p-2 mr-3 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-label="Toggle navigation">
                    <svg class="w-6 h-6" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                  </button>
                  <label for="topbar-search" class="sr-only">Поисковый запрос</label>
                  <div className='w-full'>
                    <label for="search" class="mb-2 text-sm font-medium text-gray-900 sr-only dark:text-white">Поисковый запрос</label>
                    <div class="relative mt-1 w-full">
                      <div class="absolute inset-y-0 left-0 flex items-center pl-3 pointer-events-none">
                        <svg class="w-4 h-4 text-gray-500 dark:text-gray-400" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 20 20">
                          <path stroke="currentColor" stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="m19 19-4-4m0-7A7 7 0 1 1 1 8a7 7 0 0 1 14 0Z" />
                        </svg>
                      </div>
                      <input
                        class="py-3 bg-gray-50 border border-gray-300 text-gray-900 sm:text-sm rounded-lg focus:ring-primary-500 focus:border-primary-500 block w-full pl-10 p-2.5"
                        id="search"
                        type="text"
                        name="search_field"
                        placeholder="Поисковый запрос"
                        value={queryText}
                        onChange={this.changeQueryText}
                        aria-label="Search" />
                      <button type="submit" value="Найти" onClick={() => this.createTask()} class="text-white absolute right-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2">Найти</button>
                    </div>
                  </div>
                </div>
              </div>
            </nav>
          </header >
          <main>
            <div>
              <div className="container-fluid">
                <div className="row align-items-stretch b-height">
                  <aside id="sidebar" className="h-screen col-md-2 my-3 bg-white collapse show width border rounded-3 g-0">
                    <div className="accordion accordion-flush" id="accordionFlushExample">
                      <div className="accordion-item">
                        <h2 className="accordion-header" id="flush-headingOne">
                          <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseOne" aria-expanded="false" aria-controls="flush-collapseOne">
                            Период поиска
                          </button>
                        </h2>
                        <div id="flush-collapseOne" className="collapse show multi-collapse" aria-labelledby="flush-headingOne" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault1"
                                checked={queryDate === 1}
                                onChange={() => this.changeQueryDate(1)}
                              />
                              <label class="form-check-label" for="flexRadioDefault1">
                                1 год
                              </label>
                            </div>
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault2"
                                checked={queryDate === 3}
                                onChange={() => this.changeQueryDate(3)}
                              />
                              <label class="form-check-label" for="flexRadioDefault2">
                                3 года
                              </label>
                            </div>
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault3"
                                checked={queryDate === 5}
                                onChange={() => this.changeQueryDate(5)}
                              /><label class="form-check-label" for="flexRadioDefault3">
                                5 лет
                              </label>
                            </div>
                            <div class="form-check">
                              <input
                                className="form-check-input"
                                type="radio"
                                name="flexRadioDefault"
                                id="flexRadioDefault1"
                                checked={!queryDate}
                                onChange={() => this.changeQueryDate(null)}
                              />
                              <label class="form-check-label" for="flexRadioDefault4">
                                            > 5 лет
                              </label>
                            </div>
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 class="accordion-header" id="flush-headingThree">
                          <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFour" aria-expanded="false" aria-controls="flush-collapseThree">
                            Тип статьи
                          </button>
                        </h2>
                        <div id="flush-collapseFour" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxClinicalTrial"
                                name="CheckBoxClinicalTrial"
                                checked={this.state.queryTypes.has('Clinical Trial')}
                                onChange={() => this.changeQueryTypes('Clinical Trial')}
                              />
                              <label className="form-check-label" for="CheckboxClinicalTrial">Clinical Trial</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxMetaAnalysys"
                                name="CheckboxMetaAnalysys"
                                checked={this.state.queryTypes.has('Meta Analysys')}
                                onChange={() => this.changeQueryTypes('Meta Analysys')}
                              />
                              <label className="form-check-label" for="CheckboxMetaAnalysys">Meta Analysys</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxRandomizedControlledTrial"
                                name="CheckboxRandomizedControlledTrial"
                                checked={this.state.queryTypes.has('Randomized Controlled Trial')}
                                onChange={() => this.changeQueryTypes('Randomized Controlled Trial')}
                              />
                              <label className="form-check-label" for="CheckboxRandomizedControlledTrial">Randomized Controlled Trial</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxReview"
                                name="CheckboxReview"
                                checked={this.state.queryTypes.has('Review')}
                                onChange={() => this.changeQueryTypes('Review')}
                              />
                              <label className="form-check-label" for="CheckboxReview">Review</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxSystematicReview"
                                name="CheckboxSystematicReview"
                                checked={this.state.queryTypes.has('Systematic Review')}
                                onChange={() => this.changeQueryTypes('Systematic Review')}
                              />
                              <label className="form-check-label" for="CheckboxSystematicReview">Systematic Review</label>
                            </div>
                            <div className="form-check form-check-inline">
                              <input
                                className="form-check-input"
                                type="checkbox"
                                id="CheckboxJournalArticle"
                                name="CheckboxJournalArticle"
                                checked={this.state.queryTypes.has('Journal Article')}
                                onChange={() => this.changeQueryTypes('Journal Article')}
                              />
                              <label className="form-check-label" for="CheckboxJournalArticle">Journal Article</label>
                            </div>

                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 class="accordion-header" id="flush-headingThree">
                          <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseThree" aria-expanded="false" aria-controls="flush-collapseThree">
                            Точность
                          </button>
                        </h2>
                        <div id="flush-collapseThree" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            <p>Требуемая точность = {queryScore.toFixed(2)}</p>
                            <Slider
                              axis="x"
                              x={queryScore}
                              xmax={1}
                              xmin={0}
                              xstep={0.01}
                              onChange={({ x }) => this.setState({ queryScore: x })}
                            />
                          </div>
                        </div>
                      </div>
                      <div className="accordion-item">
                        <h2 class="accordion-header" id="flush-headingThree">
                          <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFive" aria-expanded="false" aria-controls="flush-collapseThree">
                            Выделенные сущности
                          </button>
                        </h2>
                        <div id="flush-collapseFive" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                          <div className="accordion-body">
                            {Object.entries(obj_color).map(tag => <p class="pb-2 mb-3 border-bottom" style={{ color: `${tag[1]}` }}>{tag[0]}.</p>)}
                          </div>
                        </div>
                      </div>
                    </div>
                  </aside>
                  <section class="col p-3 m-3 border rounded-3 bg-white overflow-auto">
                    <div class="accordion accordion-flush" id="accordion">
                      <div class="accordion-item">
                        <h2 class="accordion-header" id="">
                          <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSeven" aria-expanded="false" aria-controls="flush-collapseSeven">
                            Запросы. {message}
                          </button>
                        </h2>
                        <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">
                          <div class="accordion-body">
                            {query_list?.map((query, index) =>
                              query.status === 2 ?
                                <p class="pb-2 mb-3 border-bottom" style={{ color: 'red' }}>{index + 1} - {query.query}.</p>
                                : query.status === 1 ?
                                  <p class="pb-2 mb-3 border-bottom" style={{ color: 'green' }}>{index + 1} - {query.query}.</p>
                                  :
                                  <p class="pb-2 mb-3 border-bottom" style={{ color: 'black' }}>{index + 1} - {query.query}.</p>
                            )}
                          </div>
                          <div class="accordion-body">
                            <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Очистить" onClick={() => this.clearTask()} />
                          </div>
                        </div>
                      </div>
                      <div>
                        {summarise ?
                          <>
                            <p>Summarise</p>
                            <p>{summarise}</p>
                          </>
                          : <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Суммаризовать" onClick={() => this.createSummariseQuery()} />}
                      </div>
                    </div>
                    <div>
                      <div class="bd-example">
                        <div class="tab-content" id="myTabContent">
                          <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                            <div class="container-fluid g-0">
                              <div className="ag-theme-alpine ag-theme-acmecorp" style={{ height: 700 }}>
                                <AgGridReact
                                  ref={this.gridRef}
                                  rowData={articles}
                                  columnDefs={articlesInfo}
                                  pagination={true}
                                  rowSelection={'single'}
                                  onSelectionChanged={this.onSelectionChanged}
                                  suppressCutToClipboard={this.suppressCutToClipboard}
                                  onCellValueChanged={this.onCellValueChanged}
                                  onCutStart={this.onCutStart}

                                  onCutEnd={this.onCutEnd}
                                  sideBar={{
                                    toolPanels: [
                                      {
                                        id: 'columns',
                                        labelDefault: 'Columns',
                                        labelKey: 'columns',
                                        iconKey: 'columns',
                                        toolPanel: 'agColumnsToolPanel',
                                        minWidth: 225,
                                        width: 225,
                                        maxWidth: 225,
                                      },
                                      {
                                        id: 'filters',
                                        labelDefault: 'Filters',
                                        labelKey: 'filters',
                                        iconKey: 'filter',
                                        toolPanel: 'agFiltersToolPanel',
                                        minWidth: 180,
                                        maxWidth: 400,
                                        width: 250,
                                      },
                                    ],
                                    position: 'left',

                                  }}
                                >
                                </AgGridReact>
                              </div>
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </section>

                  <aside id="sidebar2" class="col-md-4 h-screen collapse show width col p-3 my-3 border rounded-3 bg-white">
                    <h3 class="pb-2 mb-3 border-bottom">Подробное описание</h3>
                    <nav class="small" id="toc">
                      {DetailArticle ?
                        <div class="card mb-3">
                          <div class="card-body">
                            <a href={DetailArticle.url} class="card-title link-primary text-decoration-none h5" target="_blank"> {DetailArticle.titl} </a>
                            <p class="card-text">---------------------------------- </p>
                            <p class="card-text">Авторы :  {DetailArticle.auth} </p>
                            <p class="card-text">---------------------------------- </p>
                            <p class="card-text">Аннотация :  </p>
                            <p class="card-text" dangerouslySetInnerHTML={{ __html: markup_text(DetailArticle.tiab, DetailArticle.annotations) }} />
                            <p class="card-text">---------------------------------- </p>
                            <p class="card-text"><small class="text-success">Дата публикации : {DetailArticle.pdat} </small></p>
                            <p class="card-text"><small class="text-success">Издание : {DetailArticle.jour}</small></p>
                            <p class="card-text"><small class="text-success">Вид публикации : {DetailArticle.pt}</small></p>
                            <p class="card-text"><small class="text-success">Страна : {DetailArticle.pl} </small></p>
                            <p class="card-text"><small class="text-success">{DetailArticle.mesh} </small></p>
                            {summarise ?
                              <>
                                <p>Summarise</p>
                                <p>{summarise}</p>
                              </>
                              : loading ?
                                <p>Loading...</p>
                                : <input className="text-white right-2.5 my-4 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Разметить" onClick={() => this.markUpArticle(DetailArticle)} />
                            }
                            <input className="text-white right-2.5 my-4 bottom-2.5 bg-red-700 hover:bg-red-800 focus:ring-4 focus:outline-none focus:ring-red-300 font-medium rounded-lg text-sm px-4 py-2" type="submit" value="Удалить" onClick={() => this.onRemoveSelected()} />
                          </div>
                        </div>
                        : null}
                    </nav>
                  </aside>

                </div>
              </div>
            </div>
          </main>
        </>
      )
    }
  }
}