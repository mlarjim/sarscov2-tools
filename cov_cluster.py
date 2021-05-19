#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Para obtener los datos para recrear la gráfica de proporcion de numero total de secuencias, con el tiempo, que entran dentro de mutaciones definidas


# In[2]:


#para la d614


# In[1]:


import pandas as pd
import glob as g
from collections import Counter
import argparse
import re


# In[56]:


#tiene que leer:
# metadata
# tsv de ivar
# muestras de españa del gisaid


# In[ ]:


parser = argparse.ArgumentParser(description = "Devuelve una tabla con la frecuencia de las mutaciones ecpecificadas por semana. Se necesita una tabla de metadatos de GISAID por cada mutacion que deseemos incorporar al estudio (esa tabla solo contiene las muestras con esa mutacion), y cada tabla debera llamarse con el numero de la posición + '_spain.tsv' Issue: solo alcanza hasta fechas de la tercera semana de enero 2021. Issue2: solo para aminoacidos de la spike")
parser.add_argument("metadata", help= "tsv with (at least) columns 'Strain' 'date' and 'Province'")
parser.add_argument("directory1", help= "'/complete/directory/to/*_analysis_report.csv'")
#parser.add_argument("directory2", help = "'/complete/directory/to/analysis_MN908947/'")
parser.add_argument("position", help = "comma separated aminoacid positions of protein s")
parser.add_argument("--metadata_spain", help = "tsv with (at least) column 'Collection date' from GISAID spanish samples")

args = parser.parse_args()


# In[ ]:





# In[ ]:





# In[ ]:





# In[57]:


####             ANDALUCIA


# In[ ]:





# In[58]:


#tabla = pd.read_csv("metadata.tsv", sep = "\t")
tabla = pd.read_csv(args.metadata, sep = "\t")


# In[59]:


#qc = g.glob("/home/marialara/gattaca/all_run_QC_results/*_analysis_report.csv")
qc = g.glob(args.directory1)


# In[60]:


#text = g.glob("/home/marialara/gattaca/gtc/*.tsv")
#mn = args.directory2
mn = "../covid_gattaca_ola1/"

# In[61]:


yes = []
for w in range(len(qc)):
    csv = pd.read_csv(qc[w], sep=",") #lee una de las tablas (csv)
    
    for i in range(len(csv)):
        if csv["selected_for_nextstrain"][i] == "yes": #si está seleccionada para nextstrain
            run = qc[w].split("/")[-1].split("_")[0]	
            num = qc[w].split("/")[-1].split("_")[1]	
            tog = [run, num]
            rute = "_".join(tog)
            ap = mn + rute + "_batch*/variants/ivar/" + csv["sample"][i] +".tsv" #ruta de muestras que han pasado el qc
            yes.append(ap) #guarda las rutas

gisaid = g.glob("/home/mlara/covid_gattaca_ola1/run_gisaid*/variants/ivar/*.tsv") #guarda directorio de los falsos .tsv de gisaid

yes = yes + gisaid
total = len(yes)


# In[ ]:


countm = 0
for i in range(len(tabla)):
    for y in yes:
        id = g.glob(y)[0].split("/")[-1].split(".")[0].split("_")[-1]
        if bool(re.search(id, tabla["Strain"][i]))== True:
            countm = countm +1


# In[ ]:





# In[3]:


### mutaciones de interes en la spike: N501, V1122L
posi = args.position

muts = []
for i in posi.split(","):
    muts.append(int(i))


# In[73]:


week = ["ult_Febre",
                    "seg_Marzo",
                    "terc_marzo",
                    "cuar_marzo",
                    "quint_marzo",
                    "sext_marzo",
                    "seg_abril",
                    "terc_abril",
                    "cuar_abril",
                    "quint_abril",
                    "seg_mayo",
                    "terc_mayo",
                    "cuar_mayo",
                    "quint_mayo",
                    "pri_junio",
                  "seg_junio",
                  "terc_junio",
                  "cuar_junio",
                  "pri_julio",
                  "seg_julio",
                  "terc_julio",
                  "cuar_julio",
                  "quint_julio",
                  "pri_agosto",
                  "seg_agosto",
                  "terc_agosto",
                  "cuar_agosto",
                  "pri_sept",
                  "seg_sept",
                  "terc_sept",
                  "cuar_sept",
                  "pri_octu",
                  "seg_octu",
                  "terc_octu",
                  "cuar_octu",
                  "quint_octu",
                  "pri_nov",
                  "seg_nov",
                  "terc_nov",
                  "cuar_nov",
                  "pri_dici",
                  "seg_dici",
                  "terc_dici",
                  "cuar_dici",
                  "quint_dici",
                  "pri_ene_2021",
                  "seg_ene_2021",
                  "terc_ene_2021"
       ]


# In[ ]:





# In[ ]:





# In[75]:


data = {"week" : week
       }

df = pd.DataFrame (data, columns = ["week"])
for m in muts: 
    samples = [] #muestras con la mutacion 
    for a in yes:
        id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta del directorio
        #for b in yes: #recorre los id de yes
            #if b == id: #si esa muestra ha pasado el control
        #ruta = str(mn) + str(a)
        if bool(re.match("HUVR_(\d*)UK", id)) == False: #si son de la primera ola
            tsv = pd.read_csv(g.glob(a)[0], sep="\t") #lee el tsv de esa muestra
            for i in range(len(tsv)):
                if  tsv["ALT_FREQ"][i]>=0.75 and int((tsv["POS"][i] - 21562 +2)/3) == m and tsv["REF_AA"][i] != tsv["ALT_AA"][i]: 
                    #si tienen mutación en el aa deseado en la spike (SOLO!!!)
                    samples.append(g.glob(a)[0].split("/")[-1].split(".")[0]) #extrae el id de la muestra de la ruta del directorio
        

    ult_Febre= 0
    seg_Marzo = 0
    cuar_abril = 0
    cuar_marzo = 0
    quint_abril = 0
    quint_marzo = 0
    seg_abril = 0
    seg_mayo = 0
    sext_marzo = 0
    terc_abril = 0
    terc_marzo = 0

    terc_mayo = 0
    cuar_mayo = 0
    quint_mayo = 0
    pri_junio = 0
    seg_junio= 0
    terc_junio= 0
    cuar_junio= 0
    pri_julio= 0
    seg_julio= 0
    terc_julio= 0
    cuar_julio= 0
    quint_julio= 0
    pri_agosto= 0
    seg_agosto= 0
    terc_agosto= 0
    cuar_agosto= 0
    pri_sept= 0
    seg_sept= 0
    terc_sept= 0
    cuar_sept= 0
    pri_octu= 0
    seg_octu= 0
    terc_octu= 0
    cuar_octu= 0
    quint_octu= 0
    pri_nov= 0
    seg_nov= 0
    terc_nov= 0
    cuar_nov= 0
    pri_dici= 0
    seg_dici= 0
    terc_dici= 0
    cuar_dici= 0
    quint_dici= 0
    pri_ene_2021= 0
    seg_ene_2021= 0
    terc_ene_2021= 0
    
    

    for i in range(len(tabla)):
        for a in samples:
                    id = a.split("_")[-1]
                    
                    if bool(re.search(id, tabla["Strain"][i]))== True:
                        if 20200224 <= int("".join(tabla["date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
                            ult_Febre+=1
                        elif 20200302 <= int("".join(tabla["date"][i].split("-"))) <= 20200308:
                            seg_Marzo+=1
                        elif 20200309 <= int("".join(tabla["date"][i].split("-"))) <= 20200315:
                            terc_marzo+=1        
                        elif 20200316 <= int("".join(tabla["date"][i].split("-"))) <= 20200322:
                            cuar_marzo+=1    
                        elif 20200323 <= int("".join(tabla["date"][i].split("-"))) <= 20200329:
                            quint_marzo+=1
                        elif 20200330 <= int("".join(tabla["date"][i].split("-"))) <= 20200405:
                            sext_marzo+=1
                        elif 20200406 <= int("".join(tabla["date"][i].split("-"))) <= 20200412:
                            seg_abril+=1 
                        elif 20200413 <= int("".join(tabla["date"][i].split("-"))) <= 20200419:
                            terc_abril+=1     
                        elif 20200420 <= int("".join(tabla["date"][i].split("-"))) <= 20200426:
                            cuar_abril+=1   
                        elif 20200427 <= int("".join(tabla["date"][i].split("-"))) <= 20200503:
                            quint_abril+=1
                        elif 20200504 <= int("".join(tabla["date"][i].split("-"))) <= 20200510:
                            seg_mayo+=1
                        elif 20200511 <= int("".join(tabla["date"][i].split("-"))) <= 20200517:
                            terc_mayo+=1
                        elif 20200518 <= int("".join(tabla["date"][i].split("-"))) <= 20200524:
                            cuar_mayo+=1
                        elif 20200525 <= int("".join(tabla["date"][i].split("-"))) <= 20200531:
                            quint_mayo+=1
                        elif 20200601 <= int("".join(tabla["date"][i].split("-"))) <= 20200607:
                            pri_junio+=1
                        elif 20200608 <= int("".join(tabla["date"][i].split("-"))) <= 20200614:
                            seg_junio+=1
                        elif 20200615 <= int("".join(tabla["date"][i].split("-"))) <= 20200621:
                            terc_junio+=1
                        elif 20200622 <= int("".join(tabla["date"][i].split("-"))) <= 20200628:
                            cuar_junio+=1
                        elif 20200629 <= int("".join(tabla["date"][i].split("-"))) <= 20200705:
                            pri_julio+=1
                        elif 20200706 <= int("".join(tabla["date"][i].split("-"))) <= 20200712:
                            seg_julio+=1
                        elif 20200713 <= int("".join(tabla["date"][i].split("-"))) <= 20200719:
                            terc_julio+=1
                        elif 20200720 <= int("".join(tabla["date"][i].split("-"))) <= 20200726:
                            cuar_julio+=1
                        elif 20200727 <= int("".join(tabla["date"][i].split("-"))) <= 20200802:
                            quint_julio+=1
                        elif 20200803 <= int("".join(tabla["date"][i].split("-"))) <= 20200809:
                            pri_agosto+=1
                        elif 20200810 <= int("".join(tabla["date"][i].split("-"))) <= 20200816:
                            seg_agosto+=1
                        elif 20200817 <= int("".join(tabla["date"][i].split("-"))) <= 20200823:
                            terc_agosto+=1
                        elif 20200824 <= int("".join(tabla["date"][i].split("-"))) <= 20200830:
                            cuar_agosto+=1
                        elif 20200831 <= int("".join(tabla["date"][i].split("-"))) <= 20200906:
                            pri_sept+=1
                        elif 20200907 <= int("".join(tabla["date"][i].split("-"))) <= 20200913:
                            seg_sept+=1
                        elif 20200914 <= int("".join(tabla["date"][i].split("-"))) <= 20200920:
                            terc_sept+=1
                        elif 20200921 <= int("".join(tabla["date"][i].split("-"))) <= 20200927:
                            cuar_sept+=1
                        elif 20200928 <= int("".join(tabla["date"][i].split("-"))) <= 20201004:
                            pri_octu+=1
                        elif 20201005 <= int("".join(tabla["date"][i].split("-"))) <= 20201011:
                            seg_octu+=1
                        elif 20201012 <= int("".join(tabla["date"][i].split("-"))) <= 20201018:
                            terc_octu+=1
                        elif 20201019 <= int("".join(tabla["date"][i].split("-"))) <= 20201025:
                            cuar_octu+=1
                        elif 20201026 <= int("".join(tabla["date"][i].split("-"))) <= 20201101:
                            quint_octu+=1
                        elif 20201102 <= int("".join(tabla["date"][i].split("-"))) <= 20201108:
                            pri_nov+=1
                        elif 20201109 <= int("".join(tabla["date"][i].split("-"))) <= 20201115:
                            seg_nov+=1
                        elif 20201116 <= int("".join(tabla["date"][i].split("-"))) <= 20201122:
                            terc_nov+=1
                        elif 20201123 <= int("".join(tabla["date"][i].split("-"))) <= 20201129:
                            cuar_nov+=1
                        elif 20201130 <= int("".join(tabla["date"][i].split("-"))) <= 20201206:
                            pri_dici+=1
                        elif 20201207 <= int("".join(tabla["date"][i].split("-"))) <= 20201213:
                            seg_dici+=1
                        elif 20201214 <= int("".join(tabla["date"][i].split("-"))) <= 20201220:
                            terc_dici+=1
                        elif 20201221 <= int("".join(tabla["date"][i].split("-"))) <= 20201227:
                            cuar_dici+=1
                        elif 20201228 <= int("".join(tabla["date"][i].split("-"))) <= 20210103:
                            quint_dici+=1
                        elif 20210104 <= int("".join(tabla["date"][i].split("-"))) <= 20210110:
                            pri_ene_2021+=1
                        elif 20210111 <= int("".join(tabla["date"][i].split("-"))) <= 20210117:
                            seg_ene_2021+=1
                        elif 20210118 <= int("".join(tabla["date"][i].split("-"))) <= 20210124:
                            terc_ene_2021+=1
                            


    tult_Febre= 0
    tseg_Marzo = 0
    tcuar_abril = 0
    tcuar_marzo = 0
    tquint_abril = 0
    tquint_marzo = 0
    tseg_abril = 0
    tseg_mayo = 0
    tsext_marzo = 0
    tterc_abril = 0
    tterc_marzo = 0

    tterc_mayo = 0
    tcuar_mayo = 0
    tquint_mayo = 0
    tpri_junio = 0
    tseg_junio = 0
    tterc_junio = 0
    tcuar_junio = 0
    tpri_julio = 0
    tseg_julio = 0
    tterc_julio = 0
    tcuar_julio = 0
    tquint_julio = 0
    tpri_agosto = 0
    tseg_agosto = 0
    tterc_agosto = 0
    tcuar_agosto = 0
    tpri_sept = 0
    tseg_sept = 0
    tterc_sept = 0
    tcuar_sept = 0
    tpri_octu = 0
    tseg_octu = 0
    tterc_octu = 0
    tcuar_octu = 0
    tquint_octu = 0
    tpri_nov = 0
    tseg_nov = 0
    tterc_nov = 0
    tcuar_nov = 0
    tpri_dici = 0
    tseg_dici = 0
    tterc_dici = 0
    tcuar_dici = 0
    tquint_dici = 0
    tpri_ene_2021 = 0
    tseg_ene_2021 = 0
    tterc_ene_2021 = 0


    for i in range(len(tabla)):
        #if tabla["Source"][i] == "SAS":
            if 20200224 <= int("".join(tabla["date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
                tult_Febre+=1
            elif 20200302 <= int("".join(tabla["date"][i].split("-"))) <= 20200308:
                tseg_Marzo+=1
            elif 20200309 <= int("".join(tabla["date"][i].split("-"))) <= 20200315:
                tterc_marzo+=1        
            elif 20200316 <= int("".join(tabla["date"][i].split("-"))) <= 20200322:
                tcuar_marzo+=1    
            elif 20200323 <= int("".join(tabla["date"][i].split("-"))) <= 20200329:
                tquint_marzo+=1
            elif 20200330 <= int("".join(tabla["date"][i].split("-"))) <= 20200405:
                tsext_marzo+=1
            elif 20200406 <= int("".join(tabla["date"][i].split("-"))) <= 20200412:
                tseg_abril+=1 
            elif 20200413 <= int("".join(tabla["date"][i].split("-"))) <= 20200419:
                tterc_abril+=1     
            elif 20200420 <= int("".join(tabla["date"][i].split("-"))) <= 20200426:
                tcuar_abril+=1   
            elif 20200427 <= int("".join(tabla["date"][i].split("-"))) <= 20200503:
                tquint_abril+=1
            elif 20200504 <= int("".join(tabla["date"][i].split("-"))) <= 20200510:
                tseg_mayo+=1
            elif 20200511 <= int("".join(tabla["date"][i].split("-"))) <= 20200517:
                tterc_mayo+=1
            elif 20200518 <= int("".join(tabla["date"][i].split("-"))) <= 20200524:
                tcuar_mayo+=1
            elif 20200525 <= int("".join(tabla["date"][i].split("-"))) <= 20200531:
                tquint_mayo+=1
            elif 20200601 <= int("".join(tabla["date"][i].split("-"))) <= 20200607:
                tpri_junio+=1
            elif 20200608 <= int("".join(tabla["date"][i].split("-"))) <= 20200614:
                tseg_junio+=1
            elif 20200615 <= int("".join(tabla["date"][i].split("-"))) <= 20200621:
                tterc_junio+=1
            elif 20200622 <= int("".join(tabla["date"][i].split("-"))) <= 20200628:
                tcuar_junio+=1
            elif 20200629 <= int("".join(tabla["date"][i].split("-"))) <= 20200705:
                tpri_julio+=1
            elif 20200706 <= int("".join(tabla["date"][i].split("-"))) <= 20200712:
                tseg_julio+=1
            elif 20200713 <= int("".join(tabla["date"][i].split("-"))) <= 20200719:
                tterc_julio+=1
            elif 20200720 <= int("".join(tabla["date"][i].split("-"))) <= 20200726:
                tcuar_julio+=1
            elif 20200727 <= int("".join(tabla["date"][i].split("-"))) <= 20200802:
                tquint_julio+=1
            elif 20200803 <= int("".join(tabla["date"][i].split("-"))) <= 20200809:
                tpri_agosto+=1
            elif 20200810 <= int("".join(tabla["date"][i].split("-"))) <= 20200816:
                tseg_agosto+=1
            elif 20200817 <= int("".join(tabla["date"][i].split("-"))) <= 20200823:
                tterc_agosto+=1
            elif 20200824 <= int("".join(tabla["date"][i].split("-"))) <= 20200830:
                tcuar_agosto+=1
            elif 20200831 <= int("".join(tabla["date"][i].split("-"))) <= 20200906:
                tpri_sept+=1
            elif 20200907 <= int("".join(tabla["date"][i].split("-"))) <= 20200913:
                tseg_sept+=1
            elif 20200914 <= int("".join(tabla["date"][i].split("-"))) <= 20200920:
                tterc_sept+=1
            elif 20200921 <= int("".join(tabla["date"][i].split("-"))) <= 20200927:
                tcuar_sept+=1
            elif 20200928 <= int("".join(tabla["date"][i].split("-"))) <= 20201004:
                tpri_octu+=1
            elif 20201005 <= int("".join(tabla["date"][i].split("-"))) <= 20201011:
                tseg_octu+=1
            elif 20201012 <= int("".join(tabla["date"][i].split("-"))) <= 20201018:
                tterc_octu+=1
            elif 20201019 <= int("".join(tabla["date"][i].split("-"))) <= 20201025:
                tcuar_octu+=1
            elif 20201026 <= int("".join(tabla["date"][i].split("-"))) <= 20201101:
                tquint_octu+=1
            elif 20201102 <= int("".join(tabla["date"][i].split("-"))) <= 20201108:
                tpri_nov+=1
            elif 20201109 <= int("".join(tabla["date"][i].split("-"))) <= 20201115:
                tseg_nov+=1
            elif 20201116 <= int("".join(tabla["date"][i].split("-"))) <= 20201122:
                tterc_nov+=1
            elif 20201123 <= int("".join(tabla["date"][i].split("-"))) <= 20201129:
                tcuar_nov+=1
            elif 20201130 <= int("".join(tabla["date"][i].split("-"))) <= 20201206:
                tpri_dici+=1
            elif 20201207 <= int("".join(tabla["date"][i].split("-"))) <= 20201213:
                tseg_dici+=1
            elif 20201214 <= int("".join(tabla["date"][i].split("-"))) <= 20201220:
                tterc_dici+=1
            elif 20201221 <= int("".join(tabla["date"][i].split("-"))) <= 20201227:
                tcuar_dici+=1
            elif 20201228 <= int("".join(tabla["date"][i].split("-"))) <= 20210103:
                tquint_dici+=1
            elif 20210104 <= int("".join(tabla["date"][i].split("-"))) <= 20210110:
                tpri_ene_2021+=1
            elif 20210111 <= int("".join(tabla["date"][i].split("-"))) <= 20210117:
                tseg_ene_2021+=1
            elif 20210118 <= int("".join(tabla["date"][i].split("-"))) <= 20210124:
                tterc_ene_2021+=1

                    
    total_and = [
    tult_Febre,
    tseg_Marzo,
    tcuar_abril,
    tcuar_marzo,
    tquint_abril,
    tquint_marzo,
    tseg_abril,
    tseg_mayo,
    tsext_marzo,
    tterc_abril,
    tterc_marzo,

    tterc_mayo,
    tcuar_mayo,
    tquint_mayo,
    tpri_junio,
    tseg_junio,
    tterc_junio,
    tcuar_junio,
    tpri_julio,
    tseg_julio,
    tterc_julio,
    tcuar_julio,
    tquint_julio,
    tpri_agosto,
    tseg_agosto,
    tterc_agosto,
    tcuar_agosto,
    tpri_sept,
    tseg_sept,
    tterc_sept,
    tcuar_sept,
    tpri_octu,
    tseg_octu,
    tterc_octu,
    tcuar_octu,
    tquint_octu,
    tpri_nov,
    tseg_nov,
    tterc_nov,
    tcuar_nov,
    tpri_dici,
    tseg_dici,
    tterc_dici,
    tcuar_dici,
    tquint_dici,
    tpri_ene_2021,
    tseg_ene_2021,
    tterc_ene_2021
    ]

    parc_and = [
    ult_Febre,
    seg_Marzo,
    cuar_abril,
    cuar_marzo,
    quint_abril,
    quint_marzo,
    seg_abril,
    seg_mayo,
    sext_marzo,
    terc_abril,
    terc_marzo,

    terc_mayo,
    cuar_mayo,
    quint_mayo,
    pri_junio,
    seg_junio,
    terc_junio,
    cuar_junio,
    pri_julio,
    seg_julio,
    terc_julio,
    cuar_julio,
    quint_julio,
    pri_agosto,
    seg_agosto,
    terc_agosto,
    cuar_agosto,
    pri_sept,
    seg_sept,
    terc_sept,
    cuar_sept,
    pri_octu,
    seg_octu,
    terc_octu,
    cuar_octu,
    quint_octu,
    pri_nov,
    seg_nov,
    terc_nov,
    cuar_nov,
    pri_dici,
    seg_dici,
    terc_dici,
    cuar_dici,
    quint_dici,
    pri_ene_2021,
    seg_ene_2021,
    terc_ene_2021
    ]

    frec_and = []

    for i in range(len(total_and)):
        if total_and[i] == 0:
            frec_and.append("-")
        else:
            frec_and.append(parc_and[i]/total_and[i])


    df["frec_%s" %m] = frec_and
                    
df["n"] = total_and                    
df.to_csv ("14.cov_cluster_and%s.tsv" %countm, sep="\t", index = False)


# In[ ]:





# In[ ]:





# In[ ]:





# In[76]:


#####               ESPAÑA


# In[ ]:





# In[77]:


if args.metadata_spain:
    data = {"week" : week
       }

    df = pd.DataFrame (data, columns = ["week"])

    tabla_spall = pd.read_csv(args.metadata_spain, sep = "\t") #todas las muestras españolas

    tult_Febre= 0
    tseg_Marzo = 0
    tcuar_abril = 0
    tcuar_marzo = 0
    tquint_abril = 0
    tquint_marzo = 0
    tseg_abril = 0
    tseg_mayo = 0
    tsext_marzo = 0
    tterc_abril = 0
    tterc_marzo = 0

    tterc_mayo = 0
    tcuar_mayo = 0
    tquint_mayo = 0
    tpri_junio = 0
    tseg_junio = 0
    tterc_junio = 0
    tcuar_junio = 0
    tpri_julio = 0
    tseg_julio = 0
    tterc_julio = 0
    tcuar_julio = 0
    tquint_julio = 0
    tpri_agosto = 0
    tseg_agosto = 0
    tterc_agosto = 0
    tcuar_agosto = 0
    tpri_sept = 0
    tseg_sept = 0
    tterc_sept = 0
    tcuar_sept = 0
    tpri_octu = 0
    tseg_octu = 0
    tterc_octu = 0
    tcuar_octu = 0
    tquint_octu = 0
    tpri_nov = 0
    tseg_nov = 0
    tterc_nov = 0
    tcuar_nov = 0
    tpri_dici = 0
    tseg_dici = 0
    tterc_dici = 0
    tcuar_dici = 0
    tquint_dici = 0
    tpri_ene_2021 = 0
    tseg_ene_2021 = 0
    tterc_ene_2021 = 0


    for i in range(len(tabla_spall)):
        if 20200224 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
            tult_Febre+=1
        elif 20200302 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200308:
            tseg_Marzo+=1
        elif 20200309 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200315:
            tterc_marzo+=1        
        elif 20200316 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200322:
            tcuar_marzo+=1    
        elif 20200323 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200329:
            tquint_marzo+=1
        elif 20200330 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200405:
            tsext_marzo+=1
        elif 20200406 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200412:
            tseg_abril+=1 
        elif 20200413 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200419:
            tterc_abril+=1     
        elif 20200420 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200426:
            tcuar_abril+=1   
        elif 20200427 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200503:
            tquint_abril+=1
        elif 20200504 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200510:
            tseg_mayo+=1
        elif 20200511 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200517:
            tterc_mayo+=1
        elif 20200518 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200524:
            tcuar_mayo+=1
        elif 20200525 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200531:
            tquint_mayo+=1
        elif 20200601 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200607:
            tpri_junio+=1
        elif 20200608 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200614:
            tseg_junio+=1
        elif 20200615 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200621:
            tterc_junio+=1
        elif 20200622 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200628:
            tcuar_junio+=1
        elif 20200629 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200705:
            tpri_julio+=1
        elif 20200706 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200712:
            tseg_julio+=1
        elif 20200713 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200719:
            tterc_julio+=1
        elif 20200720 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200726:
            tcuar_julio+=1
        elif 20200727 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200802:
            tquint_julio+=1
        elif 20200803 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200809:
            tpri_agosto+=1
        elif 20200810 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200816:
            tseg_agosto+=1
        elif 20200817 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200823:
            tterc_agosto+=1
        elif 20200824 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200830:
            tcuar_agosto+=1
        elif 20200831 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200906:
            tpri_sept+=1
        elif 20200907 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200913:
            tseg_sept+=1
        elif 20200914 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200920:
            tterc_sept+=1
        elif 20200921 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20200927:
            tcuar_sept+=1
        elif 20200928 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201004:
            tpri_octu+=1
        elif 20201005 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201011:
            tseg_octu+=1
        elif 2020112 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201018:
            tterc_octu+=1
        elif 20201019 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201025:
            tcuar_octu+=1
        elif 20201026 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201101:
            tquint_octu+=1
        elif 20201102 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201108:
            tpri_nov+=1
        elif 20201109 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201115:
            tseg_nov+=1
        elif 20201116 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201122:
            tterc_nov+=1
        elif 20201123 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201129:
            tcuar_nov+=1
        elif 20201130 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201206:
            tpri_dici+=1
        elif 20201207 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201213:
            tseg_dici+=1
        elif 20201214 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201220:
            tterc_dici+=1
        elif 20201221 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20201227:
            tcuar_dici+=1
        elif 20201228 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20210103:
            tquint_dici+=1
        elif 20210104 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20210110:
            tpri_ene_2021+=1
        elif 20210111 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20210117:
            tseg_ene_2021+=1
        elif 20210118 <= int("".join(tabla_spall["Collection date"][i].split("-"))) <= 20210124:
            tterc_ene_2021+=1


    for m in muts:

        tabla_sp = pd.read_csv("%s_spain.tsv" %m, sep = "\t") #solo las que tienen la mutacion especificada


        ult_Febre= 0
        seg_Marzo = 0
        cuar_abril = 0
        cuar_marzo = 0
        quint_abril = 0
        quint_marzo = 0
        seg_abril = 0
        seg_mayo = 0
        sext_marzo = 0
        terc_abril = 0
        terc_marzo = 0

        terc_mayo = 0
        cuar_mayo = 0
        quint_mayo = 0
        pri_junio = 0
        seg_junio= 0
        terc_junio= 0
        cuar_junio= 0
        pri_julio= 0
        seg_julio= 0
        terc_julio= 0
        cuar_julio= 0
        quint_julio= 0
        pri_agosto= 0
        seg_agosto= 0
        terc_agosto= 0
        cuar_agosto= 0
        pri_sept= 0
        seg_sept= 0
        terc_sept= 0
        cuar_sept= 0
        pri_octu= 0
        seg_octu= 0
        terc_octu= 0
        cuar_octu= 0
        quint_octu= 0
        pri_nov= 0
        seg_nov= 0
        terc_nov= 0
        cuar_nov= 0
        pri_dici= 0
        seg_dici= 0
        terc_dici= 0
        cuar_dici= 0
        quint_dici= 0
        pri_ene_2021= 0
        seg_ene_2021= 0
        terc_ene_2021= 0

        for i in range(len(tabla_sp)):
            if 20200224 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
                ult_Febre+=1
            elif 20200302 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200308:
                seg_Marzo+=1
            elif 20200309 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200315:
                terc_marzo+=1        
            elif 20200316 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200322:
                cuar_marzo+=1    
            elif 20200323 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200329:
                quint_marzo+=1
            elif 20200330 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200405:
                sext_marzo+=1
            elif 20200406 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200412:
                seg_abril+=1 
            elif 20200413 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200419:
                terc_abril+=1     
            elif 20200420 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200426:
                cuar_abril+=1   
            elif 20200427 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200503:
                quint_abril+=1
            elif 20200504 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200510:
                seg_mayo+=1
            elif 20200511 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200517:
                terc_mayo+=1
            elif 20200518 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200524:
                cuar_mayo+=1
            elif 20200525 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200531:
                quint_mayo+=1
            elif 20200601 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200607:
                pri_junio+=1
            elif 20200608 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200614:
                seg_junio+=1
            elif 20200615 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200621:
                terc_junio+=1
            elif 20200622 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200628:
                cuar_junio+=1
            elif 20200629 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200705:
                pri_julio+=1
            elif 20200706 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200712:
                seg_julio+=1
            elif 20200713 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200719:
                terc_julio+=1
            elif 20200720 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200726:
                cuar_julio+=1
            elif 20200727 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200802:
                quint_julio+=1
            elif 20200803 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200809:
                pri_agosto+=1
            elif 20200810 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200816:
                seg_agosto+=1
            elif 20200817 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200823:
                terc_agosto+=1
            elif 20200824 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200830:
                cuar_agosto+=1
            elif 20200831 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200906:
                pri_sept+=1
            elif 20200907 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200913:
                seg_sept+=1
            elif 20200914 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200920:
                terc_sept+=1
            elif 20200921 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20200927:
                cuar_sept+=1
            elif 20200928 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201004:
                pri_octu+=1
            elif 20201005 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201011:
                seg_octu+=1
            elif 20200112 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201018:
                terc_octu+=1
            elif 20201019 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201025:
                cuar_octu+=1
            elif 20201026 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201101:
                quint_octu+=1
            elif 20201102 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201108:
                pri_nov+=1
            elif 20201109 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201115:
                seg_nov+=1
            elif 20201116 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201122:
                terc_nov+=1
            elif 20201123 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201129:
                cuar_nov+=1
            elif 20201130 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201206:
                pri_dici+=1
            elif 20201207 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201213:
                seg_dici+=1
            elif 20201214 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201220:
                terc_dici+=1
            elif 20201221 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20201227:
                cuar_dici+=1
            elif 20201228 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20210103:
                quint_dici+=1
            elif 20210104 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20210110:
                pri_ene_2021+=1
            elif 20210111 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20210117:
                seg_ene_2021+=1
            elif 20210118 <= int("".join(tabla_sp["Collection date"][i].split("-"))) <= 20210124:
                terc_ene_2021+=1



        parc_esp = [
            ult_Febre,
            seg_Marzo,
            cuar_abril,
            cuar_marzo,
            quint_abril,
            quint_marzo,
            seg_abril,
            seg_mayo,
            sext_marzo,
            terc_abril,
            terc_marzo,

            terc_mayo,
            cuar_mayo,
            quint_mayo,
            pri_junio,
            seg_junio,
            terc_junio,
            cuar_junio,
            pri_julio,
            seg_julio,
            terc_julio,
            cuar_julio,
            quint_julio,
            pri_agosto,
            seg_agosto,
            terc_agosto,
            cuar_agosto,
            pri_sept,
            seg_sept,
            terc_sept,
            cuar_sept,
            pri_octu,
            seg_octu,
            terc_octu,
            cuar_octu,
            quint_octu,
            pri_nov,
            seg_nov,
            terc_nov,
            cuar_nov,
            pri_dici,
            seg_dici,
            terc_dici,
            cuar_dici,
            quint_dici,
            pri_ene_2021,
            seg_ene_2021,
            terc_ene_2021
            ]

        total_esp = [
            tult_Febre,
            tseg_Marzo,
            tcuar_abril,
            tcuar_marzo,
            tquint_abril,
            tquint_marzo,
            tseg_abril,
            tseg_mayo,
            tsext_marzo,
            tterc_abril,
            tterc_marzo,

            tterc_mayo,
            tcuar_mayo,
            tquint_mayo,
            tpri_junio,
            tseg_junio,
            tterc_junio,
            tcuar_junio,
            tpri_julio,
            tseg_julio,
            tterc_julio,
            tcuar_julio,
            tquint_julio,
            tpri_agosto,
            tseg_agosto,
            tterc_agosto,
            tcuar_agosto,
            tpri_sept,
            tseg_sept,
            tterc_sept,
            tcuar_sept,
            tpri_octu,
            tseg_octu,
            tterc_octu,
            tcuar_octu,
            tquint_octu,
            tpri_nov,
            tseg_nov,
            tterc_nov,
            tcuar_nov,
            tpri_dici,
            tseg_dici,
            tterc_dici,
            tcuar_dici,
            tquint_dici,
            tpri_ene_2021,
            tseg_ene_2021,
            tterc_ene_2021
            ]

        frec_esp = []
        for i in range(len(total_esp)):
            if total_esp[i] == 0:
                frec_esp.append("-")
            else:
                frec_esp.append(parc_esp[i]/total_esp[i])

        df["frec_%s" %m] = frec_esp

if args.metadata_spain:
    df["n"] = total_esp
    df.to_csv ("cov_cluster_esp.tsv", sep="\t", index = False)


# 

# In[78]:


#tabla_spall["Collection date"][1]


# In[23]:


#date = dt.datetime(int(tabla_sp["Collection date"][1].split("-")[0]), int(tabla_sp["Collection date"][1].split("-")[1]), int(tabla_sp["Collection date"][1].split("-")[2]))


# In[33]:


#prueba = date.isocalendar()[1] #para saber el número de la semana


# In[31]:


#date


# In[70]:


#total_esp


# In[40]:


#dt.datetime.strptime(tabla_spall["Collection date"][1], "%Y-%m-%d").strftime("%b-%Y")


# In[ ]:




