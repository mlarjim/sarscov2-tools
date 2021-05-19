#!/usr/bin/env python
# coding: utf-8

# In[1]:


## extractor de datos para representar la frecuencia de una mutación respecto al tiempo por provincia. Frecuencia calculada nº muestras provincia con la mut/ nº de muestras de esa provincia en esa semana


# In[48]:


import pandas as pd
import glob as g
import argparse
from collections import Counter
import re

# In[49]:


# tiene que leer:
# metadata
# las mutaciones por muestra
# el control de calidad de nextstrain


# In[1]:


parser = argparse.ArgumentParser(description = "Devuelve una tabla con la frecuencia de una determinada mutación por semana por provincia. Issue: solo alcanza hasta fechas de la tercera semana de enero 2021. Issue2: solo para aminoacidos de la spike")
parser.add_argument("metadata", help= "tsv with (at least) columns 'Strain' 'date' and 'Province'")
parser.add_argument("directory1", help= "'/complete/directory/to/*_analysis_report.csv'")
#parser.add_argument("directory2", help = "'/complete/directory/to/analysis_MN908947/'")
parser.add_argument("position", help = "aminoacid position of protein s")
parser.add_argument("-n", action = "store_true", help = "devuelve tabla con el numero de muestras total por semana por provincia")
args = parser.parse_args()


# In[50]:


#tabla = pd.read_csv("metadata.tsv", sep="\t")
tabla = pd.read_csv(args.metadata, sep = "\t")


# In[51]:


#qc = g.glob("/home/marialara/gattaca/all_run_QC_results/*_analysis_report.csv")
qc = g.glob(args.directory1)


# In[52]:


#text = g.glob("/home/marialara/gattaca/gtc/*.tsv")
#text = g.glob(args.directory2)
mn = "../covid_gattaca_ola1/"


# In[ ]:


aa = args.position

print ("Analisis de las mutaciones en la posición %s de la spike por provincia andaluza por semana" %aa)

# In[53]:


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

print("%s secuencias han pasado el control de calidad" %total)
# In[54]:

countm = 0
#sasm = 0
#metaivar =[]
#meta = []
for i in range(len(tabla)):
    for y in yes:
        id = g.glob(y)[0].split("/")[-1].split(".")[0]
        if bool(re.match("HUVR_(\d*)UK", id)) == False:
            countm = countm +1 
               
               # if tabla["Source"][i] == "SAS":
                    #sasm += 1
                    
#seto = set(meta) - set(metaivar)
#seto = list(seto)
#for i in seto:
#    print(i)
#quit()
print ("De las muestras de la primera ola del cluster, %s han pasado el control de calidad y tienen metadatos" %countm) 
#print ("De las cuales %s son del SAS" %sasm)            

primola = 0
samples = [] #muestras con la mutacion 
for a in yes:
    id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta del directorio
    #for b in yes: #recorre los id de yes
        #if b == id: #si esa muestra ha pasado el control
    #ruta = str(mn) + str(a)
    if bool(re.match("HUVR_(\d*)UK", id)) == False: #si son de la primera ola
        primola += 1
        tsv = pd.read_csv(g.glob(a)[0], sep="\t") #lee el tsv de esa muestra
        for i in range(len(tsv)):
            if  tsv["ALT_FREQ"][i]>=0.75 and int((tsv["POS"][i] - 21562 +2)/3) == int(aa) and tsv["REF_AA"][i] != tsv["ALT_AA"][i]: #si tienen mutación en ese aa
                samples.append(g.glob(a)[0].split("/")[-1].split(".")[0]) #extrae el id de la muestra de la ruta del directorio

print( "%s muestras tienen la mutacion y han pasado el control de calidad" %len(samples))

print("%s muestras de la primera ola (gisaid + sas)" %primola)
# In[57]:


ult_Febre= []
seg_Marzo = []
cuar_abril = []
cuar_marzo = []
quint_abril = []
quint_marzo = []
seg_abril = []
seg_mayo = []
sext_marzo = []
terc_abril = []
terc_marzo = []

terc_mayo = []
cuar_mayo = []
quint_mayo = []
pri_junio = []
seg_junio = []
terc_junio = []
cuar_junio = []
pri_julio = []
seg_julio = []
terc_julio = []
cuar_julio = []
quint_julio = []
pri_agosto = []
seg_agosto = []
terc_agosto = []
cuar_agosto = []
pri_sept = []
seg_sept = []
terc_sept = []
cuar_sept = []
pri_octu = []
seg_octu = []
terc_octu = []
cuar_octu = []
quint_octu = []
pri_nov = []
seg_nov = []
terc_nov = []
cuar_nov = []
pri_dici = []
seg_dici = []
terc_dici = []
cuar_dici = []
quint_dici = []
pri_ene_2021 = []
seg_ene_2021 = []
terc_ene_2021 = []


cont = 0
for i in range(len(tabla)):
    for a in samples:
        id = a.split("_")[-1]
        if bool(re.search(id, tabla["Strain"][i])) == True:
            cont += 1
            if 20200224 <= int("".join(tabla["date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
                ult_Febre.append(tabla["Province"][i])
            elif 20200302 <= int("".join(tabla["date"][i].split("-"))) <= 20200308:
                seg_Marzo.append(tabla["Province"][i])
            elif 20200309 <= int("".join(tabla["date"][i].split("-"))) <= 20200315:
                terc_marzo.append(tabla["Province"][i])        
            elif 20200316 <= int("".join(tabla["date"][i].split("-"))) <= 20200322:
                cuar_marzo.append(tabla["Province"][i])    
            elif 20200323 <= int("".join(tabla["date"][i].split("-"))) <= 20200329:
                quint_marzo.append(tabla["Province"][i]) 
            elif 20200330 <= int("".join(tabla["date"][i].split("-"))) <= 20200405:
                sext_marzo.append(tabla["Province"][i])
            elif 20200406 <= int("".join(tabla["date"][i].split("-"))) <= 20200412:
                seg_abril.append(tabla["Province"][i])  
            elif 20200413 <= int("".join(tabla["date"][i].split("-"))) <= 20200419:
                terc_abril.append(tabla["Province"][i])       
            elif 20200420 <= int("".join(tabla["date"][i].split("-"))) <= 20200426:
                cuar_abril.append(tabla["Province"][i])    
            elif 20200427 <= int("".join(tabla["date"][i].split("-"))) <= 20200503:
                quint_abril.append(tabla["Province"][i])
            elif 20200504 <= int("".join(tabla["date"][i].split("-"))) <= 20200510:
                seg_mayo.append(tabla["Province"][i])
            elif 20200511 <= int("".join(tabla["date"][i].split("-"))) <= 20200517:
                terc_mayo.append(tabla["Province"][i])
            elif 20200518 <= int("".join(tabla["date"][i].split("-"))) <= 20200524:
                cuar_mayo.append(tabla["Province"][i])
            elif 20200525 <= int("".join(tabla["date"][i].split("-"))) <= 20200531:
                quint_mayo.append(tabla["Province"][i])
            elif 20200601 <= int("".join(tabla["date"][i].split("-"))) <= 20200607:
                pri_junio.append(tabla["Province"][i])
            elif 20200608 <= int("".join(tabla["date"][i].split("-"))) <= 20200614:
                seg_junio.append(tabla["Province"][i])
            elif 20200615 <= int("".join(tabla["date"][i].split("-"))) <= 20200621:
                terc_junio.append(tabla["Province"][i])
            elif 20200622 <= int("".join(tabla["date"][i].split("-"))) <= 20200628:
                cuar_junio.append(tabla["Province"][i])
            elif 20200629 <= int("".join(tabla["date"][i].split("-"))) <= 20200705:
                pri_julio.append(tabla["Province"][i])
            elif 20200706 <= int("".join(tabla["date"][i].split("-"))) <= 20200712:
                seg_julio.append(tabla["Province"][i])
            elif 20200713 <= int("".join(tabla["date"][i].split("-"))) <= 20200719:
                terc_julio.append(tabla["Province"][i])
            elif 20200720 <= int("".join(tabla["date"][i].split("-"))) <= 20200726:
                cuar_julio.append(tabla["Province"][i])
            elif 20200727 <= int("".join(tabla["date"][i].split("-"))) <= 20200802:
                quint_julio.append(tabla["Province"][i])
            elif 20200803 <= int("".join(tabla["date"][i].split("-"))) <= 20200809:
                pri_agosto.append(tabla["Province"][i])
            elif 20200810 <= int("".join(tabla["date"][i].split("-"))) <= 20200816:
                seg_agosto.append(tabla["Province"][i])
            elif 20200817 <= int("".join(tabla["date"][i].split("-"))) <= 20200823:
                terc_agosto.append(tabla["Province"][i])
            elif 20200824 <= int("".join(tabla["date"][i].split("-"))) <= 20200830:
                cuar_agosto.append(tabla["Province"][i])
            elif 20200831 <= int("".join(tabla["date"][i].split("-"))) <= 20200906:
                pri_sept.append(tabla["Province"][i])
            elif 20200907 <= int("".join(tabla["date"][i].split("-"))) <= 20200913:
                seg_sept.append(tabla["Province"][i])
            elif 20200914 <= int("".join(tabla["date"][i].split("-"))) <= 20200920:
                terc_sept.append(tabla["Province"][i])
            elif 20200921 <= int("".join(tabla["date"][i].split("-"))) <= 20200927:
                cuar_sept.append(tabla["Province"][i])
            elif 20200928 <= int("".join(tabla["date"][i].split("-"))) <= 20201004:
                pri_octu.append(tabla["Province"][i])
            elif 20201005 <= int("".join(tabla["date"][i].split("-"))) <= 20201011:
                seg_octu.append(tabla["Province"][i])
            elif 20201012 <= int("".join(tabla["date"][i].split("-"))) <= 20201018:
                terc_octu.append(tabla["Province"][i])
            elif 20201019 <= int("".join(tabla["date"][i].split("-"))) <= 20201025:
                cuar_octu.append(tabla["Province"][i])
            elif 20201026 <= int("".join(tabla["date"][i].split("-"))) <= 20201101:
                quint_octu.append(tabla["Province"][i])
            elif 20201102 <= int("".join(tabla["date"][i].split("-"))) <= 20201108:
                pri_nov.append(tabla["Province"][i])
            elif 20201109 <= int("".join(tabla["date"][i].split("-"))) <= 20201115:
                seg_nov.append(tabla["Province"][i])
            elif 20201116 <= int("".join(tabla["date"][i].split("-"))) <= 20201122:
                terc_nov.append(tabla["Province"][i])
            elif 20201123 <= int("".join(tabla["date"][i].split("-"))) <= 20201129:
                cuar_nov.append(tabla["Province"][i])
            elif 20201130 <= int("".join(tabla["date"][i].split("-"))) <= 20201206:
                pri_dici.append(tabla["Province"][i])
            elif 20201207 <= int("".join(tabla["date"][i].split("-"))) <= 20201213:
                seg_dici.append(tabla["Province"][i])
            elif 20201214 <= int("".join(tabla["date"][i].split("-"))) <= 20201220:
                terc_dici.append(tabla["Province"][i])
            elif 20201221 <= int("".join(tabla["date"][i].split("-"))) <= 20201227:
                cuar_dici.append(tabla["Province"][i])
            elif 20201228 <= int("".join(tabla["date"][i].split("-"))) <= 20210103:
                quint_dici.append(tabla["Province"][i])
            elif 20210104 <= int("".join(tabla["date"][i].split("-"))) <= 20210110:
                pri_ene_2021.append(tabla["Province"][i])
            elif 20210111 <= int("".join(tabla["date"][i].split("-"))) <= 20210117:
                seg_ene_2021.append(tabla["Province"][i])
            elif 20210118 <= int("".join(tabla["date"][i].split("-"))) <= 20210124:
                terc_ene_2021.append(tabla["Province"][i])

    
print("%s muestras de calidad con la mutacion tienen metadatos" %cont)

# In[76]:


cad = []
mal = []
gra = []
sev = []
cor = []
jae = []
alm = []
hue = []


# In[77]:


week = [ult_Febre,
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
        
        
        


# In[78]:


pr = Counter(week[1])


# In[79]:


pr


# In[80]:


pr["Málaga"]


# In[81]:


for j in week:
    count = Counter(j)
    alm.append(count["Almería"])
    cad.append(count["Cádiz"])
    cor.append(count["Córdoba"])
    gra.append(count["Granada"])
    hue.append(count["Huelva"])
    jae.append(count["Jaén"])
    mal.append(count["Málaga"])
    sev.append(count["Sevilla"])

    
    


# In[84]:


sev


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[85]:


ult_Febre= []
seg_Marzo = []
cuar_abril = []
cuar_marzo = []
quint_abril = []
quint_marzo = []
seg_abril = []
seg_mayo = []
sext_marzo = []
terc_abril = []
terc_marzo = []

terc_mayo = []
cuar_mayo = []
quint_mayo = []
pri_junio = []
seg_junio = []
terc_junio = []
cuar_junio = []
pri_julio = []
seg_julio = []
terc_julio = []
cuar_julio = []
quint_julio = []
pri_agosto = []
seg_agosto = []
terc_agosto = []
cuar_agosto = []
pri_sept = []
seg_sept = []
terc_sept = []
cuar_sept = []
pri_octu = []
seg_octu = []
terc_octu = []
cuar_octu = []
quint_octu = []
pri_nov = []
seg_nov = []
terc_nov = []
cuar_nov = []
pri_dici = []
seg_dici = []
terc_dici = []
cuar_dici = []
quint_dici = []
pri_ene_2021 = []
seg_ene_2021 = []
terc_ene_2021 = []

print("En la tabla input hay metadatos para %s secuencias" %len(tabla))

#sas = 0
for i in range(len(tabla)):
    #if tabla["Source"][i]=="SAS":
        #sas +=1 
        if 20200224 <= int("".join(tabla["date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
            ult_Febre.append(tabla["Province"][i])
        elif 20200302 <= int("".join(tabla["date"][i].split("-"))) <= 20200308:
            seg_Marzo.append(tabla["Province"][i])
        elif 20200309 <= int("".join(tabla["date"][i].split("-"))) <= 20200315:
            terc_marzo.append(tabla["Province"][i])        
        elif 20200316 <= int("".join(tabla["date"][i].split("-"))) <= 20200322:
            cuar_marzo.append(tabla["Province"][i])    
        elif 20200323 <= int("".join(tabla["date"][i].split("-"))) <= 20200329:
            quint_marzo.append(tabla["Province"][i]) 
        elif 20200330 <= int("".join(tabla["date"][i].split("-"))) <= 20200405:
	        sext_marzo.append(tabla["Province"][i])
        elif 20200406 <= int("".join(tabla["date"][i].split("-"))) <= 20200412:
	        seg_abril.append(tabla["Province"][i])  
        elif 20200413 <= int("".join(tabla["date"][i].split("-"))) <= 20200419:
	        terc_abril.append(tabla["Province"][i])       
        elif 20200420 <= int("".join(tabla["date"][i].split("-"))) <= 20200426:
	        cuar_abril.append(tabla["Province"][i])    
        elif 20200427 <= int("".join(tabla["date"][i].split("-"))) <= 20200503:
	        quint_abril.append(tabla["Province"][i])
        elif 20200504 <= int("".join(tabla["date"][i].split("-"))) <= 20200510:
            seg_mayo.append(tabla["Province"][i])
        elif 20200511 <= int("".join(tabla["date"][i].split("-"))) <= 20200517:
            terc_mayo.append(tabla["Province"][i])
        elif 20200518 <= int("".join(tabla["date"][i].split("-"))) <= 20200524:
            cuar_mayo.append(tabla["Province"][i])
        elif 20200525 <= int("".join(tabla["date"][i].split("-"))) <= 20200531:
            quint_mayo.append(tabla["Province"][i])
        elif 20200601 <= int("".join(tabla["date"][i].split("-"))) <= 20200607:
            pri_junio.append(tabla["Province"][i])
        elif 20200608 <= int("".join(tabla["date"][i].split("-"))) <= 20200614:
            seg_junio.append(tabla["Province"][i])
        elif 20200615 <= int("".join(tabla["date"][i].split("-"))) <= 20200621:
            terc_junio.append(tabla["Province"][i])
        elif 20200622 <= int("".join(tabla["date"][i].split("-"))) <= 20200628:
            cuar_junio.append(tabla["Province"][i])
        elif 20200629 <= int("".join(tabla["date"][i].split("-"))) <= 20200705:
            pri_julio.append(tabla["Province"][i])
        elif 20200706 <= int("".join(tabla["date"][i].split("-"))) <= 20200712:
            seg_julio.append(tabla["Province"][i])
        elif 20200713 <= int("".join(tabla["date"][i].split("-"))) <= 20200719:
            terc_julio.append(tabla["Province"][i])
        elif 20200720 <= int("".join(tabla["date"][i].split("-"))) <= 20200726:
            cuar_julio.append(tabla["Province"][i])
        elif 20200727 <= int("".join(tabla["date"][i].split("-"))) <= 20200802:
            quint_julio.append(tabla["Province"][i])
        elif 20200803 <= int("".join(tabla["date"][i].split("-"))) <= 20200809:
            pri_agosto.append(tabla["Province"][i])
        elif 20200810 <= int("".join(tabla["date"][i].split("-"))) <= 20200816:
            seg_agosto.append(tabla["Province"][i])
        elif 20200817 <= int("".join(tabla["date"][i].split("-"))) <= 20200823:
            terc_agosto.append(tabla["Province"][i])
        elif 20200824 <= int("".join(tabla["date"][i].split("-"))) <= 20200830:
            cuar_agosto.append(tabla["Province"][i])
        elif 20200831 <= int("".join(tabla["date"][i].split("-"))) <= 20200906:
            pri_sept.append(tabla["Province"][i])
        elif 20200907 <= int("".join(tabla["date"][i].split("-"))) <= 20200913:
            seg_sept.append(tabla["Province"][i])
        elif 20200914 <= int("".join(tabla["date"][i].split("-"))) <= 20200920:
            terc_sept.append(tabla["Province"][i])
        elif 20200921 <= int("".join(tabla["date"][i].split("-"))) <= 20200927:
            cuar_sept.append(tabla["Province"][i])
        elif 20200928 <= int("".join(tabla["date"][i].split("-"))) <= 20201004:
            pri_octu.append(tabla["Province"][i])
        elif 20201005 <= int("".join(tabla["date"][i].split("-"))) <= 20201011:
            seg_octu.append(tabla["Province"][i])
        elif 20201012 <= int("".join(tabla["date"][i].split("-"))) <= 20201018:
            terc_octu.append(tabla["Province"][i])
        elif 20201019 <= int("".join(tabla["date"][i].split("-"))) <= 20201025:
            cuar_octu.append(tabla["Province"][i])
        elif 20201026 <= int("".join(tabla["date"][i].split("-"))) <= 20201101:
            quint_octu.append(tabla["Province"][i])
        elif 20201102 <= int("".join(tabla["date"][i].split("-"))) <= 20201108:
            pri_nov.append(tabla["Province"][i])
        elif 20201109 <= int("".join(tabla["date"][i].split("-"))) <= 20201115:
            seg_nov.append(tabla["Province"][i])
        elif 20201116 <= int("".join(tabla["date"][i].split("-"))) <= 20201122:
            terc_nov.append(tabla["Province"][i])
        elif 20201123 <= int("".join(tabla["date"][i].split("-"))) <= 20201129:
            cuar_nov.append(tabla["Province"][i])
        elif 20201130 <= int("".join(tabla["date"][i].split("-"))) <= 20201206:
            pri_dici.append(tabla["Province"][i])
        elif 20201207 <= int("".join(tabla["date"][i].split("-"))) <= 20201213:
            seg_dici.append(tabla["Province"][i])
        elif 20201214 <= int("".join(tabla["date"][i].split("-"))) <= 20201220:
            terc_dici.append(tabla["Province"][i])
        elif 20201221 <= int("".join(tabla["date"][i].split("-"))) <= 20201227:
            cuar_dici.append(tabla["Province"][i])
        elif 20201228 <= int("".join(tabla["date"][i].split("-"))) <= 20210103:
            quint_dici.append(tabla["Province"][i])
        elif 20210104 <= int("".join(tabla["date"][i].split("-"))) <= 20210110:
            pri_ene_2021.append(tabla["Province"][i])
        elif 20210111 <= int("".join(tabla["date"][i].split("-"))) <= 20210117:
            seg_ene_2021.append(tabla["Province"][i])
        elif 20210118 <= int("".join(tabla["date"][i].split("-"))) <= 20210124:
            terc_ene_2021.append(tabla["Province"][i])

#print("%s muestras de la tabla de metadatos son del SAS" %sas)
# In[86]:


ult_Febre


# In[87]:


week = [ult_Febre,
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
   

tcad = []
tmal = []
tgra = []
tsev = []
tcor = []
tjae = []
talm = []
thue = []


for j in week:
    count = Counter(j)
    talm.append(count["Almería"])
    tcad.append(count["Cádiz"])
    tcor.append(count["Córdoba"])
    tgra.append(count["Granada"])
    thue.append(count["Huelva"])
    tjae.append(count["Jaén"])
    tmal.append(count["Málaga"])
    tsev.append(count["Sevilla"])


# In[89]:


tmal


# In[90]:


mal


# In[104]:


fcad = []
fmal = []
fgra = []
fsev = []
fcor = []
fjae = []
falm = []
fhue = []


# In[105]:


for t in range(len(tcad)):
    if tcad[t] == 0:
        fcad.append("-")
    else:
        fcad.append(str(cad[t]/tcad[t]))

for t in range(len(tmal)):
    if tmal[t] == 0:
        fmal.append("-")
    else:
        fmal.append(str(mal[t]/tmal[t]))
        
for t in range(len(tgra)):
    if tgra[t] == 0:
        fgra.append("-")
    else:
        fgra.append(str(gra[t]/tgra[t]))
        
for t in range(len(tsev)):
    if tsev[t] == 0:
        fsev.append("-")
    else:
        fsev.append(str(sev[t]/tsev[t]))
        
for t in range(len(tcor)):
    if tcor[t] == 0:
        fcor.append("-")
    else:
        fcor.append(str(cor[t]/tcor[t]))
    
for t in range(len(tjae)):
    if tjae[t] == 0:
        fjae.append("-")
    else:
        fjae.append(str(jae[t]/tjae[t]))
    
for t in range(len(talm)):
    if talm[t] == 0:
        falm.append("-")
    else:
        falm.append(str(alm[t]/talm[t]))
        
for t in range(len(thue)):
    if thue[t] == 0:
        fhue.append("-")
    else:
        fhue.append(str(hue[t]/thue[t]))


# In[107]:



# In[109]:


data = {"week" : ["ult_Febre",
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
],
"almeria" : falm,
"cadiz" : fcad,
"cordoba" : fcor,
"granada" : fgra,
"huelva" : fhue,
"jaen" : fjae,
"malaga" : fmal,
"sevilla" : fsev

}


# In[110]:


df = pd.DataFrame(data, columns =  [ "week", "almeria", "cadiz", "cordoba", "granada", "huelva", "jaen", "malaga", "sevilla" ])
df.to_csv("13.%(num)scov_provincias_spike%(pos)s.csv" %{"num" : countm, "pos": aa}, sep="\t", index = False)

if args.n:
	data = {"week" : ["ult_Febre",
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
	],
	"almeria" : talm,
	"cadiz" : tcad,
	"cordoba" : tcor,
	"granada" : tgra,
	"huelva" : thue,
	"jaen" : tjae,
	"malaga" : tmal,
	"sevilla" : tsev

	}


# In[110]:


	df = pd.DataFrame(data, columns =  [ "week", "almeria", "cadiz", "cordoba", "granada", "huelva", "jaen", "malaga", "sevilla" ])
	df.to_csv("13.%(num)scov_provincias_nspike%(pos)s.csv" %{"num" : cont, "pos": aa}, sep="\t", index = False)

# In[111]:


hue


# In[112]:


thue


# In[ ]:




