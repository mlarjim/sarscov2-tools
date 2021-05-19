#!/usr/bin/env python
#coding: utf-8

# In[1]:


## Script para recrear la tabla de Zekri et al. 2020


# In[1]:


import pandas as pd
import os
from collections import Counter
import glob as g
import argparse
import re

# In[2]:

parser = argparse.ArgumentParser(description = "Analizador del tipo de mutación")
parser.add_argument("report", help= "'/complete/path/to/*_analysis_report.csv'")
#parser.add_argument("variants", help= "'/complete/path/to/variants/ivar/*.tsv'")
parser.add_argument("-a", action = "store_true", help = "incluye una columna con las alteraciones")
args = parser.parse_args()

#qc = g.glob("/home/marialara/gattaca/all_run_QC_results/*_analysis_report.csv")
qc = g.glob(args.report)


# In[3]:

mn = "../covid_gattaca_ola1/"

yes = []
for w in range(len(qc)):
    tabla = pd.read_csv(qc[w], sep=",") #lee una de las tablas (csv)
    
    for i in range(len(tabla)):
        if tabla["selected_for_nextstrain"][i] == "yes": #si está seleccionada para nextstrain
            run = qc[w].split("/")[-1].split("_")[0]	
            num = qc[w].split("/")[-1].split("_")[1]	
            tog = [run, num]
            rute = "_".join(tog)
            ap = mn + rute + "_batch*/variants/ivar/" + tabla["sample"][i] +".tsv" #ruta de muestras que han pasado el qc
            yes.append(ap) #guarda las rutas
            
gisaid = g.glob("/home/mlara/covid_gattaca_ola1/run_gisaid*/variants/ivar/*.tsv") #guarda directorio de los falsos .tsv de gisaid

yes = yes + gisaid            


# In[4]:


total = len(yes) #aqui van incluidas tambien las de la segunda ola
print("%s secuencias han pasado el control de calidad" %total)

# In[5]:


#archivos = g.glob("/home/marialara/gattaca/gtc/*.tsv")
#archivos = g.glob(args.variants)

archivos = yes
# In[49]:


utr5 = []
orf1ab= []
spike = []
orf3a= []
e = []
me = []
orf6 = []
orf7a = []
orf7b=[]
orf8 = []
nu = []
orf10 = []
utr3 = []
intergen = []
orf7ab = []

mutorf1ab= []
mutspike = []
mutorf3a= []
mute = []
mutme = []
mutorf6 = []
mutorf7a = []
mutorf7b=[]
mutorf8 = []
mutnu = []
mutorf10 = []
mutorf7ab = []

fw = 0

for a in archivos:
    id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta del directorio
    #for b in yes: #recorre los id de yes
    if bool(re.match("HUVR_(\d*)UK", id)) == False: #si es de la primera ola
        fw += 1
        #if b == id: #si esa muestra ha pasado el control
        tabla = pd.read_csv(g.glob(a)[0], sep="\t") #lee el tsv
        for i in range(len(tabla)):
           # if tabla["ALT_FREQ"][i] >= 0.75: #si la mutación tiene una frecuencia mayor o igual a 0.75
                
                
                if tabla["POS"][i]< 266 and tabla["ALT_FREQ"][i] >= 0.75:
                    utr5.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                    
                elif 266<=tabla["POS"][i]<13468 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 265 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf1ab.append(mut)
                    mutorf1ab.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 13468<=tabla["POS"][i]<=21555 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 264 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf1ab.append(mut)
                    mutorf1ab.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))    
                    
                elif 21563<=tabla["POS"][i]<=25384 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 21562 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    spike.append(mut)
                    mutspike.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))  

                elif 25393<=tabla["POS"][i]<=26220 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 25392 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf3a.append(mut)
                    mutorf3a.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 26245<=tabla["POS"][i]<=26472 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 26244 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    e.append(mut)
                    mute.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 26523<=tabla["POS"][i]<=27191 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 26522 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    me.append(mut)
                    mutme.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 27202<=tabla["POS"][i]<=27387 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27201 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf6.append(mut)
                    mutorf6.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 27394<=tabla["POS"][i]<=27759 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27393 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf7a.append(mut)
                    mutorf7a.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 27756<=tabla["POS"][i]<=27759 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27393 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf7ab.append(mut)
                    mutorf7ab.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 27756<=tabla["POS"][i]<=27887 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27755 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf7b.append(mut)
                    mutorf7b.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 27894<=tabla["POS"][i]<=28259 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27893 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf8.append(mut)
                    mutorf8.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 28274<=tabla["POS"][i]<=29533 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 28273 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    nu.append(mut)
                    mutnu.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif 29558<=tabla["POS"][i]<=29674 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 29557 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf10.append(mut)
                    mutorf10.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                    
                elif tabla["POS"][i] > 29674 and tabla["ALT_FREQ"][i]>=0.75:
                    utr3.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
                
                elif tabla["ALT_FREQ"][i]>=0.75: 
                    intergen.append(str(tabla["REF"][i]) + str(tabla["POS"][i]) + str(tabla["ALT"][i]))
            
        


# In[50]:

mutgenes = [mutorf1ab, mutspike, mutorf3a, mute, mutme, mutorf6, mutorf7a, mutorf7ab, mutorf7b, mutorf8, mutnu, mutorf10]
genes2 = [orf1ab, spike, orf3a, e, me, orf6, orf7a, orf7ab, orf7b, orf8, nu, orf10]

for g in range(len(mutgenes)): #recorre lista con las listas de mutaciones de nucleotidos
    for j in range(len(mutgenes[g])):
        for i in range(len(mutgenes[g])):
            if mutgenes[g][j] == mutgenes[g][i]: #si la mutacion es la misma a nivel ntds
                if genes2[g][j][0:3] == "non": #si alguna del gisaid no se conoce
                    genes2[g][j] = genes2[g][i] #se cambia por la conocida
                if genes2[g][i][0:3] == "non":
                    genes2[g][i] = genes2[g][j]
                
        
    

regiones = [utr5, orf1ab, spike, orf3a, e, me, orf6, orf7a, orf7ab, orf7b, orf8, nu, orf10, utr3, intergen]
orf1ab = list(set(orf1ab))
spike = list(set(spike))
orf3a = list(set(orf3a))
e = list(set(e))
me = list(set(me))
orf6 = list(set(orf6))
orf7a = list(set(orf7a))
orf7ab = list(set(orf7ab))
orf8 = list(set(orf8))
nu = list(set(nu))
orf10 = list(set(orf10))
genes = [orf1ab, spike, orf3a, e, me, orf6, orf7a, orf7ab, orf7b, orf8, nu, orf10]

#for x in genes1:
#    for m in range(len(x)):
#        for n in range(len(x)):
#            if re.findall("\d+", x[m]) == re.findall("\d+", x[n]):
#                if x[m][0:3] == "non":
#                    x[m] = x[n]
 #               if x[n][0:3] == "non":
#                    x[n] = x[m]

#genes = genes2     

# In[51]:


indels = []


for i in genes:
    count = 0
    for j in i:
        if j[0:3] == "nan":
            count = count + 1
    indels.append(count)
            


# In[52]:


indels


# In[53]:


for g in genes:
    quitar = []
    for j in g:
        if j[0:3] == "nan": #se quitan las mutaciones de inserción/deleción
            quitar.append(j)
    for k in quitar:
        g.remove(k)
        
            
    


# In[54]:


sin = []
for g in genes:
    count = 0
    for j in g: #mutaciones en ese gen
        if j[0] == j[-1] and j[0:3] != "nan": #si el aa de ref y la alteracion son el mismo
            count = count + 1
    sin.append(count)


# In[55]:


nosin = []
for g in genes:
    count = 0
    for j in g:
        if j[0] != j[-1]:
            count = count + 1
    nosin.append(count)


# In[56]:


sin


# In[57]:


names = ["orf1ab", "spike", "orf3a", "e", "me", "orf6", "orf7a", "orf7ab", "orf7b", "orf8", "nu", "orf10"]


# In[58]:


genes_c = []
for i in genes:
    genes_c.append(list(Counter(i).keys()))


# In[59]:

if args.a:
	data = { "GEN" : names,
		        "SYNONYMOUS": sin,
		        "MISSENSE": nosin,
		        "INDELS" : indels,
		        "MUTATIONS" : genes_c
		           }
	df = pd.DataFrame(data, columns =  ["GEN", "SYNONYMOUS", "MISSENSE", "INDELS", "MUTATIONS"])
	df.to_csv("10.tabla_zekri_mutaa%s.tsv" %fw, sep="\t", index=False)

else:
	data = { "GEN" : names,
		        "SYNONYMOUS": sin,
		        "MISSENSE": nosin,
		        "INDELS" : indels,
		           }
	df = pd.DataFrame(data, columns =  ["GEN", "SYNONYMOUS", "MISSENSE", "INDELS"])
	df.to_csv("10.tabla_zekri%s.tsv" %fw, sep="\t", index=False)


# In[60]:


ncr = [list(set(utr5)), list(set(utr3)), list(set(intergen))]
ncr_c = [] 


# In[61]:


for i in ncr:
    ncr_c.append(list(Counter(i).keys()))
    


# In[ ]:





# In[62]:


num = []
for i in ncr:
    num.append(len(i))
    


# In[63]:


num


# In[64]:


names_ncr = ["5UTR", "3UTR", "intergen"]


# In[65]:


data = { "REGION" : names_ncr,
                "NUM_MUT": num,
                "MUTATIONS": ncr_c
                   }
df = pd.DataFrame(data, columns =  ["REGION", "NUM_MUT", "MUTATIONS"])
df.to_csv("10.tabla_zekri_mutncr%s.tsv" %fw, sep="\t", index=False)


# In[ ]:





# In[ ]:




