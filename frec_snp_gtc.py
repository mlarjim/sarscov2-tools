#!/usr/bin/env python
# coding: utf-8

# In[34]:


## Tabla mutaciones por muestra (solo muestras con el yes de nextstrain) y tabla frecuencia de las mutaciones


# In[35]:


import pandas as pd
import os
from collections import Counter
import glob as g
import argparse
import re


# In[36]:


parser = argparse.ArgumentParser(description = "Extractor de mutaciones con una frecuencia >=0.75 de los tsv de iVar de las muestras que han pasado el control de calidad de Nextstrain. Devuelve, por cada tsv, un csv con el identificador de la muestra como nombre, y contiene las mutaciones en formato de nucleótido. Devuelve estos csv fusionados en una tabla. Devuelve un tsv con la frecuencia de las mutaciones en la poblacion de muestras")
parser.add_argument("Directory1", help= "'/complete/path/to/*_analysis_report.csv'")
#parser.add_argument("Directory2", help = "'/complete/path/to/variants/ivar/*.tsv'")
parser.add_argument("-t", action = "store_true", help = "Devuelve un tsv con el recuento de los tipos de sustituciones de nucleotidos")
args = parser.parse_args()


# In[37]:


#qc=open("qc.txt", "r") #qc.txt es un archivo de texto plano con el directorio de todos los report de gattaca


# In[38]:


#qc = qc.read() #para leer el txt


# In[39]:


#qc=qc.split() #dividirlo por directorios


# In[40]:


#qc = g.glob("/home/marialara/gattaca/all_run_QC_results/*_analysis_report.csv")
qc = g.glob(args.Directory1)


# In[41]:


# qc[0]
mn = "../covid_gattaca_ola1/"

# In[42]:


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
            yes.append(ap) #guarda el nombre de las muestras seleccionadas
        
gisaid = g.glob("/home/mlara/covid_gattaca_ola1/run_gisaid*/variants/ivar/*.tsv") #guarda directorio de los falsos .tsv de gisaid

yes = yes + gisaid
# In[43]:


#total = len(yes) #numero de muestras aceptadas


archivos = yes
# In[44]:


#archivos=open("anotacion.txt", "r") #anotacion.txt es un archivo de texto plano que guarda los directorios de los archivos tsv

# archivos = g.glob("/home/marialara/gattaca/gtc/*.tsv")
#archivos = g.glob(args.Directory2)


# In[ ]:





# In[ ]:





# In[46]:


#archivos= archivos.read().split() #se lee el txt y se separan los directorios
# archivos[0]

# archivos[0].split("/")[5].split(".")[0]


# In[48]:
total = 0
for a in archivos:
    id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta 
    if bool(re.match("HUVR_(\d*)UK", id)) == False: #si son de la primera ola
        total += 1

#len(archivos)
print ("%s muestras de la primera ola han pasado el control de calidad de Nextstrain"%total)  

os.mkdir("9.mut%s" %total)


# In[49]:


snps = []
for a in archivos:
    id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta del directorio
    #print (id)
    #for b in yes: #recorre los id de yes
        #if b == id: #si esa muestra ha pasado el control
    if bool(re.match("HUVR_(\d*)UK", id)) == False: #si son de la primera ola
        sample = []
        reg = []
        pos = []
        ref = []
        alt = []
        tabla = pd.read_csv(g.glob(a)[0], sep="\t") #lee el tsv
        for c in range(len(tabla)):
            if tabla["ALT_FREQ"][c] >= 0.75 : #si la mutación tiene una frecuencia mayor o igual a 0.75
                sample.append(id) 
                pos.append(tabla["POS"][c])
                ref.append(tabla["REF"][c])
                alt.append(tabla["ALT"][c])
               
                
                if len(tabla["REF"][c]) == 1 and len(tabla["ALT"][c]) == 1: #solo snps
                    snps.append(str(tabla["REF"][c]) + str(tabla["POS"][c]) + str(tabla["ALT"][c]))
                    
                if tabla["POS"][c]< 266:
                    reg.append("5UTR")
                if 266<=int(tabla["POS"][c])<=21555:
                    reg.append("ORF1ab")
                if 21556<=tabla["POS"][c]<=21562:
                    reg.append("ds_ORF1ab")
                if 21563<=tabla["POS"][c]<=25384:
                    reg.append("spike")
                if 25385<=tabla["POS"][c]<=25392:
                    reg.append("ds_spike")    
                if 25393<=tabla["POS"][c]<=26220:
                    reg.append("ORF3a")
                if 26221<=tabla["POS"][c]<=26244:
                    reg.append("ds_ORF3a")
                if 26245<=tabla["POS"][c]<=26472:
                    reg.append("E")
                if 26473<=tabla["POS"][c]<=26522:
                    reg.append("ds_E")
                if 26523<=tabla["POS"][c]<=27191:
                    reg.append("M")
                if 27192<=tabla["POS"][c]<=27201:
                    reg.append("ds_M")
                if 27202<=tabla["POS"][c]<=27387:
                    reg.append("ORF6")
                if 27388<=tabla["POS"][c]<=27393:
                    reg.append("ds_ORF6")
                if 27394<=tabla["POS"][c]<=27755:
                    reg.append("ORF7a")
                if 27756<=tabla["POS"][c]<=27759:
                    reg.append("ORF7a/ORF7b")
                if 27760<=tabla["POS"][c]<=27887:
                    reg.append("ORF7b")
                if 27888<=tabla["POS"][c]<=27893:
                    reg.append("ds_ORF7b")
                if 27894<=tabla["POS"][c]<=28259:
                    reg.append("ORF8")
                if 28260<=tabla["POS"][c]<=28273:
                    reg.append("ds_ORF8")
                if 28274<=tabla["POS"][c]<=29533:
                    reg.append("N")
                if 29534<=tabla["POS"][c]<=29557:
                    reg.append("ds_N")
                if 29558<=tabla["POS"][c]<=29674:
                    reg.append("ORF10")
                if tabla["POS"][c] > 29674:
                    reg.append("3UTR")
                
               
                
        
        data = { "SAMPLE" : sample,
            "REGION" : reg,
            "POS": pos,
            "REF": ref,
            "ALT" : alt
               }
        df = pd.DataFrame(data, columns =  ["SAMPLE", "REGION", "POS", "REF", "ALT"])
        df.to_csv(r"9.mut%(total)s/%(id)s_muts.tsv" %{ "total" : total, "id": id}, sep="\t", index=False)

                
            
            
         


# In[51]:


snps = list(Counter(snps).keys())


# In[56]:


#snps[0][0]


# In[67]:


if args.t:
    tg, tc, ta, gt, gc, ga, ct, cg, ca, at, ag, ac = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 #contadores
    for i in snps:
        if i[0] == "T" and i[-1] == "G":
            tg = tg +1 
        elif i[0] == "T" and i[-1] == "C":
            tc = tc +1
        elif i[0] == "T" and i[-1] == "A":
            ta = ta +1
        elif i[0] == "G" and i[-1] == "T":
            gt = gt +1
        elif i[0] == "G" and i[-1] == "C":
            gc = gc +1
        elif i[0] == "G" and i[-1] == "A":
            ga = ga +1
        elif i[0] == "C" and i[-1] == "T":
            ct = ct +1
        elif i[0] == "C" and i[-1] == "G":
            cg = cg +1
        elif i[0] == "C" and i[-1] == "A":
            ca = ca +1
        elif i[0] == "A" and i[-1] == "T":
            at = at +1
        elif i[0] == "A" and i[-1] == "G":
            ag = ag +1
        elif i[0] == "A" and i[-1] == "C":
            ac = ac +1

    data = { "MUT" : ["T / G", "T / C", "T / A", "G / T", "G / C", "G / A", "C / T", "C / G", "C / A", "A / T", "A / G", "A / C"],
            "COUNT" : [tg, tc, ta, gt, gc, ga, ct, cg, ca, at, ag, ac]   
    }

    df = pd.DataFrame(data, columns =  ["MUT", "COUNT"])
    df.to_csv(r"9.snp_frec_mut%s.tsv" %total, sep="\t", index=False)


# In[30]:


nombresFicheros = g.glob("9.mut%s/*_muts.tsv" %total)


# In[256]:


array = []
for i in nombresFicheros:
    df = pd.read_csv(i, sep='\t')
    array.append(df)

dfgrande = pd.concat(array)


# In[257]:


dfgrande.to_csv("9.all_mut%s.tsv" %total, sep="\t", index=False)


# In[258]:


tabla = pd.read_csv("9.all_mut%s.tsv" %total, sep="\t")


# In[259]:


muts = []
for i in range(len(tabla)):
    mut = []
    mut = str(tabla["REGION"][i]) + ":" + " " + str(tabla["POS"][i]) + "-" + str(tabla["REF"][i]) + "-" + str(tabla["ALT"][i])
    muts.append(mut)
    


# In[260]:


total


# In[261]:


c = Counter(muts)


# In[262]:


mut = list(c.keys())
frec = list(c.values())


# In[263]:


for j in range(len(frec)):
    frec[j] = frec[j]/total


# In[ ]:





# In[264]:


data = { "MUT" : mut,
        "FREC" : frec,   
}

df = pd.DataFrame(data, columns =  ["MUT", "FREC"])
df.to_csv("9.frec_mut%s.tsv" %total, sep="\t", index=False)


# In[ ]:





# In[ ]:


# array = []
# for i in nombresficheros:
#    df = pd.read_csv(i, sep='\t')
 #   df['SAMPLE'] = i
  #  array.append(df)

# dfgrande = pd.concat(array)

# dfgrande['MUT'] = dfgrande['POS'] + '-' + dfgrande['REF'] + '-' + dfgrande['ALT']

    

