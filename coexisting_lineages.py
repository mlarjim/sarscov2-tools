#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Lineages coexistentes por semana


# In[2]:


import pandas as pd
from collections import Counter
import argparse
import datetime as dt


# In[51]:


parser = argparse.ArgumentParser(description = "Devuelve una tabla con la frecuencia de los lineages por semana. Issue: solo alcanza hasta las primeras semanas de enero 2021")
parser.add_argument("Metadata", help= "tsv with (at least) columns 'Lineage' and 'date'")
parser.add_argument("-f", action = "store_true", help = "Devuelve tsv con los linajes y la fecha en la que se reportaron por primera vez")
parser.add_argument("-p", action = "store_true", help = "Devuelve tsv con los linajes por provincia")
args = parser.parse_args()


# In[3]:


tabla = pd.read_csv(args.Metadata, sep="\t")


# In[4]:


tabla=pd.DataFrame(tabla)


# In[5]:


tabla


# In[6]:


"".join(tabla["date"][0].split("-"))


# In[7]:


tabla["date"][0]


# In[8]:



num = len(tabla)

# In[9]:


lineages = []
sinlin = 0
for i in range(len(tabla)):
    lineages.append(tabla["Lineage"][i])
    if pd.isnull(tabla["Lineage"][i]): 
        sinlin +=1 #cuenta las secuencias que no tienen linaje asignado

total = num - sinlin

# In[10]:


len(lineages)


# In[11]:


lineages = list(set(lineages))
#print(lineages)

# In[12]:


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


# In[13]:


for i in range(len(tabla)):
    if 20200224 <= int("".join(tabla["date"][i].split("-"))) <= 20200301: #se transforma la fecha en formato aaaa-mm-dd en aaaammdd
        ult_Febre.append(tabla["Lineage"][i])
    elif 20200302 <= int("".join(tabla["date"][i].split("-"))) <= 20200308:
        seg_Marzo.append(tabla["Lineage"][i])
    elif 20200309 <= int("".join(tabla["date"][i].split("-"))) <= 20200315:
        terc_marzo.append(tabla["Lineage"][i])        
    elif 20200316 <= int("".join(tabla["date"][i].split("-"))) <= 20200322:
        cuar_marzo.append(tabla["Lineage"][i])    
    elif 20200323 <= int("".join(tabla["date"][i].split("-"))) <= 20200329:
        quint_marzo.append(tabla["Lineage"][i]) 
    elif 20200330 <= int("".join(tabla["date"][i].split("-"))) <= 20200405:
        sext_marzo.append(tabla["Lineage"][i])
    elif 20200406 <= int("".join(tabla["date"][i].split("-"))) <= 20200412:
        seg_abril.append(tabla["Lineage"][i])  
    elif 20200413 <= int("".join(tabla["date"][i].split("-"))) <= 20200419:
        terc_abril.append(tabla["Lineage"][i])       
    elif 20200420 <= int("".join(tabla["date"][i].split("-"))) <= 20200426:
        cuar_abril.append(tabla["Lineage"][i])    
    elif 20200427 <= int("".join(tabla["date"][i].split("-"))) <= 20200503:
        quint_abril.append(tabla["Lineage"][i])
    elif 20200504 <= int("".join(tabla["date"][i].split("-"))) <= 20200510:
        seg_mayo.append(tabla["Lineage"][i])
    elif 20200511 <= int("".join(tabla["date"][i].split("-"))) <= 20200517:
        terc_mayo.append(tabla["Lineage"][i])
    elif 20200518 <= int("".join(tabla["date"][i].split("-"))) <= 20200524:
        cuar_mayo.append(tabla["Lineage"][i])
    elif 20200525 <= int("".join(tabla["date"][i].split("-"))) <= 20200531:
        quint_mayo.append(tabla["Lineage"][i])
    elif 20200601 <= int("".join(tabla["date"][i].split("-"))) <= 20200607:
        pri_junio.append(tabla["Lineage"][i])
    elif 20200608 <= int("".join(tabla["date"][i].split("-"))) <= 20200614:
        seg_junio.append(tabla["Lineage"][i])
    elif 20200615 <= int("".join(tabla["date"][i].split("-"))) <= 20200621:
        terc_junio.append(tabla["Lineage"][i])
    elif 20200622 <= int("".join(tabla["date"][i].split("-"))) <= 20200628:
        cuar_junio.append(tabla["Lineage"][i])
    elif 20200629 <= int("".join(tabla["date"][i].split("-"))) <= 20200705:
        pri_julio.append(tabla["Lineage"][i])
    elif 20200706 <= int("".join(tabla["date"][i].split("-"))) <= 20200712:
        seg_julio.append(tabla["Lineage"][i])
    elif 20200713 <= int("".join(tabla["date"][i].split("-"))) <= 20200719:
        terc_julio.append(tabla["Lineage"][i])
    elif 20200720 <= int("".join(tabla["date"][i].split("-"))) <= 20200726:
        cuar_julio.append(tabla["Lineage"][i])
    elif 20200727 <= int("".join(tabla["date"][i].split("-"))) <= 20200802:
        quint_julio.append(tabla["Lineage"][i])
    elif 20200803 <= int("".join(tabla["date"][i].split("-"))) <= 20200809:
        pri_agosto.append(tabla["Lineage"][i])
    elif 20200810 <= int("".join(tabla["date"][i].split("-"))) <= 20200816:
        seg_agosto.append(tabla["Lineage"][i])
    elif 20200817 <= int("".join(tabla["date"][i].split("-"))) <= 20200823:
        terc_agosto.append(tabla["Lineage"][i])
    elif 20200824 <= int("".join(tabla["date"][i].split("-"))) <= 20200830:
        cuar_agosto.append(tabla["Lineage"][i])
    elif 20200831 <= int("".join(tabla["date"][i].split("-"))) <= 20200906:
        pri_sept.append(tabla["Lineage"][i])
    elif 20200907 <= int("".join(tabla["date"][i].split("-"))) <= 20200913:
        seg_sept.append(tabla["Lineage"][i])
    elif 20200914 <= int("".join(tabla["date"][i].split("-"))) <= 20200920:
        terc_sept.append(tabla["Lineage"][i])
    elif 20200921 <= int("".join(tabla["date"][i].split("-"))) <= 20200927:
        cuar_sept.append(tabla["Lineage"][i])
    elif 20200928 <= int("".join(tabla["date"][i].split("-"))) <= 20201004:
        pri_octu.append(tabla["Lineage"][i])
    elif 20201005 <= int("".join(tabla["date"][i].split("-"))) <= 20201011:
        seg_octu.append(tabla["Lineage"][i])
    elif 20201012 <= int("".join(tabla["date"][i].split("-"))) <= 20201018:
        terc_octu.append(tabla["Lineage"][i])
    elif 20201019 <= int("".join(tabla["date"][i].split("-"))) <= 20201025:
        cuar_octu.append(tabla["Lineage"][i])
    elif 20201026 <= int("".join(tabla["date"][i].split("-"))) <= 20201101:
        quint_octu.append(tabla["Lineage"][i])
    elif 20201102 <= int("".join(tabla["date"][i].split("-"))) <= 20201108:
        pri_nov.append(tabla["Lineage"][i])
    elif 20201109 <= int("".join(tabla["date"][i].split("-"))) <= 20201115:
        seg_nov.append(tabla["Lineage"][i])
    elif 20201116 <= int("".join(tabla["date"][i].split("-"))) <= 20201122:
        terc_nov.append(tabla["Lineage"][i])
    elif 20201123 <= int("".join(tabla["date"][i].split("-"))) <= 20201129:
        cuar_nov.append(tabla["Lineage"][i])
    elif 20201130 <= int("".join(tabla["date"][i].split("-"))) <= 20201206:
        pri_dici.append(tabla["Lineage"][i])
    elif 20201207 <= int("".join(tabla["date"][i].split("-"))) <= 20201213:
        seg_dici.append(tabla["Lineage"][i])
    elif 20201214 <= int("".join(tabla["date"][i].split("-"))) <= 20201220:
        terc_dici.append(tabla["Lineage"][i])
    elif 20201221 <= int("".join(tabla["date"][i].split("-"))) <= 20201227:
        cuar_dici.append(tabla["Lineage"][i])
    elif 20201228 <= int("".join(tabla["date"][i].split("-"))) <= 20210103:
        quint_dici.append(tabla["Lineage"][i])
    elif 20210104 <= int("".join(tabla["date"][i].split("-"))) <= 20210110:
        pri_ene_2021.append(tabla["Lineage"][i])
    elif 20210111 <= int("".join(tabla["date"][i].split("-"))) <= 20210117:
        seg_ene_2021.append(tabla["Lineage"][i])
    elif 20210118 <= int("".join(tabla["date"][i].split("-"))) <= 20210124:
        terc_ene_2021.append(tabla["Lineage"][i])

    


# In[14]:


len(cuar_abril)


# In[15]:


len(cuar_abril)


# In[1]:





# In[ ]:





# In[22]:


macro = []
for i in lineages:
    macro.append([Counter(ult_Febre)[i], 
                   Counter(seg_Marzo)[i],
                  Counter(terc_marzo)[i],
                   Counter(cuar_marzo)[i], 
                   Counter(quint_marzo)[i], 
                   Counter(sext_marzo)[i],
                   Counter(seg_abril)[i],
                   Counter(terc_abril)[i],
                   Counter(cuar_abril)[i],
                   Counter(quint_abril)[i],
                   Counter(seg_mayo)[i],
                   Counter(terc_mayo)[i],
                   Counter(cuar_mayo)[i],
                   Counter(quint_mayo)[i],
                   Counter(pri_junio)[i],
                  Counter(seg_junio)[i],
                    Counter(terc_junio)[i],
                    Counter(cuar_junio)[i],
                    Counter(pri_julio)[i],
                    Counter(seg_julio)[i],
                    Counter(terc_julio)[i],
                    Counter(cuar_julio)[i],
                    Counter(quint_julio)[i],
                    Counter(pri_agosto)[i],
                    Counter(seg_agosto)[i],
                    Counter(terc_agosto)[i],
                    Counter(cuar_agosto)[i],
                    Counter(pri_sept)[i],
                    Counter(seg_sept)[i],
                    Counter(terc_sept)[i],
                    Counter(cuar_sept)[i],
                    Counter(pri_octu)[i],
                    Counter(seg_octu)[i],
                    Counter(terc_octu)[i],
                    Counter(cuar_octu)[i],
                    Counter(quint_octu)[i],
                    Counter(pri_nov)[i],
                    Counter(seg_nov)[i],
                    Counter(terc_nov)[i],
                    Counter(cuar_nov)[i],
                    Counter(pri_dici)[i],
                    Counter(seg_dici)[i],
                    Counter(terc_dici)[i],
                    Counter(cuar_dici)[i],
                    Counter(quint_dici)[i],
                    Counter(pri_ene_2021)[i],
                    Counter(seg_ene_2021)[i],
                    Counter(terc_ene_2021)[i]
                  
                 ])
    


# In[ ]:





# In[24]:


lineages[0]


# In[46]:


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
]
}
                  


# In[ ]:





# In[47]:


df = pd.DataFrame(data, columns =  [ "week" ])
df.to_csv("6.coexisting_lineages_%s.csv" %total, sep="\t", index = False)


# In[49]:


for i in range(len(lineages)):
    df[lineages[i]] = macro [i]
df.to_csv("6.coexisting_lineages_%s.csv" %total, sep="\t", index = False)


# In[ ]:

firsttime = []
if args.f:
    for l in lineages: 
        fecha = []
        for i in range(len(tabla)):
            if tabla["Lineage"][i] == l:
                #fecha.append(int("".join(tabla["date"][i].split("-"))))
                fecha.append(tabla["date"][i])
        if len(fecha) == 0:
            firsttime.append("-")
        else:
            firsttime.append(min(fecha))
    data = {"Lineage" : lineages,
            "Primera_vez": firsttime
                }
    df = pd.DataFrame (data, columns = ["Lineage", "Primera_vez"])
    df.to_csv("6.first_time_%s.tsv" %total, sep="\t", index = False)
    
    
        
if args.p:
    data = { "Provincia" : ["Almeria",
    "Cadiz",
    "Cordoba",
    "Granada",
    "Huelva",
    "Jaen",
    "Malaga",
    "Sevilla" ]
    }
    df = pd.DataFrame(data, columns =  [ "Provincia" ])
    
    alm = []
    cad = []
    cor = []
    hue = []
    gra = []
    jae = []
    mal = []
    sev = []
    for i in range(len(tabla)):
        if tabla["Province"][i] == "Almería":
            alm.append(tabla["Lineage"][i])
        if tabla["Province"][i] == "Cádiz":
            cad.append(tabla["Lineage"][i])
        if tabla["Province"][i] == "Córdoba":
            cor.append(tabla["Lineage"][i])
        if tabla["Province"][i] == "Granada":
            gra.append(tabla["Lineage"][i]) 
        if tabla["Province"][i] == "Huelva":
            hue.append(tabla["Lineage"][i]) 
        if tabla["Province"][i] == "Jaén":
            jae.append(tabla["Lineage"][i]) 
        if tabla["Province"][i] == "Málaga":
            mal.append(tabla["Lineage"][i]) 
        if tabla["Province"][i] == "Sevilla":
            sev.append(tabla["Lineage"][i]) 
    prov = []
    for l in lineages:
        prov.append([Counter(alm)[l], 
                    Counter(cad)[l],
                    Counter(cor)[l],
                    Counter(gra)[l], 
                    Counter(hue)[l], 
                    Counter(jae)[l],
                    Counter(mal)[l],
                    Counter(sev)[l]
                    ] ) 
    for i in range(len(lineages)):
        df[lineages[i]] = prov [i]
    df.to_csv("6.lineages_prov_%s.csv" %total, sep="\t", index = False)        
                
            
        

# In[ ]:




