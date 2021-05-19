#!/usr/bin/env python
# coding: utf-8

# In[76]:


## Calculador de coeficiente de variabilidad de Wu-Kabat de las muestras que han pasado el control de nextstrain a partir de los tsv de iVar (frec >= 0,75)


# In[4]:


import pandas as pd
import os
from collections import Counter
import glob as g
import argparse
import re


# In[16]:


parser = argparse.ArgumentParser(description = "Calculador de coeficiente de Wu-Kabat de las muestras que han pasado el control de nextstrain a partir de los tsv de iVar (frec >= 0.75). Devuelve un tsv para cada gen")
parser.add_argument("Directory1", help= "'/complete/directory/to/*_analysis_report.csv'")
#parser.add_argument("Directory2", help = "'/complete/directory/to/variants/ivar/*.tsv'")
parser.add_argument("-m", action="store_true", help = "Devuelve un tsv con las mutaciones en formato de aminoacido de todos los genes, aunque esten repetidas")
parser.add_argument("-p", action="store_true", help = "Devuelve un tsv por cada proteína, con los valores de posición de aa y variabilidad de Wu-Kabat listos para su representación gráfica")
args = parser.parse_args()



# In[9]:


#qc = g.glob("/home/marialara/gattaca/all_run_QC_results/*_analysis_report.csv")
qc = g.glob(args.Directory1)
#text = g.glob("/home/marialara/gattaca/gtc/*.tsv")
#text = g.glob(args.Directory2)


# In[7]:

#todas = len(text)
#print ("Hay %s tsvs" %todas) 





# lista=[]
# for i in text:
#    lista.append(i.split("/")[7])

    
mn = "../covid_gattaca_ola1/"

# In[40]:


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

#total = len(yes)       

text = yes
# In[41]:


#print ("%s muestras han pasado el control de calidad de Nextstrain"%total) 


# In[42]:


#lista[0]
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


# In[ ]:
total = 0
for a in text:
        id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta del directorio
        if bool(re.match("HUVR_(\d*)UK", id)) == False: #si son de la primera ola
            total += 1
os.mkdir("7.wu_kabat%s" %total)

print(total)

# In[43]:

for a in text:
        id = g.glob(a)[0].split("/")[-1].split(".")[0] #extrae el id de la muestra de la ruta del directorio
    #for b in yes: #recorre los id de yes
        #if b == id: #si esa muestra ha pasado el control
        if bool(re.match("HUVR_(\d*)UK", id)) == False: #si son de la primera ola
            tabla = pd.read_csv(g.glob(a)[0], sep="\t")



            ## ORF1ab
            ## join(266..13468,13468..21555)

            # Para convertir la posición a posición de aa, si es el primer nucleotido del codon: (nt - 265+2)/3
            # Con esa formula, el segundo nt del codon daria ,333. Para quitarle los decimales: int((nt - 265+2)/3)
            # A partir de la posicion 13468, la formula anterior no sirve. Es: ins((nt - 264+2)/3)


            for i in range(len(tabla)):
                if 266<=tabla["POS"][i]<13468 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 265 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf1ab.append(mut)
                if 13468<=tabla["POS"][i]<=21555 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 264 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf1ab.append(mut)



            ## Spike 
            ## 21563..25384

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt - 21562+2)/3)

            for i in range(len(tabla)):
                if 21563<=tabla["POS"][i]<=25384 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 21562 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    spike.append(mut)
                    #if pos == 490:
                        #print(tabla["POS"][i])
                        #print(g.glob(a)[0])




            ## ORF3a
            ## 25393..26220

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-25392+2)/3)

            for i in range(len(tabla)):
                if 25393<=tabla["POS"][i]<=26220 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 25392 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf3a.append(mut)




            ## E
            ## 26245..26472

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-26244+2)/3)

            for i in range(len(tabla)):
                if 26245<=tabla["POS"][i]<=26472 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 26244 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    e.append(mut)





            ## M
            ## 26523..27191

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-26522+2)/3)

            for i in range(len(tabla)):
                if 26523<=tabla["POS"][i]<=27191 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 26522 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    me.append(mut)

            ## orf6
            ## 27202..27387

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-27201+2)/3)

            for i in range(len(tabla)):
                if 27202<=tabla["POS"][i]<=27387 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27201 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf6.append(mut)

            ## orf7a
            ## 27394..27759

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-27393+2)/3)

            for i in range(len(tabla)):
                if 27394<=tabla["POS"][i]<=27759 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27393 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf7a.append(mut)

            ## orf7b
            ## 27756..27887

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-27755+2)/3)

            for i in range(len(tabla)):
                if 27756<=tabla["POS"][i]<=27887 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27755 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf7b.append(mut)

            ## orf8
            ## 27894..28259

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-27892+2)/3)

            for i in range(len(tabla)):
                if 27894<=tabla["POS"][i]<=28259 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 27892 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf8.append(mut)




            ## N
            ## 28274..29533

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-28273+2)/3)

            for i in range(len(tabla)):
                if 28274<=tabla["POS"][i]<=29533 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 28273 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    nu.append(mut)



            ## ORF10
            ## 29558..29674

            # Para convertir la posición a posición de aa, para cualquier nucleotido del codon: int((nt-29557+2)/3)

            for i in range(len(tabla)):
                if 29558<=tabla["POS"][i]<=29674 and tabla["ALT_FREQ"][i]>=0.75:
                    pos = int((tabla["POS"][i] - 29557 +2)/3)
                    mut = str(tabla["REF_AA"][i]) + str(pos) + str(tabla["ALT_AA"][i])
                    orf10.append(mut)




# In[44]:


if args.m :
	data = { "GENE" : ["orf1ab"]*len(orf1ab) + ["spike"]*len(spike) + ["orf3a"]*len(orf3a) + ["envelope"]*len(e) + ["M"]*len(me) + ["orf6"]*len(orf6) 
	       + ["orf7a"]*len(orf7a) + ["orf7b"]*len(orf7b) + ["orf8"]*len(orf8) + ["n"]*len(nu) + ["orf10"]*len(orf10),
		"MUT" : orf1ab + spike + orf3a + e + me + orf6 + orf7a + orf7b + orf8 + nu + orf10
		}

	df = pd.DataFrame(data, columns =  ["GENE", "MUT"])
	df.to_csv(r"7.wu_kabat%s/mut_aas_gtc.tsv" %total, sep="\t", index=False)


# In[48]:


len(orf10)


# In[49]:


# Variability=N∗k/n
# N = total number of sequences in the alignment
# k = number of different aa at a given position
# n = frequency of the most common aa at that position
name = 0
nom =["orf1ab", "spike", "orf3a", "e", "m", "orf6", "orf7a", "orf7b", "orf8", "nu", "orf10"]
superl = [orf1ab, spike, orf3a, e, me, orf6, orf7a, orf7b, orf8, nu, orf10]
lon = [ 7096, 1273, 275, 75, 222, 61, 121, 43, 121, 419, 38 ] #en realidad el orf7b parece no existir en el genoma de referencia MN908947, se ha cogido la lontigutd de NC_045512.2
#nom = ["spike"]
#superl= [spike]


N = len(yes) #numero de muestras que han pasado el control de nextstrain

for g in superl: 


    # k = number of different aa at a given position



    quitar = []
    gk=g[:] #modificacion de la lista con el objetivo de obtener k, se crea un duplicado independiente
    for i in range(len(gk)):
        if gk[i][0] == gk[i][-1]: #si es sinonima
            quitar.append(gk[i])


    for j in quitar: #se quitan los sinonimos porque no hay un cambio de nucleotido
        gk.remove(j)



    for x in range(len(gk)):
        gk[x] = gk[x][1:len(gk[x])] #se le quita el primer caracter para poder ordenarlos

    #gk.sort() #se ordenan
    c =Counter(gk) # para la frecuencia
    

    c1 = list(c.keys()) #mutaciones unficadas sin el aa de referencia
    c2 = []
    for h in range(len(c1)):
        c2.append(c1[h][0:-1]) #se le quita el aa alterado para saber el numero de aas diferentes
    c2.sort()
    ter = Counter(c2)
    c3 = list(ter.keys()) #posiciones con cambio de aa unificadas
    c4 =  list(ter.values()) #lo normal en posiciones mutadas es que k sea 2
    k = []
    for v in c4:
    	k.append(v+1) #+1 por el aa de referencia
   

    #f = 0
    #for i in range(len(c1)-1): 
    #    if c1[i][0:len(c1[i])-1] == c1[i+1][0:len(c1[i+1])-1]: #si la posicion es la misma
    #        f = f +1  #asi se sabe cuantos elementos hay que quitar
  


   # if f != 0:      
   #     for i in range(len(c1)-f-1): #se le resta el numero de elementos a eliminar
   #         if c1[i][0:len(c1[i])-1] == c1[i+1][0:len(c1[i+1])-1]:
   #             c1.pop(i+1)
   #             k[i] = k[i] + 1
   #             k.pop(i+1)
    #if len(c1)>1 and c1[-1][0:-1] == c1[-2][0:-1]:
    #    c1.pop(-1)
     #   k[i] = k[i] + 1
      #  k.pop(-1)
        
      #  22d 22a 22c 23a    f=2 range 2
 

    for p in range(len(c1)):
        c1[p] = c1[p][0:len(c1[p])-1]

        
        
        
        
        
    # n = frequency of the most common aa at that position



    gn = g[:] #duplicado independiente

    quitar = []
    for i in range(len(gn)):
        if gn[i][0] == gn[i][-1]: #si es sinonima
            quitar.append(gn[i])
    for j in quitar: #se quitan los sinonimos porque no hay un cambio de nucleotido
        gn.remove(j)
    gn2 = gn[:] #otro duplicado independiente, esta vez para determinar si el aa de referencia es o no el más frecuente


    for x in range(len(gn)):
        gn[x] = gn[x][1:len(gn[x])-1] #se le quita el primer y ultimo caracter para dejar la posicion

    #Se hace un recuento de las posiciones con un aminoacido diferente a la referencia (mutacion no sinonima)

    gn.sort() #se ordenan


    casin = list(Counter(gn).values())
    p = list(Counter(gn).keys())


    n = []
    for m in casin:
        n.append(N - m)
    
    noref=[] #por si hay un aminoacido mas frecuente que el de la referencia
    for x in range(len(gn2)):
        noref.append(gn2[x][1:]) #se le quita el primer caracter a los elementos de la lista inicial
    noref.sort() #se ordenan por numero
    cnoref = Counter(noref)
    vcnoref = list(cnoref.values()) #la frecuencia de la mutacion
    kcnoref = list(cnoref.keys()) #la posicion + el aminoacido mutado
    
    for i in range(len(p)):
        for h in range(len(vcnoref)):
            if p[i] == kcnoref[h][0:-1] and vcnoref[h] > n[i]: #si la frecuencia del aa mutado es mayor a la del aa de referencia
                n[i] = vcnoref[h]


    # Variability=N∗k/n

    wk = []

    for i in range(len(p)):
        wk.append((N*k[i])/n[i])
        

    if args.p:
        aas = [] #aqui se guardan las posiciones de todos los aas de esa proteina
        wk2 = [1] * lon[name] #lista con los valores de variabilidad para aas no mutados
        for d in range(1,lon[name]+1):
            aas.append(d)
            for f in range(len(p)): #itera sobre posiciones mutadas
                if str(d) == p[f]: #si la posicion tiene mutaciones
                    wk2[d-1] = wk[f] #d-1 porque d empieza a contar en 1
	    
        data2 = { "POS_AA" : aas,
        "Wu-Kabat" : wk2
        }
        
		     
		    
        df = pd.DataFrame(data2, columns =  ["POS_AA", "Wu-Kabat"])
        df.to_csv(r"7.wu_kabat%(total)s/plot_wk_%(prot)s.csv" %{"total": total, "prot": nom[name]}, sep="\t", index=False)
    
    
    

    data = { "POS_AA" : p,
            "n": n,
            "k": k,
            "Wu-Kabat" : wk
            
           }

    	

    df = pd.DataFrame(data, columns =  ["POS_AA", "n", "k", "Wu-Kabat"])
    df.to_csv(r"7.wu_kabat%(total)s/wk_%(prot)s.csv" %{"total": total, "prot": nom[name]}, sep="\t", index=False)
    name = name + 1


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




