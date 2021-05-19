#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import argparse


# In[2]:

parser = argparse.ArgumentParser(description = "Calculador de media y desviacion estandar de una matriz simetrica")
parser.add_argument("matrix", help = "Matriz simetrica de distancias geneticas (tsv, recomendacion: output de snp-dists)")
args=parser.parse_args()

matrix = pd.read_csv(args.matrix, sep= "\t")


# In[3]:


df = pd.DataFrame(matrix)


# In[4]:


diago = len(matrix)


# In[10]:


df.columns[0]


# In[11]:


df_1=df.drop(df.columns[0], axis=1)


# In[12]:


sumatorio = df_1.sum().sum()


# In[13]:


total = diago*diago-diago


# In[14]:
txt = open("2.average_distance%s.txt" %diago, "w")

media = sumatorio/total
print ("Distancia genetica media: %s" %media)
txt.write("Distancia genetica media: %s\n" %media)

# In[20]:


## desviacion estandar


# In[21]:


df_2=np.triu(df_1, k=1)


# In[22]:


updf=list(df_2[np.triu_indices(len(df_2), k=1)])


# In[23]:


np.mean(updf)


# In[24]:


print ("Desviacion estandar: %s" %np.std(updf))
txt.write("Desviacion estandar: %s" %np.std(updf))


# In[ ]:




