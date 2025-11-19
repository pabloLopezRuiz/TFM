import scanpy as sc
import plotly.express as ex
import pandas as pd

x = sc.read_h5ad("test.h5ad")
x.obs["Cell.type"].value_counts()
x.obs.columns
x.obs["nFeature_RNA"]
sum(x.obs["num_genes_expressed"] != x.obs["nFeature_RNA"]) # Las dos columnas son iguales

# número de conteos celulares, con neuronas marcadas
df_counts = x.obs["Cell.type"].value_counts()
is_neuron = x.obs[["Cell.type", "Tissue"]]
is_neuron = is_neuron.reset_index().drop_duplicates(subset="Cell.type")[["Cell.type", "Tissue"]]
is_neuron["is.neuron"] = is_neuron["Tissue"] == "Neuron"

df_counts_merged = pd.merge(df_counts, is_neuron, on="Cell.type")
# Ya está sorted
df_counts_merged.to_pickle("data-taylor/df_counts.pkl")

# Para el violin plot de ncount y ngenes
df_ng_nc = x.obs[["Cell.type", "nCount_RNA", "nFeature_RNA", "Tissue"]]
df_ng_nc["is.neuron"] = (df_ng_nc["Tissue"] == "Neuron")
 
df_ng_nc.to_pickle("data-taylor/df_ng_nc.pkl")


nameneuron =  pd.DataFrame(
    data = { 
        "cell_name": ,
        "is_neuron": 
    },
    
)



cell_name = x.obs["Cell.type"].reset_index()
is_neuron = (x.obs["Tissue"] == "Neuron").reset_index()

name_curated = pd.merge(cell_name, is_neuron, on="index")
name_curated.drop("index", axis = 1, inplace=True)
name_curated = name_curated[name_curated["Tissue"] == True]
print(list(name_curated["Cell.type"]))
