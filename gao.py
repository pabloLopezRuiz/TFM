import scanpy as sc
import plotly.express as ex
import plotly.graph_objects as go
import pandas as pd

x = sc.read_h5ad("daty/GSE992099.h5ad")
x.obs.columns
print(x.obs.head())
# Extraigo los datos de los counts, features, genotype, cb, tissue
df_gao = x.obs[["nCount_RNA", "nFeature_RNA", "genotype", "cb", "tissue"]]
df_gao.to_pickle("data-gao/df_gao.pkl")


# Cargo el df

df_gao = pd.read_pickle("data-gao/df_gao.pkl")

# Cantidad de neuronas
sum(df_gao.tissue == "Neuron")

# añadimos timepoint
days = [x for x in df_gao.genotype.astype(str).str.split("D")]
df_gao["timepoint"] = df_gao.genotype.astype(str).str.split("D").str[-1]

# Strains
# no se si debo de juntar las strains N2 y N2HT, -> NO, no debjo de juntar ya que distintos medios puede que
# distinta expresión 
# wildtype_mask = df_gao.genotype.astype(str).str.startswith("N2")
# df_gao["is_wt"] = wildtype_mask
# df_gao[df_gao.is_wt].genotype.unique()

# Aislo solamente las WT N2
df_gao_wt = df_gao[df_gao.genotype.str.startswith("N2D")]


# Quito la categoría 1early de timepoint
df_gao_wt = df_gao_wt[df_gao_wt.timepoint != "1early"]

# Sorting
df_gao_sorted = df_gao_wt
df_gao_sorted.timepoint = df_gao_wt.timepoint.astype(int)
df_gao_sorted = df_gao_wt.sort_values(by="timepoint")

# Saving
df_gao_sorted["is_neuron"] = (df_gao.tissue == "Neuron")
df_gao_sorted = df_gao_sorted.drop("tissue", axis= 1)
df_gao_sorted.columns = ['n_counts', 'n_genes', 'genotype', 'cb', 'timepoint', 'is_neuron']
df_gao_sorted.to_pickle("data-gao/df_gao_srt.pkl")

# Violin de counts y genes por celula por timepoint
mask = df_gao_sorted.is_wt
vl_ngenes = ex.violin(
    df_gao_sorted, 
    x="timepoint",
    y="nFeature_RNA",
    box=True,
    points="outliers",
    color="timepoint",
)
vl_ngenes.show()
vl_ngenes.write_html("viloin_ngenes_by_timepoint.html")

# COMPARACION OP Y HT

# Aislo N2, para diferenciar OP y HT
df_gao_compwt = df_gao[df_gao.genotype.str.startswith("N2")]

# Como solo hay OP y HT en d1, asilo d1
df_gao_compwt = df_gao_compwt[df_gao_compwt.timepoint == "1"]
df_gao_compwt["wt_strain"] = df_gao_compwt.genotype.str.contains("HT").map({True: "HT", False: "OP"})

# Parece que hay 3 n2 
n2_type = [x[0] for x in df_gao_compwt.index.str.split("D")]
set(n2_type)
df_gao_compwt["n2_type"] = n2_type

df_n2 = df_gao_compwt.groupby("n2_type")[["nCount_RNA", "nFeature_RNA"]].mean()
df_n2 = df_n2.reset_index()
df_n2.to_pickle("data-gao/df_n2.pkl")

# Barplot n2_type
fig_n2 = ex.bar(
    df_n2,
    x = "n2_type",
    y = "nCount_RNA",
    color = "n2_type",
)
fig_n2.show()

fig_n2 = go.Figure()
fig_n2.add_trace(go.Bar(
    x = df_n2["n2_type"],
    y = df_n2["nCount_RNA"],
    name="count",
    marker_color = "green",
    ))
fig_n2.add_trace(go.Bar(
    x = df_n2["n2_type"],
    y = df_n2["nFeature_RNA"],
    name="features",
    marker_color = "red",
    ))
