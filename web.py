from shiny import App, render, ui
from shinywidgets import output_widget, render_widget  
import plotly.express as ex
import plotly.graph_objects as go
import pandas as pd

# data-roux load
df_counts_r = pd.read_pickle("data-roux/df_counts.pkl")
df_timepoints_g_r = pd.read_pickle("data-roux/df_timepoints_g.pkl")
df_ng_nc_time = pd.read_pickle("data-roux/df_ng_nc_time_rtg.pkl")

# data-taylor load
df_counts_t = pd.read_pickle("data-taylor/df_counts.pkl")
# No hace falta porque está integrado en df_ng_nc_time_r
# df_ng_nc_time_t = pd.read_pickle("data-taylor/df_ng_nc.pkl")

# data-gao load
# No hace falta porque está integrado en df_ng_nc_time_rt
# df_counts_g = pd.read_pickle("data-gao/df_gao_srt.pkl")

app_ui = ui.page_fluid(
    ui.card(
        ui.h1("Dataset de Roux", style="text-align: center;"),
        ui.a("Página web con los datos", href="https://c.elegans.aging.atlas.research.calicolabs.com/"),
        ui.p("Dataset con datos single-cell de todo el organismo de C.elegans, cuentan con 6 timepoints (day 1, day 3, day 5, day 8, day 11 y day 15)"),
    ),
    ui.card(
        ui.card_header("Conteos celulares en dÍa 1"),
        ui.layout_sidebar(
            ui.sidebar(
                "Sidebar",
                ui.input_switch("specific", "Show only neurons", value=False),
                ui.input_slider("slider", "Slider", min=0, max=df_counts_r.shape[0], value=[0, df_counts_r.shape[0]]), 
            ),
            output_widget("roux_cell_counts_t1_barp"),
        ),
        ui.p(
            "Información obtenida agrupando por nombre las células observadas en el dia 1 de adultez",
        ),
        full_screen=True,
    ),
    ui.card(
        ui.card_header("Conteos celulares y de neuronas por día"),
        output_widget("roux_cell_counts_tx_barp"),
        ui.p(
            "Información obtenida agrupando por nombre las células observadas para cada día. " \
            "En gris se representan todas las células, en rojo solamente las neuronas",
        ),

        full_screen=True,
    ),
    ui.card(
        ui.input_switch("specific3", "Show only neurons", value=False),
        ui.input_selectize(  
            "var",  
            "Select an option below:",  
            {"n_genes": "Genes", "n_counts": "Counts"},  
        ), 
        output_widget("roux_ng_nc_tx_violin"),
        ui.p(
            "Información obtenida con plotly.express.violin.",
        ),

        ui.p(
            "Mediana de genes y conteos expresados por célula y timepoint:",
        ),
        ui.p(
            "day 1: genes: 529, conteos: 1259",
        ),
        ui.p(
            "day 3: genes: 582, conteos: 1413",
        ),
        ui.p(
            "day 5: genes: 575, conteos: 1358",
        ),
        ui.p(
            "day 8: genes: 568, conteos: 1411",
        ),
        ui.p(
            "day 11: genes: 434, conteos: 934",
        ),
        ui.p(
            "day 15: genes: 393, conteos: 891",
        ),

        full_screen=True,
    ),
    ui.card(
        ui.h5("Cantidad de categorías de neuronas detectadas por timepoint:"),
        ui.p("t1: 124"),
        ui.p("t3: 123"),
        ui.p("t5: 123"),
        ui.p("t8: 122"),
        ui.p("t11: 98"),
        ui.p("t15: 5"),
    ),

    ui.card(
        ui.h1("Dataset de Taylor - Adult hermaphrodite", style="text-align: center;"),
        ui.a("Enlace de descarga", href="https://cengen.org/storage/032224_L4_all_cells_Seurat5.rds"),
        ui.p("Dataset con datos single-cell específicos de neuronas en adultos hermafroditas"),
    ),
    ui.card(
        ui.h5("Conteos celulares"),
        ui.layout_sidebar(
            ui.sidebar(
                "Sidebar",
                ui.input_switch("specifict", "Show only neurons", value=False),
            ),
            output_widget("taylor_cell_counts"),
        ),
        ui.p("Información obtenida con el mismo procedimiento que en Roux. Cantidad de categorías de neuronas encontradas: 130")
    ),
    ui.card(
        ui.card_header("Comparación entre los datasets de Taylor y el day 1 de Roux"),
        ui.input_switch("specifict3", "Show only neurons", value=False),
        ui.input_selectize(  
            "vart",  
            "Select an option below:",  
            {"n_genes": "Genes", "n_counts": "Counts"},  
        ), 
        output_widget("taylor_ng_nc_violin"),
        ui.p("En azul (R) los datos de Roux, en rojo (T) los de taylor.")
    ),




    ui.card(
        ui.h5("Violin plot"),
        ui.input_switch("specificg3", "Show only neurons", value=False),
        ui.input_selectize(  
            "varg",  
            "Select an option below:",  
            {"n_genes": "Genes", "n_counts": "Counts"},  
        ), 
        output_widget("gao_ng_nc_tx_violin"),

        full_screen=True,
    ),
)


def server(input, output, session):
    @render_widget
    def roux_cell_counts_t1_barp():
        # switch neuron specific
        if input.specific():
            df1 = df_counts_r[df_counts_r["is_neuron"]]
        else:
            df1 = df_counts_r

        fig1 = ex.bar(
            df1, 
            x="annotate_name", 
            y="count"
            )
        fig1.update_traces(marker_color=df1["is_neuron"].map({True: "red", False: "gray"}))
        fig1.update_xaxes(range=[input.slider()[0], input.slider()[1]],
                          tickangle=-45,
                          tickfont=dict(size=10))
        return fig1


    @render_widget
    def roux_cell_counts_tx_barp():
        # switch neuron specific 
        fig2 = go.Figure()
        fig2.add_trace(go.Bar(
            x = df_timepoints_g_r["timepoint"],
            y = df_timepoints_g_r["count"],
            name = "cel",
            marker_color="lightgrey",
        ))
        fig2.add_trace(go.Bar(
            x = df_timepoints_g_r["timepoint"],
            y = df_timepoints_g_r["neuron_count"],
            name = "neurons",
            marker_color="indianred"
        ))
        fig2.update_layout(barmode="group")

        return fig2


    @render_widget
    def roux_ng_nc_tx_violin():
        if input.var() == "n_genes": 
            str_var = "genes"
        else: 
            str_var = "conteos"

        dfr = df_ng_nc_time[df_ng_nc_time["art"] == "R"]
        if input.specific3():
            df3 = dfr[dfr["is_neuron"]]
            title = "Número de " + str_var + " por neurona para cada timepoint"
        else:
            df3 = dfr
            title = "Número de " + str_var + " por célula para cada timepoint"

        fig3 = ex.violin(
            df3,
            x="timepoint", # no se porque commente esto
            y=input.var(),
            color="timepoint",
            box=True,
            points="outliers",
            title=title,
        )
        fig3.update_layout(showlegend=False)

        return fig3


    @render_widget
    def taylor_cell_counts():
        # switch neuron specific
        if input.specifict():
            df1 = df_counts_t[df_counts_t["is_neuron"]]
        else:
            df1 = df_counts_t

        fig1 = ex.bar(
            df1, 
            x="Cell.type", 
            y="count"
            )
        fig1.update_traces(marker_color=df1["is_neuron"].map({True: "red", False: "gray"}))
        fig1.update_xaxes(tickangle=-45,
                          tickfont=dict(size=10))

        return fig1

    @render_widget
    def taylor_ng_nc_violin():
        if input.vart() == "n_genes": 
            str_var = "genes"
        else: 
            str_var = "conteos"

        if input.specifict3():
            df3 = df_ng_nc_time[df_ng_nc_time["is_neuron"]]
            title = "Número de " + str_var + " por neurona "
        else:
            df3 = df_ng_nc_time
            title = "Número de " + str_var + " por célula"

        fig3 = ex.violin(
            df3,
            x="art",
            y=input.var(),
            box=True,
            points="outliers",
            color="art",
            title=title,
        )
        fig3.update_layout(showlegend=False)

        return fig3      


    @render_widget
    def gao_ng_nc_tx_violin():
        if input.varg() == "n_genes": 
            str_var = "genes"
        else: 
            str_var = "conteos"

        dfg = df_ng_nc_time[df_ng_nc_time["art"] == "G"]
        if input.specificg3():
            df3 = dfg[dfg["is_neuron"]]
            title = "Número de " + str_var + " por neurona para cada timepoint"
        else:
            df3 = dfg
            title = "Número de " + str_var + " por célula para cada timepoint"

        fig3 = ex.violin(
            df3,
            x="timepoint",
            y=input.var(),
            color="timepoint",
            box=True,
            points="outliers",
            title=title,
        )
        fig3.update_layout(showlegend=False)

        return fig3
     


app = App(app_ui, server)