from shiny import App, render, ui
from shinywidgets import output_widget, render_widget  
import plotly.express as ex
import plotly.graph_objects as go
import pandas as pd

# data-roux load
df_counts_r = pd.read_pickle("data-roux/df_counts.pkl")
df_timepoints_g_r = pd.read_pickle("data-roux/df_timepoints_g.pkl")
df_ng_nc_time = pd.read_pickle("data-roux/df_ng_nc_time_rt.pkl")

# data-taylor load
df_counts_t = pd.read_pickle("data-taylor/df_counts.pkl")
# df_ng_nc_time_t = pd.read_pickle("data-taylor/df_ng_nc.pkl") # No hace falta porque la info está en df_ng_nc_time_rt


app_ui = ui.page_fluid(
    ui.card(
        ui.card_header("Cell counts (d1)"),
        ui.layout_sidebar(
            ui.sidebar(
                "Sidebar",
                ui.input_switch("specific", "Show only neurons", value=False),
                ui.input_slider("slider", "Slider", min=0, max=df_counts_r.shape[0], value=[0, df_counts_r.shape[0]]), 

            ),
            output_widget("roux_cell_counts_t1_barp"),
        ),

        full_screen=True,
    ),
    ui.card(
        ui.card_header("Cell counts by timepoint"),
        output_widget("roux_cell_counts_tx_barp"),

        full_screen=True,
    ),
    ui.card(
        ui.h5("Violin plot"),
        ui.input_switch("specific3", "Show only neurons", value=False),
        ui.input_selectize(  
            "var",  
            "Select an option below:",  
            {"n_genes": "Genes", "n_counts": "Counts"},  
        ), 
        output_widget("roux_ng_nc_tx_violin"),

        full_screen=True,
    ),
    ui.card(
        ui.h5("Tipo de neurona por timepoint:"),
        ui.h5("t1: 124"),
        ui.h5("t3: 123"),
        ui.h5("t5: 123"),
        ui.h5("t8: 122"),
        ui.h5("t11: 98"),
        ui.h5("t15: 5"),
    ),
    ui.card(
        ui.h5("Celullar counts comparation between Taylor and Roux"),
        output_widget("taylor_cell_counts"),
    ),
    ui.card(
        ui.input_switch("specifict3", "Show only neurons", value=False),
        ui.input_selectize(  
            "vart",  
            "Select an option below:",  
            {"n_genes": "Genes", "n_counts": "Counts"},  
        ), 
        ui.h5("Violin plot taylor"),
        output_widget("taylor_ng_nc_violin"),
    ),
)


def server(input, output, session):
    @render_widget
    def roux_cell_counts_t1_barp():
        # switch neuron specific
        if input.specific():
            df1 = df_counts_r[df_counts_r["is.neuron"]]
        else:
            df1 = df_counts_r

        fig1 = ex.bar(
            df1, 
            x="annotate_name", 
            y="count"
            )
        fig1.update_traces(marker_color=df1["is.neuron"].map({True: "red", False: "gray"}))
        fig1.update_xaxes(range=[input.slider()[0], input.slider()[1]])


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
            df3 = dfr[dfr["is.neuron"]]
            title = "Número de " + str_var + " por neurona para cada timepoint"
        else:
            df3 = dfr
            title = "Número de " + str_var + " por célula para cada timepoint"

        fig3 = ex.violin(
            df3,
            # x="timepoint",
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
        if input.specific():
            df1 = df_counts_t[df_counts_t["is.neuron"]]
        else:
            df1 = df_counts_t

        fig1 = ex.bar(
            df1, 
            x="Cell.type", 
            y="count"
            )
        fig1.update_traces(marker_color=df1["is.neuron"].map({True: "red", False: "gray"}))
        fig1.update_xaxes(range=[input.slider()[0], input.slider()[1]])

        return fig1

    @render_widget
    def taylor_ng_nc_violin():
        if input.vart() == "n_genes": 
            str_var = "genes"
        else: 
            str_var = "conteos"

        if input.specifict3():
            df3 = df_ng_nc_time[df_ng_nc_time["is.neuron"]]
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


app = App(app_ui, server)