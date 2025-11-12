from shiny import App, render, ui
from shinywidgets import output_widget, render_widget  
import plotly.express as ex
import plotly.graph_objects as go
import pandas as pd

# data load
df_counts = pd.read_pickle("data/df_counts.pkl")
df_timepoints_g = pd.read_pickle("data/df_timepoints_g.pkl")
df_ng_nc_time = pd.read_pickle("data/df_ng_nc_time.pkl")


app_ui = ui.page_fluid(
    ui.card(
        ui.card_header("Cell counts (d1)"),
        ui.layout_sidebar(
            ui.sidebar(
                "Sidebar",
                ui.input_switch("specific", "Show only neurons", value=False),
                ui.input_slider("slider", "Slider", min=0, max=df_counts.shape[0], value=[0, df_counts.shape[0]]), 

            ),
            output_widget("bar_counts_t1"),
        ),

        full_screen=True,
    ),
    ui.card(
        ui.card_header("Cell counts by timepoint"),
        output_widget("bar_counts"),

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
        output_widget("violin_ng_nc_t"),

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
)


def server(input, output, session):
    @render_widget
    def bar_counts_t1():
        # switch neuron specific
        if input.specific():
            df1 = df_counts[df_counts["is.neuron"]]
        else:
            df1 = df_counts

        fig1 = ex.bar(
            df1, 
            x="annotate_name", 
            y="count"
            )
        fig1.update_traces(marker_color=df1["is.neuron"].map({True: "red", False: "gray"}))
        fig1.update_xaxes(range=[input.slider()[0], input.slider()[1]])


        return fig1

    @render_widget
    def bar_counts():
        # switch neuron specific 

        fig2 = go.Figure()
        fig2.add_trace(go.Bar(
            x = df_timepoints_g["timepoint"],
            y = df_timepoints_g["count"],
            name = "cel",
            marker_color="lightgrey",
        ))
        fig2.add_trace(go.Bar(
            x = df_timepoints_g["timepoint"],
            y = df_timepoints_g["neuron_count"],
            name = "neurons",
            marker_color="indianred"
        ))
        fig2.update_layout(barmode="group")

        return fig2

    @render_widget
    def violin_ng_nc_t():
        if input.var() == "n_genes": 
            str_var = "genes"
        else: 
            str_var = "conteos"

        if input.specific3():
            df3 = df_ng_nc_time[df_ng_nc_time["is.neuron"]]
            title = "Número de " + str_var + " por neurona para cada timepoint"
        else:
            df3 = df_ng_nc_time
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