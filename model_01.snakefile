rule plot_fcs:
    input:
        fcs = "data/" + PROJECT + "_{sample}_{replicate}.fcs",
        labels = "data/" + PROJECT + "_{sample}_{replicate}_labels.csv"
    output:
        png_56_3 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_CD56_CD3.png",
        png_56_19 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_CD56_CD19.png",
        png_8_4 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_CD8_CD4.png",
        png_56 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_density_CD56.png",
        png_3 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_density_CD3.png",
        png_19 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_density_CD19.png",
        png_8 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_density_CD8.png",
        png_4 = "exprs_plots/" + PROJECT + "_{sample}_{replicate}_density_CD4.png"
    log:
        "logs/plot_fcs/" + PROJECT + "_{sample}_{replicate}.log"
    script:
        "scripts/plot_fcs_model_01.R"
