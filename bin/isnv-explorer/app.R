library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(readxl)
library(Biostrings)

# ============================================================
# DATA PROCESSING
# ============================================================

expand_data <- read_delim("/Users/jchang99/github/j23414/vcf-plots/data/expand_data.tsv", delim="\t", na=character()) %>%
  mutate(
    CHROM=factor(CHROM, levels=c("PV062510|PB2", "PV074323|PB1", "PV062508|PA", "PV062513|HA",
                                 "PV062507|MP", "PV062511|NA", "PV062509|NP", "PV062512|NS"))
  )

depth_data <- read_delim("/Users/jchang99/github/j23414/vcf-plots/data/combined_depth.tsv", delim = "\t", na = character()) %>%
  mutate(CHROM = factor(CHROM,
                        levels = c("PV062510|PB2", "PV074323|PB1", "PV062508|PA","PV062513|HA",
                                   "PV062507|MP", "PV062511|NA", "PV062509|NP", "PV062512|NS")))

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  titlePanel("iSNV Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "sample",
        "Sample",
        choices = c("All", sort(unique(expand_data$Sample))),
        selected = "All"#,
        #multiple = TRUE
      ),

      selectInput(
        "segment",
        "Segment",
        choices = c("All", levels(expand_data$CHROM)),
        selected = "All",
        multiple = TRUE
      ),

      checkboxGroupInput(
        "mutation_type",
        "Mutation type",
        choices = c("Synonymous", "Non-Synonymous", "Intergenic"),
        selected = c("Synonymous", "Non-Synonymous", "Intergenic")
      ),

      # sliderInput(
      #   "min_freq",
      #   "Minimum iSNV frequency",
      #   min = 0,
      #   max = 0.5,
      #   value = 0,
      #   step = 0.01
      # )
    ),

    mainPanel(
      h3("iSNVs across genome"),
      plotlyOutput("genome_plot", height = "800px"),
      h3("Mutational spectrum"),
      plotlyOutput("spectrum_plot"),
      h3("Variants"),
      DTOutput("variant_table")
    )
  ))


# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  # ----------------------------------------------------------
  # FILTER DATA
  # ----------------------------------------------------------

  sample_data <- reactive({
    selected_samples <- input$sample
    selected_segments <- input$segment

    # Remove "All" if other samples are selected
    if (length(selected_samples) > 1 &&
        "All" %in% selected_samples) {
      selected_samples <- selected_samples[selected_samples != "All"]
      updateSelectInput(session, "sample", selected = selected_samples)
    }

    # Remove "All" if other segments are selected
    if (length(selected_segments) > 1 &&
        "All" %in% selected_segments) {
      selected_segments <- selected_segments[selected_segments != "All"]
      updateSelectInput(session, "segment", selected = selected_segments)
    }

    x <- expand_data %>%
      filter(IsSynonymous %in% input$mutation_type,
             # minor_percentage >= input$min_freq
             )

    if (!"All" %in% c(input$sample)) {
      x <- x %>%
        filter(Sample %in% c(input$sample))
    }

    if (!"All" %in% c(input$segment)) {
      x <- x %>%
        filter(CHROM == input$segment)
    }

    x
  })

  single_sample <- reactive({
    selected <- input$sample

    if ("All" %in% selected) {
      return(NULL)
    }

    if (length(selected) == 1) {
      return(selected)
    }

    NULL
  })


  # ----------------------------------------------------------
  # GENOME PLOT
  # ----------------------------------------------------------

  output$genome_plot <- renderPlotly({
    variants <- sample_data()

    # ----------------------------------------------------------
    # MULTIPLE SAMPLES
    # ----------------------------------------------------------

    if (is.null(single_sample())) {
      p <- variants %>%
        ggplot(aes(
          x = POS,
          y = minor_percentage,
          color = IsSynonymous,
          text = paste0(
            "Sample: ", Sample,
            "<br>Position: ", POS,
            "<br>Mutation: ", mutation,
            "<br>Frequency: ",
            scales::percent(minor_percentage),
            "<br>Gene: ", Gene,
            "<br>Codon: ", minorCodon,
            "<br>Amino acid: ", minorAminoAcid
          )
        )) +
        geom_point(size = 1) +
        facet_wrap(~ CHROM, ncol = 1) +
        theme_bw() +
        labs(x = "Genomic position", y = "iSNV frequency", color = "Mutation type")

      return(ggplotly(p, tooltip = "text"))
    }


    # ----------------------------------------------------------
    # SINGLE SAMPLE
    # ----------------------------------------------------------

    sample_name <- single_sample()

    sample_depth <- depth_data %>%
      filter(Sample == sample_name)

    # Use only variants from this sample
    sample_variants <- variants %>%
      filter(Sample == sample_name)

    p <- ggplot() +
      # Depth
      geom_line(data = sample_depth,
                aes(x = POS, y = Depth),
                color = "blue") +
      # iSNVs
      geom_point(
        data = sample_variants,
        aes(
          x = POS,
          # Put iSNVs in the lower portion of the depth plot
          y = (
            max(log10(sample_depth$Depth), na.rm = TRUE) / 4 +
              minor_percentage *
              max(log10(sample_depth$Depth), na.rm = TRUE) / 2
          )^8,
          color = IsSynonymous,
          text = paste0(
            "Position: ", POS,
            "<br>Mutation: ", mutation,
            "<br>Frequency: ",
            scales::percent(minor_percentage),
            "<br>Gene: ", Gene,
            "<br>Codon: ", minorCodon,
            "<br>Amino acid: ", minorAminoAcid
          )
        ),
        size = 2
      ) +
      facet_wrap(~ CHROM, scales = "free_x", ncol = 2) +
      geom_hline(yintercept = 1000,
                 linetype = "dashed",
                 color = "red") +
      scale_y_log10() +
      labs(
        title = paste("Depth Plot for Sample:", sample_name),
        x = "Position",
        y = "Depth",
        color = "Mutation type"
      ) +
      theme_bw()

    ggplotly(p, tooltip = "text")
  })

  # ----------------------------------------------------------
  # MUTATIONAL SPECTRUM
  # ----------------------------------------------------------

  output$spectrum_plot <- renderPlotly({
    sdata <- sample_data() %>%
      count(CHROM, mutation, IsSynonymous)

    p <- ggplot(sdata, aes(x = mutation, y = n, fill = IsSynonymous)) +
      geom_col() +
      facet_wrap(~ CHROM, ncol = 4, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 7
      )) +
      labs(x = "Mutation", y = "Number of iSNVs", fill = "Mutation type")

    ggplotly(p)
  })


  # ----------------------------------------------------------
  # VARIANT TABLE
  # ----------------------------------------------------------

  output$variant_table <- renderDT({
    sample_data() %>%
      select(
        Sample,
        CHROM,
        POS,
        consensus,
        variant,
        mutation,
        minor_percentage,
        Gene,
        SNPCodonPosition,
        majorCodon,
        minorCodon,
        majorAminoAcid,
        minorAminoAcid,
        IsSynonymous
      ) %>%
      arrange(CHROM, POS)
  })
}

# ============================================================
# RUN APP
# ============================================================

shinyApp(ui = ui, server = server)