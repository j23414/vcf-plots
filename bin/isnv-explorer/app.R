library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(Biostrings)

# ============================================================
# DATA PROCESSING
# ============================================================

chromosome_order = c("PV062510|PB2", "PV074323|PB1", "PV062508|PA", "PV062513|HA",
                       "PV062509|NP", "PV062511|NA", "PV062507|MP", "PV062512|NS")

expand_data <- read_delim("/Users/jchang99/github/j23414/vcf-plots/data/expand_data_100000_0100.tsv", delim="\t", na=character()) %>%
  mutate(
    CHROM=factor(CHROM, levels=chromosome_order)
  )

depth_data <- read_delim("/Users/jchang99/github/j23414/vcf-plots/data/combined_depth.tsv", delim = "\t", na = character()) %>%
  mutate(
    CHROM = factor(CHROM, levels = chromosome_order),
    Depth = case_when(Depth < 1 ~ 1, TRUE ~ Depth)
  )

filtered_depth_data <- read_delim("/Users/jchang99/github/j23414/vcf-plots/data/combined_depth2.tsv", delim = "\t", na = character()) %>%
  mutate(
    CHROM = factor(CHROM, levels = chromosome_order),
    Depth = case_when(Depth < 1 ~ 1, TRUE ~ Depth)
)

mindepth=1000

# ============================================================
# Reusable plotting functions
# ============================================================

#' Plot variants and sequencing depth for a sample
#'
#' Creates a ggplot showing total sequencing depth, filtered sequencing
#' depth, and minor variants across one or more genome segments.
#'
#' @param samplename Character string containing the sample name.
#' @param depth_total Data frame containing sequencing depth before filtering.
#'   Must contain `POS` and `Depth` columns.
#' @param depth_filtered Data frame containing sequencing depth after filtering
#'   (e.g., primer trimming, duplicate removal, and base-quality filtering).
#'   Must contain `POS` and `Depth` columns.
#' @param sample_variants Data frame containing major and minor variants for the sample.
#'   Must contain the columns `POS`, `minor_percentage`, `IsSynonymous`,
#'   `mutation`, `Gene`, `minorCodon`, and `minorAminoAcid`.
#' @param mindepth Minimum depth threshold to display as a horizontal line.
#'
#' @return A ggplot object showing sequencing depth and minor variants. Can be passed to ggplotly
#'
#' @examples
#' plot_variants_and_depth(
#'   samplename = "sample_001",
#'   depth_total = depth_total,
#'   depth_filtered = depth_filtered,
#'   sample_variants = sample_variants, # for sample_001
#'   mindepth = 1000
#' )

plot_variants_and_depth <- function(samplename, depth_total, depth_filtered, sample_variants, mindepth = 1000){
  p <- ggplot() +
    # Depth
    geom_line(data = depth_filtered,
              aes(x = POS, y = Depth),
              color = "#D5DDE0") +
    geom_line(data = depth_total,
              aes(x = POS, y = Depth),
              color = "#71838C") +
    # iSNVs
    geom_point(
      data = sample_variants,
      aes(
        x = POS,
        # Put iSNVs in the lower portion of the depth plot
        y = (
          max(log10(depth_total$Depth), na.rm = TRUE) / 4 +
            minor_percentage *
            max(log10(depth_total$Depth), na.rm = TRUE) / 2
        )^8,
        color = IsSynonymous,
        text = paste0(
          "Position: ", POS,
          "<br>Mutation: ", mutation,
          "<br>Frequency: ", scales::percent(minor_percentage),
          "<br>Gene: ", Gene,
          "<br>Codon Position: ", CodonPosition,
          "<br>Codon: ", minorCodon,
          "<br>Amino acid: ", minorAminoAcid
        )
      ),
      size = 2
    ) +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 2) +
    geom_hline(yintercept = mindepth,
               linetype = "dashed",
               color = "red") +
    scale_y_log10() +
    labs(
      title = paste("Depth Plot for Sample:", sample),
      x = "Position",
      y = "Depth",
      color = "Mutation type"
    ) +
    theme_bw()

  p
}

# ============================================================
# Generate static visuals
# ============================================================
# samples = unique(expand_data$Sample)
#
# for (sample in samples ){
#   sample_variants <- expand_data %>%
#     filter(Sample %in% sample)
#
#   sample_depth <- depth_data %>%
#     filter(Sample == sample)
#
#   filtered_sample_depth <- filtered_depth_data %>%
#     filter(Sample == sample)
#
#   p <- plot_variants_and_depth(
#     sample,
#     sample_depth,
#     filtered_sample_depth,
#     sample_variants,
#     mindepth)
#
#   ggsave(filename=paste0("/Users/jchang99/Desktop/iSNV/Results_100000maxreads_0100mindepth/samples/",sample,".png", sep=""), plot=p, height=8, width=8)
# }


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

      sliderInput(
        "min_freq",
        "Minimum iSNV frequency",
        min = 0,
        max = 1.0,
        value = 0,
        step = 0.01
      )
    ),

    mainPanel(
      h3("iSNVs across genome"),
      plotlyOutput("genome_plot", height = "800px"),
      h3("Mutational spectrum"),
      plotlyOutput("spectrum_plot", height = "410px"),
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
             minor_percentage >= input$min_freq
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
            "<br>Codon Position: ", CodonPosition,
            "<br>Codon: ", minorCodon,
            "<br>Amino acid: ", minorAminoAcid
          )
        )) +
        geom_point(size = 1) +
        facet_wrap(~ CHROM, ncol = 1) +
        theme_bw() +
        labs(x = "Genomic position", y = "iSNV frequency", color = "Mutation type")

      return(ggplotly(p, tooltip = "text", width=800, height=800))
    }


    # ----------------------------------------------------------
    # SINGLE SAMPLE
    # ----------------------------------------------------------

    sample_name <- single_sample()

    sample_depth <- depth_data %>%
      filter(Sample == sample_name)

    filtered_sample_depth <- filtered_depth_data %>%
      filter(Sample == sample_name)

    # Use only variants from this sample
    sample_variants <- variants %>%
      filter(Sample == sample_name)

    p <- plot_variants_and_depth(
      sample_name,
      sample_depth,
      filtered_sample_depth,
      sample_variants,
      mindepth = mindepth
    )

    ggplotly(p, tooltip = "text", width=800, height=600)
  })

  # ----------------------------------------------------------
  # MUTATIONAL SPECTRUM
  # ----------------------------------------------------------

  output$spectrum_plot <- renderPlotly({
    sdata <- sample_data() %>%
      subset(consensus != variant) %>%
      count(CHROM, mutation, IsSynonymous)

    p <- ggplot(sdata, aes(x = mutation, y = n, fill = IsSynonymous)) +
      geom_col() +
      facet_wrap(~ CHROM, ncol = 4,
                 scales = "free_y"
                 ) +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 7
      )) +
      labs(
        x = "Mutation",
        y = "Number of iSNVs",
        fill = "Mutation type")

    ggplotly(p, width=800, height=400)
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
