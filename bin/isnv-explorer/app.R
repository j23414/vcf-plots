library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(readxl)

# ============================================================
# DATA PROCESSING
# ============================================================

data <- read_excel(
  "/Users/jchang99/github/j23414/vcf-plots/bin/data/samples_variant_results_longer.xlsx"
)

# codon_to_aa <- ...

complete_data <- data %>%
  select(
    -IsSynonymous,
    -AltAminoAcid,
    -AltCodon,
    -IsTransition,
    -AminoAcidChange
  ) %>%
  filter(
    VariantType == "SNP",
    Gene != "PB1-F2"
  ) %>%
  mutate(
    Percentage = as.numeric(Percentage),
    CHROM = factor(
      CHROM,
      levels = c(
        "PV062510_PB2",
        "PV074323_PB1",
        "PV062508_PA",
        "PV062513_HA",
        "PV062507_MP",
        "PV062511_NA",
        "PV062509_NP",
        "PV062512_NS"
      )
    )
  ) %>%
  pivot_wider(
    names_from = ALT,
    values_from = Percentage
  ) %>%
  rowwise() %>%
  mutate(
    total = sum(C, A, G, T, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    A = case_when(
      REF == "A" ~ 1 - total,
      TRUE ~ A
    ),
    C = case_when(
      REF == "C" ~ 1 - total,
      TRUE ~ C
    ),
    G = case_when(
      REF == "G" ~ 1 - total,
      TRUE ~ G
    ),
    T = case_when(
      REF == "T" ~ 1 - total,
      TRUE ~ T
    )
  ) %>%
  rowwise() %>%
  mutate(
    major_percentage = max(A, C, G, T, na.rm = TRUE),
    
    consensus = case_when(
      A == major_percentage ~ "A",
      G == major_percentage ~ "G",
      C == major_percentage ~ "C",
      T == major_percentage ~ "T"
    ),
    
    SNPCodonPosition = as.numeric(SNPCodonPosition),
    
    majorCodon = case_when(
      IsGenic ~ paste0(
        substr(RefCodon, 1, SNPCodonPosition),
        consensus,
        substr(RefCodon, SNPCodonPosition + 2, 3)
      ),
      TRUE ~ NA
    ),
    
    majorAminoAcid = codon_to_aa[majorCodon]
  ) %>%
  ungroup()


expand_data <- complete_data %>%
  pivot_longer(
    cols = c(A, C, G, T),
    names_to = "variant",
    values_to = "minor_percentage"
  ) %>%
  filter(
    !is.na(minor_percentage),
    consensus != variant,
    minor_percentage != 0
  ) %>%
  mutate(
    mutation = paste0(consensus, "->", variant),
    
    minorCodon = case_when(
      IsGenic ~ paste0(
        substr(majorCodon, 1, SNPCodonPosition),
        variant,
        substr(majorCodon, SNPCodonPosition + 2, 3)
      ),
      TRUE ~ NA
    ),
    
    minorAminoAcid = codon_to_aa[minorCodon],
    
    IsSynonymous = case_when(
      is.na(majorAminoAcid) ~ "Intergenic",
      majorAminoAcid == minorAminoAcid ~ "Synonymous",
      TRUE ~ "Non-Synonymous"
    )
  )


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
        choices = c("All",sort(unique(expand_data$Sample))),
        selected="All"
      ),
      
      selectInput(
        "segment",
        "Segment",
        choices = c("All", levels(expand_data$CHROM))
      ),
      
      checkboxGroupInput(
        "mutation_type",
        "Mutation type",
        choices = c(
          "Synonymous",
          "Non-Synonymous",
          "Intergenic"
        ),
        selected = c(
          "Synonymous",
          "Non-Synonymous",
          "Intergenic"
        )
      ),
      
      sliderInput(
        "min_freq",
        "Minimum iSNV frequency",
        min = 0,
        max = 0.5,
        value = 0,
        step = 0.01
      )
    ),
    
    mainPanel(
      
      h3("iSNVs across genome"),
      
      plotlyOutput(
        "genome_plot",
        height = "800px"
      ),
      
      h3("Mutational spectrum"),
      
      plotlyOutput(
        "spectrum_plot"
      ),
      
      h3("Variants"),
      
      DTOutput(
        "variant_table"
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  # ----------------------------------------------------------
  # FILTER DATA
  # ----------------------------------------------------------
  
  sample_data <- reactive({
    
    x <- expand_data %>%
      filter(
        IsSynonymous %in% input$mutation_type,
        minor_percentage >= input$min_freq
      )
    
    if (input$sample != "All") {
      x <- x %>%
        filter(Sample == input$sample)
    }
    
    if (input$segment != "All") {
      x <- x %>%
        filter(CHROM == input$segment)
    }
    
    x
  })
  
  
  # ----------------------------------------------------------
  # GENOME PLOT
  # ----------------------------------------------------------
  
  output$genome_plot <- renderPlotly({
    
    p <- sample_data() %>%
      ggplot(
        aes(
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
        )
      ) +
      geom_point(size = 1) +
      facet_wrap(
        ~CHROM,
        ncol = 1
      ) +
      theme_bw() +
      labs(
        x = "Genomic position",
        y = "iSNV frequency",
        color = "Mutation type"
      )
    
    ggplotly(
      p,
      tooltip = "text"
    )
  })
  
  
  # ----------------------------------------------------------
  # MUTATIONAL SPECTRUM
  # ----------------------------------------------------------
  
  output$spectrum_plot <- renderPlotly({
    
    sdata <- sample_data() %>%
      count(
        CHROM,
        mutation,
        IsSynonymous
      )
    
    p <- ggplot(
      sdata,
      aes(
        x = mutation,
        y = n,
        fill = IsSynonymous
      )
    ) +
      geom_col() +
      facet_wrap(
        ~CHROM,
        ncol = 4
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          size = 7
        )
      ) +
      labs(
        x = "Mutation",
        y = "Number of iSNVs",
        fill = "Mutation type"
      )
    
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

shinyApp(
  ui = ui,
  server = server
)