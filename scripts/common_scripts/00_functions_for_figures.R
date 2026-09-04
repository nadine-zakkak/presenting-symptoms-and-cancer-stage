# Helper functions to produce final figures 

# Stage by symptom visuals ------
plot_stage_sx <- function(df, df_sx_order, filter_by_sex, colours) {
  # only keep relevant sex
  df <- df |> filter(sex_lb == filter_by_sex)
  # only visualise if N>=100
  df <- df |> filter(N_tumour >= 100) 
  # Add labels to be displayed
  df <- df |> mutate(label = str_glue("N={N_tumour}"))
  
  # get symptom order and map to shorter symptom names
  order       <- df_sx_order[[filter_by_sex]]
  order_short <- setNames(sx_fig_rename$New, sx_fig_rename$Original)
  final       <- order_short[order]
  final       <- c("N/A", "N/K", setdiff(final, c("N/A", "N/K")))
  
  # only keep relevant symptoms
  df <- df |> filter(symptom %in% order)
  
  df |>
    mutate(symptom_short = factor(symptom_short, levels = final)) |>
    arrange(symptom_short) |>
    ggplot() +
    geom_bar(aes(y    = prop100, 
                 x    = symptom_short, 
                 fill = reorder(stage_lb, desc(stage_lb))),
             stat  = "identity", 
             width = .5) +
    ylim(-12, 101) +
    geom_text(aes(x = symptom_short, y = -12, label = label),
              size   = 3.5, 
              hjust  = 0, 
              colour = "black") +
    guides(fill = guide_legend(reverse = T, title = "Stage")) +
    scale_fill_manual(values = rev(unname(unlist(colours)))) +
    coord_flip() +
    facet_grid(~sex_lb) +
    theme_minimal() +
    theme(text               = element_text(size = 14, colour = "black"),
          axis.text.y        = element_text(size = 14, colour = "black"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    xlab("Symptom") +
    ylab("Proportion (%)")
}

# Stage distribution by symptom, stratified by cancer site -----
plot_stage_sx_site <- function(df, y_indent = 1.5, colours) {
  
  df <- df |> mutate(label   = str_glue("N={N_tumour}"),
                     y_label = y_indent)
  
  cancer_sites <- names(cancer_redflag)
  
  # Create dataframe to be used for blank "canvas" of patchwork plot
  nonspecific <- 
    df |>
    filter(symptom %in% nonspecific_symptoms) |>
    bind_rows(data.frame = expand.grid(symptom = "Organ-specific alarm", 
                                       cancer_site_desc = cancer_sites, 
                                       prop100 = 2)) |>
    mutate(symptom = factor(symptom)) |>
    mutate(symptom = factor(symptom, levels = c("Abdominal pain (NOS)",
                                                "Nausea and/or vomiting",
                                                "Weight loss",
                                                "Organ-specific alarm"))) |>
    arrange(symptom)
  
  # consistent plot margin to be used across all plots
  add_plot_margin <-  theme(plot.margin  = margin(t = .1, 
                                                  r = .3, 
                                                  b = .1, 
                                                  l = .3, 
                                                  unit = "cm"))
  
  # alpha values to be used for transparency of bar plots
  alpha_vals <- c(1, rep(0.3, length(colours) - 1))
  
  # function to plot cancer site and all corresponding red flag symptoms
  plot_redflag <- function(data, cancer_site) {
    
    # Filter data to specific cancer site and order symptoms
    df_cancer <- 
      data |>
      filter(cancer_site_desc == cancer_site,
             symptom %in% c(cancer_redflag[[cancer_site]])) |>
      mutate(symptom = factor(symptom)) |>
      mutate(symptom = factor(symptom, levels = c(cancer_redflag[[cancer_site]]))) |>
      mutate(symptom_numeric = as.numeric(symptom)) |>
      arrange(symptom) 
    
    # Conditional bar width so it appears consistent in final combined figure
    barplt_width <- ifelse(length(unique(df_cancer$symptom)) == 1, .3, .5)
    
    df_cancer <-
      df_cancer |>
      rowwise() |>
      mutate(symptom_numeric = symptom_numeric + barplt_width/4) |>
      ungroup()
    
    # Bar plot with lines at 0 and 100
    p <- 
      ggplot() +
      geom_hline(yintercept = 0,   colour = "#3D3D3D", linewidth = .2) +
      geom_hline(yintercept = 100, colour = "#3D3D3D", linewidth = .2) +
      geom_bar(aes(y     = prop100, 
                   x     = symptom, 
                   fill  = reorder(stage_bin, desc(stage_bin)),
                   alpha = reorder(stage_bin, desc(stage_bin))),
               stat  = "identity", 
               width = barplt_width, 
               data  = df_cancer)
    
    # Add error bars to proportion stage I-III
    p <- 
      p +
      geom_errorbar(aes(ymin = prop100_lb,
                        ymax = prop100_ub,
                        x    = symptom),
                    width = barplt_width/4,
                    data = df_cancer |> filter(stage_bin == "1-3")) 
    
    # Add labels within bars and above bars. 
    # filter data to only include once since it is unique between the 2
    p <-
      p + 
      geom_text(size   = 9/.pt, 
                hjust  = 0, 
                colour = "white",
                aes(x  = symptom, y = y_label, label = label),
                data   = df_cancer) +
      geom_text(size = 10/.pt,
                aes(label = symptom, 
                    x     = symptom_numeric, 
                    y     = 1),
                hjust = 0,
                vjust = -1.7,
                data = df_cancer)
    
    
    # Adjust scales, colours,  flip coordinates and add facets
    p <-
      p +
      ylim(0, 102) +
      scale_fill_manual(values = rev(unname(unlist(colours)))) +
      scale_alpha_manual(NULL, values = rev(alpha_vals)) +
      scale_x_discrete(drop = FALSE) +
      coord_flip() +
      facet_grid2(. ~ cancer_site_desc)
    
    # Adjust theme
    p <-
      p +
      theme_minimal() +
      theme(axis.text          = element_blank(),
            axis.ticks         = element_blank(),
            axis.title         = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            strip.text         = element_text(face = "bold", size = 15, 
                                              margin = margin(t = .25, 
                                                              r = 0, 
                                                              b = .5, 
                                                              l = 0, 
                                                              unit = "cm")),
            legend.position    = "none") +
      add_plot_margin
    
    return(p)
  }
  
  # function to plot cancer site and all corresponding red flag symptoms
  plot_nonspecific <- function(data, cancer_site) {
    
    # Filter data to specific cancer site and order symptoms
    df_cancer <-
      data |>
      filter(cancer_site_desc == cancer_site,
             symptom %in% c(nonspecific_symptoms)) |>
      mutate(symptom = factor(symptom)) |>
      mutate(symptom = factor(symptom, levels = c(nonspecific_symptoms))) |>
      arrange(symptom) 
    
    # Set bar plot width
    barplt_width <- .4
    
    # Bar plot with lines at 0 and 100
    p <- 
      ggplot() +
      geom_hline(yintercept = 0,   colour = "#3D3D3D", linewidth = barplt_width/2) +
      geom_hline(yintercept = 100, colour = "#3D3D3D", linewidth = barplt_width/2) +
      geom_bar(aes(y     = prop100, 
                   x     = symptom, 
                   fill  = reorder(stage_bin, desc(stage_bin)),
                   alpha = reorder(stage_bin, desc(stage_bin))),
               stat  = "identity", 
               width = barplt_width, 
               data  = df_cancer) 
    
    # Add error bars for proportion of stages I-III
    p <- 
      p +
      geom_errorbar(aes(ymin = prop100_lb,
                        ymax = prop100_ub,
                        x    = symptom),
                    width = barplt_width/4,
                    data = df_cancer |> filter(stage_bin == "1-3"))
    
    # Add labels within bars. filter data to only include once since it is unique between the 2 stages
    p <- 
      p +
      geom_text(size   = 9/.pt, 
                hjust  = 0, 
                colour = "white",
                aes(x  = symptom, y = y_label, label = label),
                data   = df_cancer |> filter(stage_bin == "1-3"))
    
    # Adjust scales, colours and  flip coordinates
    p <- 
      p +
      ylim(0, 101) +
      scale_fill_manual(values = rev(unname(unlist(colours)))) +
      scale_alpha_manual(NULL, values = rev(alpha_vals)) +
      scale_x_discrete(drop = FALSE) +
      coord_flip()
    
    # Adjust theme
    p <- 
      p + 
      theme_minimal()+
      theme(axis.text.y        = element_blank(),
            axis.title         = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            legend.position    = "none") +
      add_plot_margin
    
    return(p)
  }
  
  # Plot blank "canvas" of patchwork plot
  plot_blank <- function(data) {
    
    # Building empty plot with Red flag label
    p_general_rf <- 
      ggplot() +
      geom_bar(aes(y = prop100, 
                   x = symptom, 
                   fill = reorder(stage_bin, desc(stage_bin))),
               stat = "identity", 
               width = .5, 
               alpha = 0,
               data = data |> filter(symptom == "Organ-specific alarm"))
    
    # Adjust scales, colours,  flip coordinates and add facets
    p_general_rf <-
      p_general_rf +
      ylim(0, 102) +
      scale_fill_manual(values = rep("white", 2)) +
      scale_x_reordered() +
      coord_flip() +
      facet_grid2(. ~ cancer_site_desc)
    
    # Adjust theme
    p_general_rf <-
      p_general_rf +
      theme_minimal() +
      theme(axis.text.y        = element_text(size = 14, colour = "black", 
                                              margin = margin(t = 0, 
                                                              r = 1.7, 
                                                              b = 0, 
                                                              l = 0, 
                                                              unit = "cm"), 
                                              hjust = 0),
            axis.text.x        = element_blank(),
            axis.title         = element_blank(),
            panel.grid         = element_blank(),
            strip.text         = element_blank(),
            legend.position    = "none") +
      add_plot_margin
    
    # Building empty plot with non-specific symptom labels and legend
    p_general_nonsp <-
      ggplot() +
      geom_bar(aes(y     = prop100, 
                   x     = symptom, 
                   fill  = reorder(stage_bin, desc(stage_bin)),
                   alpha = reorder(stage_bin, desc(stage_bin))),
               stat   = "identity", 
               width  = .5, 
               alpha  = 0,
               data   = data |> filter(symptom != "Organ-specific alarm")) 
    
    # Adjust scales, colours,  flip coordinates and add facets
    p_general_nonsp <- 
      p_general_nonsp +
      ylim(0, 102) +
      guides(fill = guide_legend(reverse = T, title = "Stage", 
                                 override.aes = list(alpha = alpha_vals,
                                                     fill = unname(unlist(colours))
                                 )
      )
      ) +
      coord_flip() +
      facet_grid2(. ~ cancer_site_desc)
    
    # Adjust theme
    p_general_nonsp <-
      p_general_nonsp +
      theme_minimal() +
      theme(axis.text.y        = element_text(size   = 14, 
                                              colour = "black", 
                                              margin = margin(t = 0, 
                                                              r = 1, 
                                                              b = 0,
                                                              l = 0, 
                                                              unit = "cm"), 
                                              hjust = 0),
            axis.text.x        = element_blank(),
            axis.title         = element_blank(),
            panel.grid         = element_blank(),
            strip.text         = element_blank(),
            plot.margin        = margin(t = .1, r = .3, b = .1, l = .3, unit = "cm"),
            legend.position    = "bottom",
            legend.box         = "horizontal")
    
    return(list(rf = p_general_rf, nonsp = p_general_nonsp))
  }
  
  # Plot blank canvas for labels and legend
  p_general <- plot_blank(nonspecific)
  
  # plot stage distribution of patients with cancer for red flag symptoms and corresponding cancer sites
  p_redflags <- lapply(cancer_sites, FUN = function(x) plot_redflag(df, x))
  
  # Name list of plots similar to list of red flag symptoms for easy access
  names(p_redflags) <- names(cancer_redflag)
  
  # plot stage distribution of patients with cancer for red flag symptoms and corresponding cancer sites
  p_nonspecific <- lapply(cancer_sites, FUN = function(x) plot_nonspecific(df, x))
  
  # Name list of plots similar to list of red flag symptoms for easy access
  names(p_nonspecific) <- names(cancer_redflag)
  
  list(
    "redflags"    = p_redflags,
    "nonspecific" = p_nonspecific,
    "general"     = p_general
  )
}

## Blood test visuals ----
plot_stage_sx_bt <- function(df, filter_by_sex, colours, lancet = F) {
  
  get_sx_order <- function(df) {
    # Get order of symptoms based on decreasing stage I-III in patients without blood tests
    df |> 
      filter(grepl("[(a) ]?All cancers$", Analysis),
             blood_test == 0, 
             stage_bin  == "1-3") |> 
      arrange(sex_lb, blood_test, prop100)  |>
      mutate(sx_short = as.factor(sx_short)) |>
      mutate(symptom_short = as.factor(symptom_short)) |>
      pull(sx_short)
  }
  
  add_labels <- function(df, lancet = F) {
    df |>
      mutate(
        p_label       = format_p_val(p, lancet),
        label         = paste0("N=", N_tumour),
        sx_short      = factor(sx_short, levels = sep_sx_order),
        bt_lb         = ifelse(blood_test == 0, "No BT", "BT"),
        Analysis_lb   = str_glue("{Analysis} ({sex_lb})")
      ) |>
      mutate(sx = as.numeric(sx_short)) |>
      arrange(sx_short)
  }
  
  plot_stage_dist <- function(p, data, shift) {
    p +
      geom_bar(aes(y = prop100, 
                   x = sx + shift, 
                   fill = reorder(stage_bin, desc(stage_bin))), 
               stat  = "identity", 
               width = .3,
               data  = data) +
      geom_text(aes(y  = -13, x = sx + shift, label = bt_lb),
                size   = 3, 
                hjust  = 0, 
                colour = "black", 
                data   = data |> filter(grepl("[(a) ]?All cancers$", Analysis))) +
      geom_text(aes(y  = 0, x = sx + shift, label = label),
                size   = 9/.pt, 
                hjust  = 0, 
                colour = "white", 
                data   = data) 
  }
  
  add_p_vals <- function(p, data, shift) {
    p  +
      geom_text(aes(y  = 101, x = sx + shift, label = p_label),
                size   = 3, 
                hjust  = 0, 
                colour = "black", 
                data   = data)
  }
  
  plot_error_bars <- function(p, data, shift){
    p +
      geom_errorbar(aes(ymin = prop100_lb,
                        ymax = prop100_ub,
                        x    = sx + shift),
                    width = .3/4,
                    data = data |> filter(stage_bin == "1-3")) +
      geom_errorbar(aes(ymin = prop100_lb,
                        ymax = prop100_ub,
                        x    = sx + shift),
                    width = .3/4,
                    data = data |> filter(stage_bin == "1-3")) 
  }
  
  df <- df |> filter(sex_lb == filter_by_sex)
  df <- 
    df |> 
    mutate(p_label = format_p_val(p, lancet)) |>
    mutate(sx_short = paste0(symptom_short, "\n(", p_label, ")"))
  sep_sx_order <- df |> get_sx_order()
  df <-  df |> add_labels()
  
  # Amount to shift bar plot by and corresponding elements
  xlab_shift <- .4
  
  p <- ggplot()
  
  # Plot stage distribution for patients without a blood test (shifted upwards)
  p <- p |> plot_stage_dist(df |> filter(blood_test == 0), shift = xlab_shift)
  
  # Plot stage distriution for patients with a blood tes
  p <- p |> plot_stage_dist(df |> filter(blood_test == 1), shift = 0) 
  
  # Add in error bars for stage 1-3 for patients without a blood test (shifted upwards)
  p <- p |> plot_error_bars(df |> filter(blood_test == 0), shift = xlab_shift)
  
  p <- p |> plot_error_bars(df |> filter(blood_test == 1), shift = 0)
  
  p <- 
    p +
    guides(fill = guide_legend(reverse = T, title = "Stage")) +
    scale_fill_manual(values = rev(unname(unlist(colours)))) +
    scale_y_continuous(breaks       = seq(0, 100, by = 20), 
                       minor_breaks = seq(0, 100, by = 10), 
                       labels       = seq(0, 100, by = 20)) +
    scale_x_continuous(breaks = seq(min(df$sx), max(df$sx)) + .3, 
                       labels = unique(df$sx_short)) +
    coord_flip() +
    facet_wrap( ~ Analysis_lb)
  
  # Adjust theme
  p <- 
    p +
    theme_minimal() +
    theme(text               = element_text(size = 12, colour = "black"),
          title              = element_text(size = 13, colour = "black"),
          strip.text         = element_text(size = 10, colour = "black"),
          axis.text.y        = element_text(size = 12, colour = "black"),
          panel.spacing.x    = unit(.1, 'cm'),
          plot.margin        = margin(t = .1, r = .1, b = .1, l = .1, unit = "cm"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  # Add labels for axes
  p <- 
    p +
    xlab("Symptom") +
    ylab("Proportion (%)")
  
  return(p)
}

