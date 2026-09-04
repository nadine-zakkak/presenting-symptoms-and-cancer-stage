# order of symptom groups
sxgrp_order <- c("Non-specific", 
                 "Lump/mass/lymph node", 
                 "Ulceration",
                 "Lower abdominal", 
                 "Upper abdominal",
                 "Respiratory", 
                 "Urological", 
                 "Central nervous system",
                 "Musculoskeletal",
                 "Skin Lesion", 
                 "Breast Symptoms", 
                 "Female specific",
                 "Male specific",
                 "None recorded")

# Discordant symptoms by sex
discord_sx_by_sex <- list(
  "Men"  = c("Post-menopausal bleeding", 
             "Other vaginal bleeding", 
             "Vaginal discharge", 
             "Vulval mass",
             "Vulval bleeding",
             "Vulval ulceration",
             "Vaginal mass"),
  
  "Women" = c("Testicular lump",
              "Erectile dysfunction",
              "Testicular pain",
              "Penile ulceration")
)

# Exemplar Cancer sites and symptoms
cancer_redflag <- list("Lung"       = "Haemoptysis",
                       "Colon"      = c("Rectal bleeding", "Change in bowel habit"),
                       "Rectum"     = c("Rectal bleeding", "Change in bowel habit"),
                       "Oesophagus" = c("Dysphagia", "Dyspepsia"),
                       "Kidney"     = "Haematuria",
                       "Pancreas"   = "Jaundice")

nonspecific_symptoms <- c("Abdominal pain (NOS)",
                          "Nausea and/or vomiting",
                          "Weight loss")

# Number of permutations for t-test -----
n_permut <- 10000

# Order of referral groups
referral_order <- c("EP", "Non-EP", "Unknown")

# Colours to be used without missing stage -----
# Colours to be used for detailed stage
stage_colours <- list(stage1 = brewer.pal("Blues", n = 9)[4],
                      stage2 = brewer.pal("Blues", n = 9)[6],
                      stage3 = brewer.pal("Blues", n = 9)[8],
                      stage4 = "#DF6B1D")

# Colours to be used for categorised stage
stagebin_colours <- list(early    = brewer.pal("Blues", n = 9)[6], 
                         advanced = "#DF6B1D")


# Colours to be used with missing stage -----
# Colours to be used for detailed stage
stage_miss_colours <- append(stage_colours, list(missing = "#4e4e4e"))

# Colours to be used for categorised stage
stagebin_miss_colours <- append(stagebin_colours, list(missing = "#4e4e4e"))

# Colours to be used for cancer sites ----
cancer_site_colour <- c("Lung"  = brewer.pal("Set1", n = 9)[3],
                        "Breast" = brewer.pal("Set1", n = 9)[8])

rtd_colours <- list("EP"      = "#d08f8a", 
                    "Non-EP"  = "#0073D1", 
                    "Unknown" = "#5cb6ff")

